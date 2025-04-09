from importlib.resources import files
import pysam
from collections import defaultdict
import os
import pandas as pd
from tqdm import tqdm  # 用于显示进度条（可选）
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.impute import SimpleImputer
import numpy as np
import argparse


'''usage:
python main.py \
  -b input.bam \
  -m mutations.txt \
  -r hg38.fa \
  -i ref-info.txt \
  -t train_data.csv \
  -o predictions.txt
  '''


'''
input-mutation file
Sample1 chr1:12345 A T 0.3
Sample2 chr2:67890 C G 0.8
'''

def get_filepath(file_name):
    try:
        return files("your_package").joinpath(file_name).read_text(encoding="utf-8")
    except FileNotFoundError:
        # 降级到 pkg_resources
        from pkg_resources import resource_string
        return resource_string("your_package", file_name).decode("utf-8")


def calculate_sob(bam_path, txt_path, ref_fasta):
    # 数据结构: {样本: {位点: {'ALT_F1R2': int, 'ALT_F2R1': int}}}
    sample_data = defaultdict(lambda: defaultdict(lambda: {
        'ALT_F1R2': 0,
        'ALT_F2R1': 0,
        'ref': None,
        'alt': None
    }))

    # 读取突变文件，按样本分组
    with open(txt_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            sample = parts[0]
            chrom, pos = parts[1].split(':')
            ref = parts[2].upper()
            alt = parts[3].upper()
            pos = int(pos)
            sample_data[sample][(chrom, pos)]['ref'] = ref
            sample_data[sample][(chrom, pos)]['alt'] = alt

    # 读取参考基因组
    ref = pysam.FastaFile(ref_fasta)

    # 处理BAM文件
    bam = pysam.AlignmentFile(bam_path, 'rb')

    # 遍历每个位点
    for pileupcolumn in bam.pileup(stepper='nofilter'):
        chrom = pileupcolumn.reference_name
        pos_pileup = pileupcolumn.reference_pos + 1  # 1-based

        # 获取参考碱基
        current_ref = ref.fetch(chrom, pos_pileup - 1, pos_pileup).upper()

        # 遍历所有样本，检查该位点是否为该样本的突变位点
        for sample in sample_data:
            if (chrom, pos_pileup) not in sample_data[sample]:
                continue

            expected_ref = sample_data[sample][(chrom, pos_pileup)]['ref']
            expected_alt = sample_data[sample][(chrom, pos_pileup)]['alt']

            # 验证参考碱基
            if current_ref != expected_ref:
                print(
                    f"警告：样本 {sample} 的位点 {chrom}:{pos_pileup} 参考碱基不匹配（实际：{current_ref}，预期：{expected_ref}）")
                continue

            # 统计ALT计数
            for read in pileupcolumn.pileups:
                if read.is_del or read.is_refskip:
                    continue

                # 获取read碱基（考虑反向互补）
                qpos = read.query_position
                if qpos is None:
                    continue
                read_seq = read.alignment.query_sequence
                is_reverse = read.alignment.is_reverse
                base = read_seq[qpos].upper()
                if is_reverse:
                    base = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}.get(base, 'N')

                # 仅统计支持ALT的read
                if base != expected_alt:
                    continue

                # 判断F1R2/F2R1
                is_read1 = read.alignment.is_read1
                if (is_read1 and not is_reverse) or (not is_read1 and is_reverse):
                    sample_data[sample][(chrom, pos_pileup)]['ALT_F1R2'] += 1
                else:
                    sample_data[sample][(chrom, pos_pileup)]['ALT_F2R1'] += 1

    # 计算SOB值
    sob_results = defaultdict(dict)
    for sample in sample_data:
        for (chrom, pos), counts in sample_data[sample].items():
            f1r2 = counts['ALT_F1R2']
            f2r1 = counts['ALT_F2R1']
            total = f1r2 + f2r1
            sob = abs(f1r2 - f2r1) / total if total > 0 else 0.0
            sob_results[sample][f"{chrom}:{pos}"] = round(sob, 4)

    return dict(sob_results)
'''
{
    "Sample1": {
        "chr1:12345": 0.75,
        "chr2:67890": 0.6
    },
    "Sample2": {
        "chr3:44444": 0.33
    }
}
'''
## 调用示例
## sob_dict = calculate_sob("input.bam", "mutations.txt", "hg38.fa")
## print(sob_dict)

def load_reference_data(ref_path):
    """加载参考信息文件到字典(ref-info.txt)"""
    ref_dict = {}
    with open(ref_path, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            # 生成染色体:位点的组合键（例如 "chrM:73"）
            chrom_pos = f"{cols[0]}:{cols[1]}"
            ref_dict[chrom_pos] = {
                'ref-context': cols[2],
                '5-upstream-base': cols[3],
                '3-downstream-base': cols[4],
                '%-of-G-within-50bp': float(cols[5]),
                '%-of-C-within-50bp': float(cols[6]),
                '%-of-G-within-100bp': float(cols[7]),
                '%-of-C-within-100bp': float(cols[8])
            }
    return ref_dict


def calculate_features(mutation_file, ref_dict):
    """主处理函数"""
    # 初始化数据结构
    samples_dict = {}

    # 第一次遍历：收集所有样本的统计信息
    with open(mutation_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5: continue  # 跳过不完整行

            sample, chrom_pos, ref, alt, vaf = parts[0], parts[1], parts[2], parts[3], float(parts[4])

            # 初始化样本数据结构
            if sample not in samples_dict:
                samples_dict[sample] = {
                    'mutations': {},
                    'stats': {
                        'total_muts': 0,
                        'C>T_muts': 0, 'G>A_muts': 0,
                        'vaf_1_90_total': 0,
                        'vaf_1_90_C>T': 0, 'vaf_1_90_G>A': 0
                    }
                }

            # 统计突变类型
            is_C_T = (ref.upper() == 'C' and alt.upper() == 'T')
            is_G_A = (ref.upper() == 'G' and alt.upper() == 'A')

            samples_dict[sample]['stats']['total_muts'] += 1
            if is_C_T: samples_dict[sample]['stats']['C>T_muts'] += 1
            if is_G_A: samples_dict[sample]['stats']['G>A_muts'] += 1

            # 统计VAF 1-90%的突变
            if 0.01 <= vaf <= 0.9:
                samples_dict[sample]['stats']['vaf_1_90_total'] += 1
                if is_C_T: samples_dict[sample]['stats']['vaf_1_90_C>T'] += 1
                if is_G_A: samples_dict[sample]['stats']['vaf_1_90_G>A'] += 1

            # 存储每个突变位的原始数据
            samples_dict[sample]['mutations'][chrom_pos] = {
                'VAF': vaf,
                'ref': ref,
                'alt': alt
            }

    # 第二次遍历：计算特征并整合参考数据
    result_dict = {}
    for sample, data in samples_dict.items():
        sample_features = {}
        stats = data['stats']

        # 计算样本级特征
        total = stats['total_muts']
        vaf_total = stats['vaf_1_90_total']

        def safe_divide(a, b):
            return (a / b ) if b > 0 else 0.0

        sample_level = {
            '%-of-C>T&G>A-muts': safe_divide(
                (stats['C>T_muts'] + stats['G>A_muts']), total),
            '%-of-C>T-muts': safe_divide(stats['C>T_muts'], total),
            '%-of-G>A-muts': safe_divide(stats['G>A_muts'], total),
            '%-of-muts(1%-90%)': safe_divide(vaf_total, total),
            '%-of-muts(1%-90%-C>T&G>A)': safe_divide(
                (stats['vaf_1_90_C>T'] + stats['vaf_1_90_G>A']), vaf_total),
            '%-of-muts(1%-90%-C>T)': safe_divide(stats['vaf_1_90_C>T'], vaf_total),
            '%-of-muts(1%-90%-G>A)': safe_divide(stats['vaf_1_90_G>A'], vaf_total)
        }

        # 处理每个突变位点
        for chrom_pos, mut_data in data['mutations'].items():
            # 获取参考数据
            ref_data = ref_dict.get(chrom_pos, {})

            # 合并所有特征
            feature_dict = {
            **sample_level,  # 样本级特征
            'VAF': mut_data['VAF'],
            'ref-context': ref_data.get('ref-context', 'N/A'),
            'ref-base':mut_data['ref'],
            'alt-base':mut_data['alt'],
            '5-upstream-base': ref_data.get('5-upstream-base', 'N/A'),
            '3-downstream-base': ref_data.get('3-downstream-base', 'N/A'),
            '%-of-G-within-50bp': ref_data.get('%-of-G-within-50bp', 0.0),
            '%-of-C-within-50bp': ref_data.get('%-of-C-within-50bp', 0.0),
            '%-of-G-within-100bp': ref_data.get('%-of-G-within-100bp', 0.0),
            '%-of-C-within-100bp': ref_data.get('%-of-C-within-100bp', 0.0)
            }

            if chrom_pos not in sample_features:
                sample_features[chrom_pos] = {}
            sample_features[chrom_pos].update(feature_dict)

        result_dict[sample] = sample_features

    return result_dict

'''# 输出结构示例
{
    "Sample1": {
        "chrM:12345": {
            "%-of-C>T&G>A-muts": 25.0,
            "VAF": 0.3,
            "ref-context": "ATG",
            "%-of-G-within-50bp": 12.5,
            # ...其他特征...
        }
    }
}'''


def merge_features_to_dataframe(feature_dict, sob_dict):
    """
    将特征字典与SOB字典合并为DataFrame

    参数:
    feature_dict -- 包含突变特征的主字典
    sob_dict     -- 包含SOB值的补充字典

    返回:
    pandas.DataFrame -- 合并后的数据框，每行代表一个突变
    """
    rows = []

    # 遍历所有样本和突变位点
    for sample, mutations in tqdm(feature_dict.items()):
        # 获取该样本在SOB字典中的数据（若无则为空字典）
        sample_sob = sob_dict.get(sample, {})

        for chrom_pos, features in mutations.items():
            # 生成唯一标识符
            unique_id = f"{sample}_{chrom_pos.replace(':', '-')}"

            # 创建行数据副本避免修改原数据
            row_data = features.copy()

            # 添加标识列
            row_data.update({
                'unique_id': unique_id,
                'sample': sample,
                'chrom_pos': chrom_pos
            })

            # 添加SOB特征（若无则设为NaN）
            sob_value = sample_sob.get(chrom_pos, float('nan'))
            row_data['SOB'] = sob_value

            rows.append(row_data)

    # 创建DataFrame
    df = pd.DataFrame(rows)

    # 调整列顺序（将标识列放在前面）
    columns = ['unique_id', 'sample', 'chrom_pos', 'SOB'] + \
              [col for col in df.columns if col not in
               ['unique_id', 'sample', 'chrom_pos', 'SOB']]

    return df[columns]

'''
df = merge_features_to_dataframe(feature_dict, sob_dict)
   unique_id   sample   chrom_pos  SOB  %-of-C>T&G>A-muts  VAF ref-context  %-of-G-within-50bp
0  Sample1_chrM-12345  Sample1  chrM:12345  0.75               25.0  0.3         ATG                 12.5
'''
import chardet

def ml_pipeline_with_dummies(train_csv, predict_csv, output_txt, target_col='label'):
    """
    支持哑变量转换的机器学习流程

    参数:
    train_csv   -- 训练集CSV文件路径
    predict_csv -- 需要预测的数据集CSV文件路径
    output_txt  -- 预测结果输出路径
    target_col  -- 目标变量列名 (默认'label')

    usage:ml_pipeline_with_dummies(
        train_csv="train_data.csv",
        predict_csv="predict_data.csv",
        output_txt="predictions_with_dummies.txt",
        target_col='is_somatic'
    )
    """
    with open(train_csv, 'rb') as f:
        result = chardet.detect(f.read())
        print(f"检测到编码: {result['encoding']}")
    train_df = pd.read_csv(train_csv, encoding=result['encoding'])

    with open(predict_csv, 'rb') as f:
        result = chardet.detect(f.read())
        print(f"检测到编码: {result['encoding']}")
    
    # 读取数据
    #train_df = pd.read_csv(train_csv, encoding='')  # 或 'gbk', 'latin1'
    predict_df = pd.read_csv(predict_csv, encoding=result['encoding'])

    # 标识列定义
    identifier_cols = ['unique_id', 'sample', 'chrom_pos']

    # 定义需要显式转换为哑变量的分类列（示例）
    explicit_cat_cols = ['ref-context', 'ref-base', 'alt-base','3-downstream-base', '5-upstream-base']

    # 自动识别其他分类特征
    remaining_cat_cols = train_df.drop(identifier_cols + [target_col], axis=1) \
        .select_dtypes(include=['object', 'category']).columns.tolist()
    all_cat_cols = list(set(explicit_cat_cols + remaining_cat_cols))

    # 数值特征（排除分类特征）
    numerical_cols = train_df.drop(identifier_cols + [target_col] + all_cat_cols, axis=1) \
        .select_dtypes(include=np.number).columns.tolist()

    # 构建哑变量转换管道
    cat_transformer = Pipeline(steps=[
        ('imputer', SimpleImputer(strategy='most_frequent')),
        ('onehot', OneHotEncoder(
            drop='if_binary',  # 对二分类变量自动删除一列
            handle_unknown='ignore',
            sparse_output=False  # 输出密集矩阵
        ))
    ])

    # 完整预处理流程
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', SimpleImputer(strategy='median'), numerical_cols),
            ('cat', cat_transformer, all_cat_cols)
        ],
        remainder='drop'  # 丢弃未定义处理的列
    )

    # 构建模型管道
    pipeline = Pipeline(steps=[
        ('preprocessor', preprocessor),
        ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
    ])

    # 模型训练
    #print(f"训练模型... 使用分类特征: {all_cat_cols}")
    X_train = train_df.drop(identifier_cols + [target_col], axis=1)
    y_train = train_df[target_col]
    pipeline.fit(X_train, y_train)

    # 预测
    #print("生成预测...")
    X_pred = predict_df.drop(identifier_cols, axis=1)
    predictions = pipeline.predict(X_pred)
    proba = pipeline.predict_proba(X_pred)[:, 1]  # 假设二分类

    # 构造结果数据框
    result_df = predict_df[identifier_cols].copy()
    result_df['prediction'] = predictions
    result_df['probability'] = proba

    # 保存结果
    result_df.to_csv(output_txt, sep='\t', index=False)
    print(f"预测结果已保存至 {output_txt}")

    return result_df


def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='突变预测流程')
    parser.add_argument('-b', '--bam', required=True, help='输入BAM文件路径')
    parser.add_argument('-m', '--mutation', required=True, help='突变文件路径')
    parser.add_argument('-r', '--ref', required=True, help='参考基因组文件路径')
    parser.add_argument('-i', '--ref_info', required=True, help='参考信息文件路径')
    parser.add_argument('-t', '--train', required=True, help='训练集CSV文件路径')
    parser.add_argument('-o', '--output', required=True, help='预测结果输出路径')
    args = parser.parse_args()

    # Step 1: 计算SOB值
    print("正在计算SOB值...")
    sob_dict = calculate_sob(args.bam, args.mutation, args.ref)

    # Step 2: 加载参考数据
    print("正在加载参考数据...")
    ref_dict = load_reference_data(args.ref_info)

    # Step 3: 计算特征
    print("正在计算特征...")
    feature_dict = calculate_features(args.mutation, ref_dict)

    # Step 4: 合并特征到DataFrame
    print("正在合并特征数据...")
    df = merge_features_to_dataframe(feature_dict, sob_dict)

    # 保存临时特征文件（可选）
    temp_feature_csv = "temp_features.csv"
    df.to_csv(temp_feature_csv, index=False)

    # Step 5: 机器学习预测
    print("正在进行机器学习预测...")
    ml_pipeline_with_dummies(
        train_csv=args.train,
        predict_csv=temp_feature_csv,
        output_txt=args.output,
        target_col='label'  # 根据实际训练集的列名修改
    )
    print("流程完成！")


if __name__ == "__main__":
    main()
