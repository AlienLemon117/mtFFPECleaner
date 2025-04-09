# mtFFPECleaner
A program to remove false-positive mutations in mitochondrial DNA from FFPE samples

## Getting Started
mtFFPECleaner is written in python
- Python 3.8+
- 必需库：pysam, pandas, scikit-learn, numpy

### Installing
```shell
git clone https://github.com/AlienLemon117/mtFFPECleaner.git
cd mtFFPECleaner
```

### Usage

#### Please prepare the following inputs
1.bam file

2.mutation file
> Include 5 column —— sample | chrom:pos | ref-base | alt-base | mutation-frequency

3.reference fasta

4.ref-info.txt

5.train.csv


#### Running
```shell
python3 ../feature_extract.py \
    -b input.bam \
    -m mutation.txt \
    -r rCRS.fa \
    -i ref-info.txt \
    -t train.csv \
    -o predictions.txt
```
#### 参数说明
|参数|说明|
|:---|:---|
|-b||
|-m||
|-r||
|-i||
|-t||
|-o||

#### Outputs
A txt file
