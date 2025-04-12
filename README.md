# mtFFPECleaner
A program to remove false-positive mutations in mitochondrial DNA from FFPE samples 

## Getting Started
mtFFPECleaner is written in python
- Python 3.8+
- Dependencies：pysam, pandas, scikit-learn, numpy

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
#### Parameter description
|Parameter|Note|
|:---|:---|
|-b|input bam file|
|-m|input mutation file (txt format, include 5 column)|
|-r|ref fasta|
|-i|necessary feature input file (in the bin folder) |
|-t|train dataset (in the bin folder)|
|-o|output file|

#### Outputs
A txt file includs 5 column ：unique_mutation_id | sample | predict_label | predict_probability

> 📢 **Update Notice**  
> This documentation is **continuously updated**. Check back frequently for the latest changes!
