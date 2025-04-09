# mtFFPECleaner
A program to remove false-positive mutations in mitochondrial DNA from FFPE samples

## Getting Started
mtFFPECleaner is written in python

### Installing
```shell

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

#### Outputs
A txt file
