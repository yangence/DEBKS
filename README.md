# DEBKS: a tool to detect differentially expressed circular RNA

## Introduction

DEBKS is a convenient and user-friendly program to streamline the discovery of differentially expressed circRNA (DEC) between two RNA-seq sample groups with replicates. DEBKS includes four modules: (1) "merge" collects circRNA junction information from output file of circRNA detection software. (2) "anno" annotates circRNA based on the circRNA position. (3) "count" calcuates linear junction based on the circRNA position. (4) "dec" identifies DEC with rMATS statistical model.

## Availability

DEBKS is a free software, which can be downloaded from https://github.com/yangence/DEBKS

## Prequired Softwares and Packages

Python 3.x.x and corresponding versions of Pysam, NumPy, Pandas, and SciPy.

## Installation

Install latest release via conda
```
conda install -c colinliuzelin DEBKS -c bioconda
```

Install latest release from source codes
```
git clone https://github.com/yangence/DEBKS.git
cd DEBKS
pip install -r requirements.txt
python setup.py install
```

## Required Files:

Users can prepare the external files under the following instructions:

1) Indexed genome fasta file

```bash
samtools faidx $genome
```

2) Tabix indexed gene annotation GTF file

```bash
grep -v '#' $gtf |sort -k 1,1 -k 4,4n |bgzip >sort.gtf.gz
tabix sort.gtf.gz
```

## Usage

```bash
Usage: DEBKS (merge | count | anno | dec) [options]

Command:
    merge            Merge circRNA junction from other software
    count            Count linear junction with BAM file
    anno             Annotate circRNA with gene annotation
    dec              Detect differentially expressed circRNA with junction information
```


### Merge

```bash
Usage: DEBKS  merge (-s software | -d designate) (-n name | -f file)  [-b pos_based]  [-o output]

Example: DEBKS  merge -s CIRI2 -n s1.ciri,s2.ciri,s3.ciri

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -s software                 Name of software used to detect circRNA (e.g. ciri2, circexplorer2, find_circ).
    -d designate                Designate the column postion of chr,start,end,circ,linear in the tab format output file of detection software (e.g. 1,2,3,4,5).
                                The first four postion is required.
    -b pos_based                0- or 1-based position for designated start and end position, separated by comma [default: 1,1].
    -n name                     Names of output file of detection software, separated by comma (e.g. s1,s2,s3).
    -f file                     File includes names of output file of detection software, separated by line breaks.
    -o output                   Output prefix [default: merge].
```

"merge" can automatically collect circular and linear junction counts from output files of well-known circRNA detection software or user designated positions.
Note: "merge" adjust start and end position to 1-based.
### Count

```bash
Usage: DEBKS count -c circ_pos (-n name | -f file) [-l library_type] [-t threads]  [-a hangover_len]  [-o output]

Example: DEBKS  count -c merge_pos.txt -n s1.sort.bam,s2.sort.bam,s3.sort.bam

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -c circ_pos                 File of circRNA position in tab format (the first three columns are chr, start, end).
    -n name                     Names of sorted bam file, separated by comma (e.g. bam1,bam2,bam3).
    -f file                     File includes names of sorted bam file, separated by line breaks.
    -l library_type             0 represent unstrandard, 1 represent fr-firststrand (e.g., dUTP protocol), 2 represent fr-secondstrand [default: 0].
    -t threads                  Number of threads to detect DE circRNA [default: 4].
    -a hangover_len             The minimum hangover length for junction quantification [default: 6].
    -o output                   Output prefix [default: linear].
```

"count" calcuate the linear junction counts with indexed sorted aligment bam file produced by RNA-seq aligners, e.g., STAR, HISAT2, and Tophat2.
To compatiable with circRNA quantification, the minimum hangover length to detect circRNA can be inputed as the threshold for linear junction quantification with option '-a'.

### Annotation

```bash
Usage: DEBKS anno -c circ_pos -g gene_anno -m genome [-o output]

Example: DEBKS  anno -c merge_pos.txt -m hg19.fa -g gencode.v19.annotation.sort.gtf.gz

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -c circ_pos                 File of circRNA position in tab format (the first three columns are chr, start, end).
    -g gene_anno                Tabix indexed gtf File of gene annotation.
    -m genome                   Fasta file of genome.
    -o output                   Output prefix [default: circ].
```

"anno" annotates circRNAs with circRNA position and predicts the potential circRNA length based on gene annotation file.

### DEC

```bash
Usage: DEBKS dec -c circ (-l linear |--c2 circ2) [-p] [-n sample_num] [-f filter] [-t threads] [-d cutoff] [-r read_len] [-a hangover_len] [-e circ_len] [--e2 circ2_len] [-o output]

Example: DEBKS  dec -c merge_circ.txt -l merge_linear.txt -n 3,3 -f 12 -t 20 -e anno_len.txt

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -c circ                     Circular junction counts file in tab format, the first three columns is circRNA position (chr,start,end).
    -l linear                   Linear junction counts file, the format is same with -c.
    --c2 circ2                  Circular junction counts file in tab format, the format is same with -c.
    -p                          Sample group is paired.
    -n sample_num               Number of samples in each group, separated by comma (e.g. 2,3).
    -f filter                   Required total circular juction counts in all samples to filter out low expressed circRNAs [default: 0].
    -t threads                  Number of threads to detect DE circRNA [default: 4].
    -d cutoff                   Cutoff of DE circRNA in two groups [default: 0.05].
    -r read_len                 Length of RNA-seq reads [default: 150].
    -a hangover_len             The minimum hangover length for junction quantification [default: 6].
    -e circ_len                 File of circRNA length in tab format. The first four columns is circRNA position (chr, start, end) and circRNA length.
    --e2 circ2_len              File of circRNA2 length in tab format, only work with --c2. The format is same with -e.
    -o output                   Output file [default: dec_circRNA.txt].
```

"dec" detects DEC with circular and linear junction file which can be collected from upstream analysis. It also can detect DEC in alternative splicing level of circRNA. Options '--c2' supports input of circular junction file with same coordinate position in '-c' but different counts.
Option '-e' and '--e2' accepts circRNA length to adjust the calucation of junction ratio. Users can peform 'anno' command or perform circRNA full length analysis software, e.g., CIRI-full and circAST to get the potential length.
Note: Option '--c2' and '-l' is mutually exclusive.
#### Output files
|Field|Description|
|:------:|:------|
|chr| chromosome of circRNA|
|start| coordinate of start back-splicied site|
|end| coordinate of end back-splicied site|
|cjc_1| counts of circular junction of group1|
|cjc_2| counts of circular junction of group1|
|ljc_1 or cjc2_1| counts of linear or circular2 junction of group1|
|ljc_2 or cjc2_2| counts of linear or circular2 junction of group2|
|adj_cjc_len| length adjust for cjc|
|adj_ljc_len| length adjust for ljc|
|pbsi1| circular junction ratio of group1|
|pbsi2| circular junction ratio of group2|
|delta_pbsi| average of the difference between pbsi1 and pbsi2|
|P| the significance of delta_pbsi with user defined threshold|
|FDR| Benjamini-Hochberg corrected FDR of the above P|

## Copyright and License Information

Copyright (C) 2020 Zelin Liu (zlliu@bjmu.edu.cn). See the [LICENSE](https://github.com/yangence/DEBKS/blob/master/LICENSE) file for license rights and limitations.
