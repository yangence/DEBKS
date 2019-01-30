# DEBKS: a tool to detect differentially expressed circular RNA

## Introduction

DEBKS is a convenient and user-friendly program to streamline the discovery of differentially expressed circRNA between two RNA-seq sample groups with replicates. DEBKS combines well-known software *CIRCexplorer2*  for circRNA detection and annotation in chimeric RNA-seq reads, with rMATS statistical model for identifying differential isoform ratios using RNA-seq sequence count data with replicates.

## Availability

DEBKS is a free software, which can be downloaded from https://github.com/yecLab/DEBKS

## Prequired Softwares and Packages

1. Python 3.x.x and corresponding versions of NumPy, Pandas, and SciPy.

2. [STAR 2.6.1](https://github.com/alexdobin/STAR)

3. [gtfToGenePred](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred)

4. [SAMtools 1.9](https://github.com/samtools/samtools/releases/tag/1.9) 

5. [CIRCexplorer2](https://circexplorer2.readthedocs.io/en/latest/tutorial/setup/#installation-and-setup)

## Required Files:

Users can prepare the external files under the following instructions:

1) Genome fasta file

2) Gene annotation GTF file

## Usage

### Raw Fastq File Provided

If raw fastq file provided, DEBKS can map reads and calcuate differential PBSI in the following command:

```bash
python DEBKS_main.py -g genomeFasta -s1 s1File -s2 s2File  \
    -STARindex STARIndexDir -gtf gtfFile -o outDir \
    -read readType -len readLength [options]*
```

#### Required Files:

1) s1File contains sample_1 fastq files:

```
  Path/sample_1.Rep1.R1.fastq.gz;Path/sample_1.Rep1.R2.fastq.gz
  Path/sample_1.Rep2.R1.fastq.gz;Path/sample_1.Rep2.R2.fastq.gz
  Path/sample_1.Rep3.R1.fastq.gz;Path/sample_1.Rep3.R2.fastq.gz
```

2) s2File contains sample_2 fastq files:

```
  Path/sample_2.Rep1.R1.fastq.gz;Path/sample_2.Rep1.R2.fastq.gz
  Path/sample_2.Rep2.R1.fastq.gz;Path/sample_2.Rep2.R2.fastq.gz
  Path/sample_2.Rep3.R1.fastq.gz;Path/sample_2.Rep3.R2.fastq.gz
```

3) Genome index built by STAR

```bash
STAR --runMode genomeGenerate --runThreadN threads \
  --genomeFastaFiles genomeFasta \
  --sjdbGTFfile gtfFile \
  --sjdbOverhang readLength-1 \
  --genomeDir STARIndexDir
```

#### Example

```bash
python DEBKS_main.py -g hg19.fa -s1 sample_1.txt -s2 sample_2.txt -STARindex hg19_STAR/ \
    -gtf gencode.v19.annotation.gtf -o out_test -t 40 -read pair -len 150 -c 0.1 -a 6
```

### STAR Alignment Results Provided

In this mode, users can map RNA-seq reads by themself with the following command:

```bash
STAR ---genomeDir STARIndexDir -chimSegmentMin anchorLength \
    --runThreadN threads --outSAMtype BAM Unsorted --alignSJDBoverhangMin anchorLength \
    --alignSJoverhangMin anchorLength	--chimJunctionOverhangMin anchorLength \
    --outSJfilterOverhangMin -1 anchorLength -1 -1
```

Then, users can employ DEBKS to calcuate differential PBSI in the following command:

```bash
python DEBKS_main.py -g genomeFasta \
    -s1CJ s1CJFile -s2CJ s2CJFile -s1SJ s1SJFile -s2SJ s2SJFile \
    -gtf gtfFile -o outDir -read readType -len readLength [options]*
```

#### Required Files:

1) File contains sample_1 chimeric junction files from STAR output:

```
  Path/sample_1.Rep1.Chimeric.out.junction
  Path/sample_1.Rep2.Chimeric.out.junction
  Path/sample_1.Rep3.Chimeric.out.junction
```

2) File contains sample_2 chimeric junction files from STAR output:

```
  Path/sample_2.Rep1.Chimeric.out.junction
  Path/sample_2.Rep2.Chimeric.out.junction
  Path/sample_2.Rep3.Chimeric.out.junction
``` 

3) File contains sample_1 splicing junction files from STAR output:

```
  Path/sample_1.Rep1.SJ.out.tab
  Path/sample_1.Rep2.SJ.out.tab
  Path/sample_1.Rep3.SJ.out.tab
```

4) File contains sample_2 splicing junction files from STAR output:

```
  Path/sample_2.Rep1.SJ.out.tab
  Path/sample_2.Rep2.SJ.out.tab
  Path/sample_2.Rep3.SJ.out.tab
```

#### Example

```bash
python DEBKS_main.py -g genomeFasta -s1CJ sample_1.CJ.txt -s2CJ sample_2.CJ.txt -s1SJ sample_1.SJ.txt \
   -S2SJ sample_2.SJ.txt -gtf gencode.v19.annotation.gtf -o out_test  -t 40 -read pair -len 150 -c 0.1 -a 6
```

### Required Parameters:
	-g          <str>       Genome Fasta file

	-STARindex  <str>       STAR alignment index directory
	
	-s1         <str>       FASTQ files of sample 1, replicates in different lines, paired files are separated by semicolon
	
	-s2         <str>       FASTQ files of sample 2, replicates in different lines, paired files are separated by semicolon

	-s1CJ       <str>       Chimeric junction of sample 1 group, replicates in different lines

	-s2CJ       <str>       Chimeric junction of sample 2 group, replicates in different lines

	-s1SJ       <str>       Spliced junction of sample 1 group, replicates in different lines

	-s2SJ       <str>       Spliced junction of sample 2 group, replicates in different lines

	-gtf        <str>       GTF file

	-o          <str>       Output directory of the result files

	-read       <str>       RNA-seq reads are single- or pair-end reads. [single, pair]

	-len        <int>       Read length of RNA-seq reads

Note: parameters -STARindex -s1 -s2 are mutually exclusive with -s1CJ -s2CJ -s1SJ -s2SJ
### Optional Parameters:
	-h, --help              Show this help message and exit
	
	-t          <int>       Number of processors [1]

	-c          <float>     Required PBSI difference cutoff between the two samples [0.01]

	-a          <int>       Anchor length for counting chimeric or splicing junctions [6]

	-KeepTemp               Keep the temporary files. Disable by default.

### DEBKS Results Summary

|Field|Description|
|:------:|:------|
|chr| chromosome of circRNA|
|start| coordinate of start back-splicied site (0-based)|
|end| coordinate of end back-splicied site (1-based)|
|strand| '+' or '-'|
|exonCount| exon number of circRNA with comma-delimiter|
|exonSizes| length for each exon with comma-delimiter|
|exonOffsets| offset for each exon with comma-delimiter|
|geneID| ID of gene|
|isoformID| ID of isoform|
|flankIntron| flanking intron of back-spliced sites|
|linearExonL| coordinate of start site for left flanking exon|
|linearExonR| coordinate of end site for right flanking exon|
|SJL1| counts of spliced junction of left flanking in sample 1 group with comma-delimiter|
|SJL1| counts of spliced junction of left flanking in sample 2 group with comma-delimiter|
|SJR1| counts of spliced junction of right flanking in sample 1 group with comma-delimiter|
|SJR1| counts of spliced junction of right flanking in sample 2 group with comma-delimiter|
|inc1| counts of inclusion spliced junction of sample 1 group with comma-delimiter, equal to sum of SJL1 and SJR1|
|inc2| counts of inclusion spliced junction of sample 2 group with comma-delimiter, equal to sum of SJL2 and SJR2|
|bs1| counts of back-spliced junction in sample 1 group|
|bs2| counts of back-spliced junction in sample 2 group|
|effective_inclusion_length| length adjust for inc1 and inc2|
|effective_bs_length| length adjust for bs1 and bs2|
|PBSI1| percent back-spliced in of sample 1|
|PBSI2| percent back-spliced in of sample 2|
|P| the significance of differential PBSI with user defined threshold|
|FDR| Benjamini-Hochberg corrected FDR of the above P|

## Citation



## Copyright and License Information

Copyright (C) 2019 Zelin Liu (zlliu@bjmu.edu.cn). See the [LICENSE](https://github.com/haersliu/test/blob/master/LICENSE) file for license rights and limitations.
