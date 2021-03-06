# Example

## Required Files:

Users can prepare the external files under the following instructions:

1) Indexed genome fasta file

```bash
samtools faidx hg19.fa
```

2) Tabix indexed gene annotation GTF file

```bash
grep -v '#' gencode.v19.annotation.gtf |sort -k 1,1 -k 4,4n |bgzip >gencode.v19.annotation.sort.gtf.gz
tabix gencode.v19.annotation.sort.gtf.gz
```

3) Simulation of expression level of circRNA and linear RNA

```bash
# RNA-seq dataset1
mkdir dataset1
python sim_circ.py 0.01 dataset1 gencode.v19.annotation.sort.gtf.gz outJunction_mouse.txt
# RNA-seq dataset2
mkdir dataset2
python sim_circ2.py 0.01 dataset2 gencode.v19.annotation.sort.gtf.gz outJunction_mouse.txt
```

4) Simulation of RNA-seq

```bash
# RNA-seq dataset1
cd dataset1
for i in {1,2};do for j in {1,2,3};do python ../sim_seq.py coverage_0.01_g${i}_rep${j}.txt g${i}_rep${j} gencode.v19.annotation.sort.gtf.gz hg19.fa;done;done
# RNA-seq dataset2
cd dataset2
for i in {1,2};do for j in {1,2,3};do python ../sim_seq2.py coverage_0.01_g${i}_rep${j}.txt g${i}_rep${j} gencode.v19.annotation.sort.gtf.gz hg19.fa;done;done
```

5) Detect circRNA

In this example, we use CIRI2 to detect circRNA.

```bash
cd dataset1
mkdir CIRI
for i in {1,2};
do
for j in {1,2,3};
do
bwa mem -t 70 -T 19 hg19.fa g${i}_rep${j}_1.fq.gz g${i}_rep${j}_2.fq.gz >CIRI/g${i}_rep${j}.sam
CIRI2.pl -I  CIRI/g${i}_rep${j}.sam -O CIRI/g${i}_rep${j}.ciri -F hg19.fa -A gencode.v19.annotation.gtf -0  -T 20;
done;
done

```

6) Alignment (Optional)

In this example, users can choose other RNA-seq aligners not only STAR, such as Tophat and HISAT2.
Users should prepare STAR index first.

```bash
mkdir bam
for i in {1,2};
do 
for j in {1,2,3};
do
mkdir -p bam/g${i}_rep${j}/;
STAR --chimSegmentMin 6 --chimOutType WithinBAM  \
	--runThreadN 20 --outSAMtype BAM SortedByCoordinate \
	--genomeDir ${STAR_index} --chimJunctionOverhangMin 6 \
	--outSJfilterOverhangMin 6 6 6 6 \
	--readFilesCommand zcat \
	--outFileNamePrefix bam/g${i}_rep${j}/ \
	--readFilesIn g${i}_rep${j}_1.fq.gz g${i}_rep${j}_2.fq.gz ;
samtools index bam/g${i}_rep${j}/Aligned.sortedByCoord.out.bam;
done;
done
```

## Merge

```bash
ls CIRI/*.ciri>samplelist.txt
DEBKS merge  -s ciri2 -f samplelist.txt
```

## Count (optional)

Users can skip this step if the linear junction is available in upstream analysis.

```bash
ls bam/*/Aligned.sortedByCoord.out.bam>bamlist.txt
DEBKS count -c merge_pos.txt -f bamlist.txt -t 6
```

## Annotation (optional)

```bash
DEBKS anno -c merge_pos.txt  -m hg19.fa  -g gencode.v19.annotation.sort.gtf.gz
```

## DEC

```bash
DEBKS dec -c merge_circ.txt -l merge_linear.txt -t 20 -e anno_len.txt
DEBKS dec -c merge_circ.txt -l linear_jun.txt -t 20 -e anno_len.txt -o dec_star.txt #optional
```