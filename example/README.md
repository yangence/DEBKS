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
python sim_circ.py genome.fa 0.01 dataset1 gencode.v19.annotation.sort.gtf.gz outJunction_mouse.txt
# RNA-seq dataset2
mkdir dataset2
python sim_circ2.py genome.fa 0.01 dataset2 gencode.v19.annotation.sort.gtf.gz outJunction_mouse.txt

```

4) Simulation of RNA-seq

```bash
# RNA-seq dataset1
cd dataset1
for i in {1,2};do for j in {1,2,3};do python sim_seq.py coverage_0.01_g${i}_rep${j}.txt g${i}_rep${j} ;done;done
# RNA-seq dataset2
cd dataset2
for i in {1,2};do for j in {1,2,3};do python sim_seq.py coverage_0.01_g${i}_rep${j}.txt g${i}_rep${j} ;done;done
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
		bwa mem -t 70 -T 19 hg19.fa g${i}_rep${j}_1.fq g${i}_rep${j}_2.fq >CIRI/g${i}_rep${j}.sam
		CIRI2.pl -I  CIRI/g${i}_rep${j}.sam -O CIRI/g${i}_rep${j}.ciri -F hg19.fa -A gencode.v19.annotation.gtf -0  -T 20;
	done;
done

```

6) Alignment (Optional)

In this example, users can choose other RNA-seq aligners not only STAR, such as Tophat and hisat2.
Users should prepare STAR index first.

```bash
mkdir bam
for i in {1,2};
do 
	for j in {1,2,3};
	do 
		STAR --genomeLoad LoadAndKeep --chimSegmentMin 6 --chimOutType WithinBAM  \
			--runThreadN 20 --outSAMtype BAM Unsorted \
			--genomeDir STAR_index --chimJunctionOverhangMin 6 \
			--outSJfilterOverhangMin 6 6 6 6 \
			--outFileNamePrefix bam/g${i}_rep${j}/ 
			--readFilesIn g${i}_rep${j}_1.fq g${i}_rep${j}_2.fq ;
		samtools sort -m 10G -@ 4 bam/g${i}_rep${j}/Aligned.out.bam -o bam/g${i}_rep${j}/Aligned.out.sort.bam;
		samtools index bam/g${i}_rep${j}/Aligned.out.sort.bam;
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
ls bam/*/Aligned.out.sort.bam>bamlist.txt
DEBKS count -c merge_pos.txt -f bamlist.txt -t 6
```

## Annotation (optional)

```bash
DEBKS anno -c merge_pos.txt  -m hg19.fa  -g gencode.v19.annotation.sort.gtf.gz
```

## DEC

```bash
DEBKS dec -c merge_circ.txt -l merge_linear.txt -t 20 -e circ_len.txt
DEBKS dec -c merge_circ.txt -l linear_jun.txt -t 20 -e circ_len.txt #optional
```