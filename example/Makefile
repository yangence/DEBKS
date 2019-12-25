.PHONY: clean getGene getCount test1 test2 simu

SIGMA=0.01
GTF=/media/data4/lzl/annotation/GENCODE/humanV19/gencode.v19.annotation.gtf
GENE=./simu/gene.gtf
REF_DIR=/media/data4/lzl/genome/hg19_chr/
REF=/media/data4/lzl/genome/hg19/hg19.fa
STARINDEX=/media/data4/lzl/genome/hg19_STAR_149/
getGene:
	mkdir -p ./simu
	if [ -e $(GENE) ];then echo 'Pass extract gene tag';else perl -alne 'if($$F[2] eq "gene"){print $$_}' $(GTF)>$(GENE);fi

getCount: getGene
	if [ -e ./simu/coverage_$(SIGMA)_g1_rep1.txt ];then echo 'Pass simulation junction';else python sim_circ.py $(SIGMA) ./simu $(GENE);fi
simu: getCount
	perl revised_simulator.pl  ./simu/coverage_*_g1_rep1.txt  ./simu/g1_rep1 $(REF_DIR) $(GTF) 
	perl revised_simulator.pl  ./simu/coverage_*_g1_rep2.txt  ./simu/g1_rep2 $(REF_DIR) $(GTF) 
	perl revised_simulator.pl  ./simu/coverage_*_g1_rep3.txt  ./simu/g1_rep3 $(REF_DIR) $(GTF) 
	perl revised_simulator.pl  ./simu/coverage_*_g2_rep1.txt  ./simu/g2_rep1 $(REF_DIR) $(GTF) 
	perl revised_simulator.pl  ./simu/coverage_*_g2_rep2.txt  ./simu/g2_rep2 $(REF_DIR) $(GTF) 
	perl revised_simulator.pl  ./simu/coverage_*_g2_rep3.txt  ./simu/g2_rep3 $(REF_DIR) $(GTF) 


test1:
	DEBKS -g $(REF) \
	-STARindex $(STARINDEX) \
	-s1 g1.txt -s2 g2.txt \
	-gtf $(GTF) \
	-o DEBKS_test1/ -t 10 -read pair -len 150

mapping:
	STAR --genomeLoad LoadAndKeep --chimSegmentMin 6 --chimOutType WithinBAM --runThreadN 10  --outSAMtype BAM Unsorted --alignSJDBoverhangMin 6 --alignSJoverhangMin 6  --genomeDir $(STARINDEX) --chimJunctionOverhangMin 6 --outSJfilterOverhangMin -1 6 -1 -1 --outFileNamePrefix ./simu/g1_rep1. --readFilesIn ./simu/g1_rep1_1.fq ./simu/g1_rep1_2.fq
	STAR --genomeLoad LoadAndKeep --chimSegmentMin 6 --chimOutType WithinBAM --runThreadN 10  --outSAMtype BAM Unsorted --alignSJDBoverhangMin 6 --alignSJoverhangMin 6  --genomeDir $(STARINDEX) --chimJunctionOverhangMin 6 --outSJfilterOverhangMin -1 6 -1 -1 --outFileNamePrefix ./simu/g1_rep2. --readFilesIn ./simu/g1_rep2_1.fq ./simu/g1_rep2_2.fq
	STAR --genomeLoad LoadAndKeep --chimSegmentMin 6 --chimOutType WithinBAM --runThreadN 10  --outSAMtype BAM Unsorted --alignSJDBoverhangMin 6 --alignSJoverhangMin 6  --genomeDir $(STARINDEX) --chimJunctionOverhangMin 6 --outSJfilterOverhangMin -1 6 -1 -1 --outFileNamePrefix ./simu/g1_rep3. --readFilesIn ./simu/g1_rep3_1.fq ./simu/g1_rep3_2.fq
	STAR --genomeLoad LoadAndKeep --chimSegmentMin 6 --chimOutType WithinBAM --runThreadN 10  --outSAMtype BAM Unsorted --alignSJDBoverhangMin 6 --alignSJoverhangMin 6  --genomeDir $(STARINDEX) --chimJunctionOverhangMin 6 --outSJfilterOverhangMin -1 6 -1 -1 --outFileNamePrefix ./simu/g2_rep1. --readFilesIn ./simu/g2_rep1_1.fq ./simu/g2_rep1_2.fq
	STAR --genomeLoad LoadAndKeep --chimSegmentMin 6 --chimOutType WithinBAM --runThreadN 10  --outSAMtype BAM Unsorted --alignSJDBoverhangMin 6 --alignSJoverhangMin 6  --genomeDir $(STARINDEX) --chimJunctionOverhangMin 6 --outSJfilterOverhangMin -1 6 -1 -1 --outFileNamePrefix ./simu/g2_rep2. --readFilesIn ./simu/g2_rep2_1.fq ./simu/g2_rep2_2.fq
	STAR --genomeLoad LoadAndKeep --chimSegmentMin 6 --chimOutType WithinBAM --runThreadN 10  --outSAMtype BAM Unsorted --alignSJDBoverhangMin 6 --alignSJoverhangMin 6  --genomeDir $(STARINDEX) --chimJunctionOverhangMin 6 --outSJfilterOverhangMin -1 6 -1 -1 --outFileNamePrefix ./simu/g2_rep3. --readFilesIn ./simu/g2_rep3_1.fq ./simu/g2_rep3_2.fq
	STAR --genomeLoad Remove --outFileNamePrefix ./simu/temp. --genomeDir $(STARINDEX)

test2:
	DEBKS -g $(REF) \
	-s1CJ g1_CJ.txt -s2CJ g2_CJ.txt -s1SJ g1_SJ.txt -s2SJ g2_SJ.txt \
	-gtf $(GTF) \
	-o DEBKS_test2/ -t 10 -read pair -len 150

clean:
	rm -rf simu DEBKS_test1 DEBKS_test2
