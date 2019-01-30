#!/usr/bin/sh
## FASTQ files provided
### Space required: 95G
echo -e 'Example 1: Raw RNA-seq reads provided\nAttention: 95G storage required'
mkdir -p ./DEBKS_example_1
cd ./DEBKS_example_1
mkdir -p ./rawData ./genome ./annotation ./hg19_STAR_99
cd rawData
echo 'Download RNA-seq data from SRA'
wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660797/SRR1660797.sra
wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660798/SRR1660798.sra
wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660799/SRR1660799.sra
wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660800/SRR1660800.sra

wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660790/SRR1660790.sra
wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660794/SRR1660794.sra
wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660801/SRR1660801.sra
wget -c ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR166/SRR1660802/SRR1660802.sra
echo 'Down'
echo 'Transfrom sra format to FASTQ'
fastq-dump --gzip SRR1660797.sra &
fastq-dump --gzip SRR1660798.sra &
fastq-dump --gzip SRR1660799.sra &
fastq-dump --gzip SRR1660800.sra &

fastq-dump --gzip SRR1660790.sra &
fastq-dump --gzip SRR1660794.sra &
fastq-dump --gzip SRR1660801.sra &
fastq-dump --gzip SRR1660802.sra &

wait
echo 'Down'
rm *.sra

cd ..
cat <<EOF >sample_1.txt
./rawData/SRR1660797.fastq.gz
./rawData/SRR1660798.fastq.gz
./rawData/SRR1660799.fastq.gz
./rawData/SRR1660800.fastq.gz
EOF

cat <<EOF >sample_2.txt
./rawData/SRR1660790.fastq.gz
./rawData/SRR1660794.fastq.gz
./rawData/SRR1660801.fastq.gz
./rawData/SRR1660802.fastq.gz
EOF

cd ./genome/
echo 'Download hg19 reference genome'
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -zxf chromFa.tar.gz
cat *.fa > hg19.fa
echo 'Down'
cd ../annotation/
echo 'Download GENCODE annotation'
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gzip -d gencode.v19.annotation.gtf.gz
echo 'down'
cd ../
echo 'Build STAR alignment index'
STAR --runMode genomeGenerate --runThreadN 20 \
	--genomeFastaFiles ./genome/hg19.fa \
	--sjdbGTFfile ./annotation/gencode.v19.annotation.gtf \
	--sjdbOverhang 99 \
	--genomeDir ./hg19_STAR_99
echo 'Down'
echo 'DEBKS begin'
python ../../DEBKS_main.py -g ./genome/hg19.fa \
	-s1 sample_1.txt -s2 sample_2.txt -STARindex ./hg19_STAR_99/ \
	-gtf ./annotation/gencode.v19.annotation.gtf \
	-o ./DEBKS_Results/ -t 30 -read single -len 100 -c 0.01 -a 10
echo 'Down'
#Optional
cat <<EOF >sample_1_bam.txt
./DEBKS_Results/SAMPLE_1/REP_1/Aligned.out.bam
./DEBKS_Results/SAMPLE_1/REP_2/Aligned.out.bam
./DEBKS_Results/SAMPLE_1/REP_3/Aligned.out.bam
./DEBKS_Results/SAMPLE_1/REP_4/Aligned.out.bam
EOF

cat <<EOF >sample_2_bam.txt
./DEBKS_Results/SAMPLE_2/REP_1/Aligned.out.bam
./DEBKS_Results/SAMPLE_2/REP_2/Aligned.out.bam
./DEBKS_Results/SAMPLE_2/REP_3/Aligned.out.bam
./DEBKS_Results/SAMPLE_2/REP_4/Aligned.out.bam
EOF
echo 'Plot of differential PBSI'
python ../../src/DEBKS_plot.py DEBKS_plot.py -s1 sample_1_bam.txt -s2 sample_2_bam.txt -o ./DEBKS_Results/ \
	-BS ./DEBKS_Results/DEBKS_output/DEBKS_results.txt
echo 'Down'
