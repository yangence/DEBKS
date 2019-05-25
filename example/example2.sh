#!/usr/bin/sh
## Chimeric and spliced junctions files provided
### Space required: 4G
echo -e 'Example 1: Chimeric and spliced junction counts provided\nAttention: 4G storage required'
tar -xzf STAR_Results.tar.gz
mkdir -p ./DEBKS_example_2
cd ./DEBKS_example_2
mkdir -p ./genome ./annotation
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

cd ..
cat <<EOF >sample_1_CJ.txt
../STAR_Results/sample_1.Rep1.Chimeric.out.junction
../STAR_Results/sample_1.Rep2.Chimeric.out.junction
EOF

cat <<EOF >sample_2_CJ.txt
../STAR_Results/sample_2.Rep1.Chimeric.out.junction
../STAR_Results/sample_2.Rep2.Chimeric.out.junction
EOF

cat <<EOF >sample_1_SJ.txt
../STAR_Results/sample_1.Rep1.SJ.out.tab
../STAR_Results/sample_1.Rep2.SJ.out.tab
EOF

cat <<EOF >sample_2_SJ.txt
../STAR_Results/sample_2.Rep1.SJ.out.tab
../STAR_Results/sample_2.Rep2.SJ.out.tab
EOF
echo 'DEBKS begin'
DEBKS -g ./genome/hg19.fa \
	-s1CJ sample_1_CJ.txt -s2CJ sample_2_CJ.txt -s1SJ sample_1_SJ.txt -s2SJ sample_2_SJ.txt \
	-gtf ./annotation/gencode.v19.annotation.gtf \
	-o ./DEBKS_Results/ -t 10 -read single -len 100 -c 0.01 -a 10
echo 'Down'
