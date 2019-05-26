.PHONY: clean sample test1 test2

sample:
        tar -xf STAR_Results.tar.gz
        echo -e "./STAR_Results/sample_1.Rep1.Chimeric.out.junction\n./STAR_Results/sample_1.Rep2.Chimeric.out.junction">sample_1_CJ.txt
        echo -e "./STAR_Results/sample_2.Rep1.Chimeric.out.junction\n./STAR_Results/sample_2.Rep2.Chimeric.out.junction">sample_2_CJ.txt
        echo -e "./STAR_Results/sample_1.Rep1.SJ.out.tab\n./STAR_Results/sample_1.Rep2.SJ.out.tab">sample_1_SJ.txt
        echo -e "./STAR_Results/sample_2.Rep1.SJ.out.tab\n./STAR_Results/sample_2.Rep2.SJ.out.tab">sample_2_SJ.txt

test1:
         DEBKS -g ${fasta} \
         -s1CJ sample_1_CJ.txt -s2CJ sample_2_CJ.txt -s1SJ sample_1_SJ.txt -s2SJ sample_2_SJ.txt \
         -gtf ${gtf} \
         -o DEBKS_test1/ -t 1 -read single -len 100 -c 0.01 -a 10

test2:
         DEBKS -g ${fasta} \
         -s1CJ sample_1_CJ.txt -s2CJ sample_2_CJ.txt -s1SJ sample_1_SJ.txt -s2SJ sample_2_SJ.txt \
         -gtf ${gtf} \
         -o DEBKS_test2/ -t 10 -read single -len 100 -c 0.01 -a 10

clean:
        rm sample*.txt
        rm -fr STAR_Results
        rm -fr DEBKS_test*
