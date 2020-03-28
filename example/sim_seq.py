#!/usr/bin/env python
import pandas as pd,numpy as np,pysam,os,sys,gzip

matFile=sys.argv[1]
outPrefix=sys.argv[2]
#isAS=1
#if len(sys.argv)>3:
#    isAS=0
geneMat=pd.read_csv(matFile,header=None,sep='\t')

def readGTFfile(fileName):
    if not os.path.exists(fileName):
        sys.exit('ERROR: %s is not exist!!!' % fileName)
    try:
        tabixfile = pysam.TabixFile(fileName)
        return(tabixfile)
    except:
        sys.exit('ERROR: make sure %s is sorted and tabix indexed!!!' % gtfFile)

gtfFile='/media/data4/lzl/annotation/GENCODE/humanV19/gencode.v19.annotation.sort.gtf.gz'
faFile='/media/data4/lzl/genome/hg19/hg19.fa'
readLen=150
errRate=0.0001
fragmentLen=350
read1_file=outPrefix+'_1.fq.gz'
read2_file=outPrefix+'_2.fq.gz'
info_file=outPrefix+'_info.txt'
detail_file=outPrefix+'_detail.txt'
gtfFile=readGTFfile(gtfFile)

faFile=pysam.FastaFile(faFile)

quality=''.join(['F' for i in range(readLen)])

read1_out=gzip.open(read1_file,'wb')
read2_out=gzip.open(read2_file,'wb')
info_out=open(info_file,'w')
detail_out=open(detail_file,'w')
gene_dict={}
gene_chr={}
for gtf in gtfFile.fetch(parser=pysam.asGTF()):
    if gtf.feature=='gene':
        gene_chr[gtf.gene_id]=gtf.contig
    if gtf.feature=='exon':
        if gene_dict.__contains__(gtf.gene_id):
            if gene_dict[gtf.gene_id].__contains__(gtf.transcript_id):
                gene_dict[gtf.gene_id][gtf.transcript_id].append([gtf.start,gtf.end])
            else:
                gene_dict[gtf.gene_id][gtf.transcript_id]=[[gtf.start,gtf.end]]
        else:
            gene_dict[gtf.gene_id]={}
            gene_dict[gtf.gene_id][gtf.transcript_id]=[[gtf.start,gtf.end]]
            

exon_maxNum={}
gene_can3=[]
gene_can5=[]
for i in gene_dict.keys():
    exon_maxNum[i]=0
    for j in gene_dict[i].keys():
        exon_maxNum[i]=max(exon_maxNum[i],len(gene_dict[i][j]))
    if exon_maxNum[i]>2:
        gene_can3.append(i)
    if exon_maxNum[i]>4:
        gene_can5.append(i)

        
def getSeq(arr,contig,genome):
    seq=''
    for i in arr:
        seq+=genome.fetch(contig,i[0],i[1])
    return(seq)

def revCom(seq):
    return(seq[::-1].upper().replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper())

def getErr(read1,read2,errNum):
    base2base={'A':['T','G','C'],'T':['A','G','C'],'G':['T','A','C'],'C':['T','G','A'],'N':['A','G','T']}
    readLen=len(read1)
    pos=np.random.choice(2*readLen,errNum,replace=False).tolist()
    read1=list(read1)
    read2=list(read2)
    for i in pos:
        if i <readLen:
            read1[i]=base2base[read1[i]][np.random.randint(3)]
        else:
            index=i-readLen
            read2[index]=base2base[read2[index]][np.random.randint(3)]
    return(''.join(read1),''.join(read2))

def simulate_linear(cov,readLen,seq,errRate,fragmentLen_o,out1,out2,id):
    seqLen=len(seq)
    numRead=int(cov*seqLen/readLen/2)
    errBaseNum=int(errRate*numRead*readLen*2)
    eachReaderr_dict={}
    for i in range(errBaseNum):
        tmpErr=np.random.randint(numRead)
        if eachReaderr_dict.__contains__(tmpErr):
            eachReaderr_dict[tmpErr]+=1
        else:
            eachReaderr_dict[tmpErr]=1
    for i in range(numRead):
        fragmentLen=int(np.random.normal(fragmentLen_o,30))
        if fragmentLen<=readLen or fragmentLen>=seqLen:
            continue
        startLoci=np.random.randint(seqLen-fragmentLen+1)
        start2Loci=startLoci+fragmentLen-readLen
        if np.random.randint(2):
            read1=seq[startLoci:(startLoci+readLen)].upper()
            read2=revCom(seq[start2Loci:(start2Loci+readLen)])
        else:
            read1=revCom(seq[start2Loci:(start2Loci+readLen)])
            read2=seq[startLoci:(startLoci+readLen)].upper()
        if eachReaderr_dict.__contains__(i):
            read1,read2=getErr(read1,read2,eachReaderr_dict[i])
        readName='@'+id+'000'+str(i)
        out1.write(bytes("%s\n%s\n+\n%s\n" %(readName+'/1'+' length='+str(readLen), read1,quality),encoding='ASCII'))
        out2.write(bytes("%s\n%s\n+\n%s\n" %(readName+'/2'+' length='+str(readLen), read2,quality),encoding='ASCII'))

def simulate_circ(cov,readLen,seq,errRate,fragmentLen_o,out1,out2,id):
    seqLen=len(seq)
    numRead=int(cov*seqLen/readLen/2)
    seq=''.join([seq for i in range(10)])
    seqLen=10*seqLen
    errBaseNum=int(errRate*numRead*readLen*2)
    eachReaderr_dict={}
    for i in range(errBaseNum):
        tmpErr=np.random.randint(numRead)
        if eachReaderr_dict.__contains__(tmpErr):
            eachReaderr_dict[tmpErr]+=1
        else:
            eachReaderr_dict[tmpErr]=1
    for i in range(numRead):
        fragmentLen=int(np.random.normal(fragmentLen_o,30))
        if fragmentLen<=readLen or fragmentLen>=seqLen:
            continue
        startLoci=np.random.randint(seqLen-fragmentLen+1)
        start2Loci=startLoci+fragmentLen-readLen
        if np.random.randint(2):
            read1=seq[startLoci:(startLoci+readLen)].upper()
            read2=revCom(seq[start2Loci:(start2Loci+readLen)])
        else:
            read1=revCom(seq[start2Loci:(start2Loci+readLen)])
            read2=seq[startLoci:(startLoci+readLen)].upper()
        if eachReaderr_dict.__contains__(i):
            read1,read2=getErr(read1,read2,eachReaderr_dict[i])
        readName='@'+id+':'+str(i)
        out1.write(bytes("%s\n%s\n+\n%s\n" %(readName+'/1'+' length='+str(readLen), read1,quality),encoding='ASCII'))
        out2.write(bytes("%s\n%s\n+\n%s\n" %(readName+'/2'+' length='+str(readLen), read2,quality),encoding='ASCII'))
        
def getCirc2(exon,count):
    exonNum=len(exon)
    skipPos=count % (exonNum-2) +1
    return(exon[0:skipPos]+exon[(skipPos+1):exonNum])
def getCirc1(exon,count):
    exonNum=len(exon)
    choicePos=count % (exonNum-2) +1
    if choicePos % 2:
        return(exon[0:choicePos+1])
    else:
        return(exon[(choicePos+1):exonNum])
def exonTrans(exon):
    newExon=','.join([str(i[0]+1)+'-'+str(i[1]) for i in exon])
    return(newExon)

count=0
circ_list=[]
for i in range(geneMat.shape[0]):
    eachMat=geneMat.iloc[i,:].tolist()
    if count % 2:
        isAS=1
    else:
        isAS=0
    if isAS:
        gene_can=gene_can5
    else:
        gene_can=gene_can3
    if eachMat[0] not in gene_can:
        continue
    linearCov=eachMat[2]
    circCov=eachMat[1]
    circ1Cov=eachMat[3]
    circ2Cov=circCov-circ1Cov
    # choice transcript
    tran_dict=gene_dict[eachMat[0]]
    tranID=''
    for j in tran_dict.keys():
        if isAS:
            if len(tran_dict[j])>=5:
                tranID=j
                break
        else:
            if len(tran_dict[j])>=3:
                tranID=j
                break
    exon_arr=tran_dict[tranID]
    contig=gene_chr[eachMat[0]]
    linearSeq=getSeq(exon_arr,contig,faFile)
    if len(linearSeq)<fragmentLen+100:
        continue
    circ1_exon=[exon_arr[x] for x in range(1,len(exon_arr)-1)] # full exon
    if not isAS:
        if len(circ1_exon)>2:
            circ1_exon=getCirc1(circ1_exon,count)
    circ1Seq=getSeq(circ1_exon,contig,faFile)
    #if len(circ1Seq)<readLen:
    #    continue
    # only skip one exon for circ2
    if isAS:
        circ2_exon=getCirc2(circ1_exon,count) #skip internal exon
        circ2Seq=getSeq(circ2_exon,contig,faFile)
        #if len(circ2Seq)<readLen:
        #    continue
    circID=contig+':'+str(circ1_exon[0][0]+1)+'|'+str(circ1_exon[-1][1])
    if circID in circ_list:
        continue
    else:
        circ_list.append(circID)
    circStart=','.join([str(i[0]+1) for i in circ1_exon])
    circEnd=','.join([str(i[1]) for i in circ1_exon])
    info_out.write("%s\t%s\t%s\t%s\t%s\t%d\n" % (eachMat[0],tranID,circID,circStart,circEnd,isAS))
    detail_out.write("%s\t%s\t%s\t%d\n" % (eachMat[0],circID,exonTrans(circ1_exon),isAS))
    simulate_linear(linearCov,readLen,linearSeq,errRate,fragmentLen,read1_out,read2_out,'simulate:00:'+str(count))
    if isAS:
        simulate_circ(circ1Cov,readLen,circ1Seq,errRate,fragmentLen,read1_out,read2_out,'simulate:01:'+str(count))
        simulate_circ(circ2Cov,readLen,circ2Seq,errRate,fragmentLen,read1_out,read2_out,'simulate:02:'+str(count))
        detail_out.write("%s\t%s\t%s\t%d\n" % (eachMat[0],circID,exonTrans(circ2_exon),isAS))
    else:
        simulate_circ(circCov,readLen,circ1Seq,errRate,fragmentLen,read1_out,read2_out,'simulate:03:'+str(count))
    count+=1
    if count==10000:
        break
read1_out.close()
read2_out.close()
