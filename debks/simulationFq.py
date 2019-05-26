"""
author: Zelin Liu
email: zlliu@bjmu.edu.cn
license: GPL3
detail: Simulate ribo-zero RNA-seq data base on annotation file and candidate circRNA
"""
import sys,os,pysam,random
import pandas as pd
import numpy as np
if sys.version_info[0] < 3:
    print ("Python Version error: must use python 3 (not supporting python 2)");
    sys.exit(-1);
import subprocess;
# system command
def commandRun(cmd):
    p= subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE);
    status=p.wait()
    status=p.returncode
    outputALL=p.communicate()
    output=''
    for i in outputALL:
        output+=i.decode()
    return(status,output)
#parameters
psiDelt=0.01
readType=''
genomeFasta=''
circDB=''
CIRCanno=''
outDir=''
insertsize=200
fragmentdev=50
readLength=150
strandSpecific=True
rollCircLen=1000
Replicate_Num=3
validArgList=['-g','-anno','-bed','-read','-o','-delt','-rep']

for argIndex in range(1,len(sys.argv)): ## going through the all parameters 
    if(sys.argv[argIndex][0]=='-' and sys.argv[argIndex] not in validArgList): ## incorrect argument
        print('Not valid argument: %s' % sys.argv[argIndex]);
        print('Please provide valid arguments.');
        sys.exit(-1);
        
for paramIndex in range(1,len(sys.argv)): ## going through the all parameters
    if(sys.argv[paramIndex] == '-g'):  
        paramIndex += 1;  ## increase index
        genomeFasta = sys.argv[paramIndex];
    elif(sys.argv[paramIndex] == '-anno'): 
        paramIndex += 1;  ## increase index
        CIRCanno = sys.argv[paramIndex];
    elif(sys.argv[paramIndex] == '-bed'):  
        paramIndex += 1;  ## increase index
        circDB = sys.argv[paramIndex]; # start(0 based) end(1 based) bed file from circBase
    elif(sys.argv[paramIndex] == '-o'):  
        paramIndex += 1;  ## increase index
        outDir = sys.argv[paramIndex];
    elif (sys.argv[paramIndex] == '-read'):  ## readtype, single or paired
        paramIndex += 1;  ## increase index
        readType = sys.argv[paramIndex];
        if readType not in ['single','paired']:
            print('-read only accept "single" or "paired" value')
            sys.exit(1)
    elif(sys.argv[paramIndex] == '-delt'):  
        paramIndex += 1;  ## increase index
        psiDelt = float(sys.argv[paramIndex]);
    elif(sys.argv[paramIndex] == '-rep'):  
        paramIndex += 1;  ## increase index
        Replicate_Num = int(sys.argv[paramIndex]);

if (genomeFasta == '' or  CIRCanno ==''  or circDB=='' or readType=='' or outDir==''):
    print('Not enough arguments!!');
    print('Usage (with fastq files):\n\tpython simulationFq.py -g genomeFasta -o outDir -read readType -anno CIRCanno');
    sys.exit(-1);

SEPE = 'SE'; ## single-end or paired
if readType=="paired": ## single-end
    SEPE='PE';
    fragTargetSize=2*readLength + insertsize
else:
    fragTargetSize=readLength + insertsize
outPath = os.path.abspath(outDir)
os.system('mkdir -p '+ outPath)
tempPath=outPath + '/temp'
os.system('mkdir -p '+ tempPath)
circfile=pd.read_table(circDB,header=None,comment='#')
bsJ=pd.DataFrame('',index=range(circfile.iloc[:,0].size),columns=range(6))
bsJ.iloc[:,0:3]=circfile.iloc[:,0:3]
bsJ.iloc[:,3]=['FUSIONJUNC_'+str(i)+'/100' for i in range(circfile.iloc[:,0].size)]
bsJ.iloc[:,4]=0
bsJ.iloc[:,5]='+'
bsJ.to_csv(tempPath+'/back_spliced_junction.bed',index=False,sep='\t',header=False)

if not os.path.exists(tempPath+'/circularRNA_known.txt'):
    cmd='CIRCexplorer2 annotate -r '+CIRCanno+ ' -g '+ genomeFasta + ' -b '+tempPath+'/back_spliced_junction.bed -o '+tempPath+'/circularRNA_known.txt'
    status, output=commandRun(cmd)
    print('####  down CIRCexplorer2 annotate. status: %s ####' % status)
    if (int(status)!=0): ## it did not go well
        print("error in CIRCexplorer2 annotate: %s" % status);
        print("error detail: %s" % output);
        sys.exit(-1)

circRNA=pd.read_table(tempPath+'/circularRNA_known.txt',header=None)
annoTrans=pd.read_table(CIRCanno,header=None)
annoCirc=circRNA.merge(annoTrans,left_on=15,right_on=1)
noDupIndex=[not i for i in  annoCirc.iloc[:,14].duplicated().tolist()] # remove circRNA from same host gene
annoCirc=annoCirc[noDupIndex]
annoCirc.to_csv(tempPath+'/annoCirc.txt',index=False,sep='\t',header=False)

# read genome Fasta
genome=pysam.FastaFile(genomeFasta)

# Define function
def revC(seq):
    seq=seq[::-1].upper().replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper()
    return(seq)


def getIsoform(x):
    chr=x[0]
    circStart=x[1]
    circEnd=x[2]
    strand=x[5]
    circSize=[int(i) for i in x[10].split(',')]
    circOffset=[int(i) for i in x[11].split(',')]
    gene=x[14]
    isoform=x[15]
    transStart=[int(i) for i in x[27][0:(len(x[27])-1)].split(',')]
    transEnd=[int(i) for i in x[28][0:(len(x[28])-1)].split(',')]
    circFasta=''
    transFasta=''
    for i in range(len(circSize)):
        tmpStart=circStart+circOffset[i]
        tmpEnd=tmpStart+circSize[i]
        circFasta+=genome.fetch(reference=chr,start=tmpStart,end=tmpEnd)
    for i in range(len(transStart)):
        transFasta+=genome.fetch(reference=chr,start=transStart[i],end=transEnd[i])
    if strand == '-':
        circFasta=revC(circFasta)
        transFasta=revC(transFasta)
    return(circFasta,transFasta)
    

def getFragment(seq,targetSize):
    newSeq=''
    if len(seq) < targetSize-fragmentdev:
        return [seq]
    size=len(seq)
    fgLength=targetSize+random.randint(-fragmentdev,fragmentdev)
    if random.randint(0,1)==1:
        end=min(fgLength,size)
        reSeq=[seq[0:end]]
        if end < size:
            newSeq=seq[end:]
    else:
        start=max(size-fgLength,0)
        reSeq=[seq[start:]]
        if start > 0:
            newSeq=seq[0:start]
    return(reSeq+getFragment(newSeq,targetSize))
    
def getReadPE(seq):
    return(revC(seq[(len(seq)-readLength):len(seq)]),seq[0:readLength])

def getReadSE(seq):
    return(revC(seq[(len(seq)-readLength):len(seq)]))
    
def simuLinear(seq,n):
    seq=seq+'AAAAAAAAAAAA'
    if SEPE == 'SE':
        if len(seq) < fragTargetSize-fragmentdev:
            return([])
        else:
            read=[]
            for i in range(n):
                fragment=[i for i in getFragment(seq,fragTargetSize) if len(i) >= fragTargetSize-fragmentdev ]
                for i in fragment:
                    read.append(getReadSE(i))
            return(read)
    else:
        if len(seq) < fragTargetSize-fragmentdev:
            return([],[])
        else:
            read1=[]
            read2=[]
            for i in range(n):
                fragment=[i for i in getFragment(seq,fragTargetSize) if len(i) >= fragTargetSize-fragmentdev ]
                for i in fragment:
                    tmp1,tmp2=getReadPE(i)
                    read1.append(tmp1)
                    read2.append(tmp2)
            return(read1,read2)

def getCircRollseq(seq):
    oLen=len(seq)
    if oLen < fragTargetSize-fragmentdev:
        newSeq=''
        lStart=random.randint(0,oLen-1)
        rEnd=random.randint(0,oLen-1)
        newrollCircLen=rollCircLen-(oLen-lStart)-rEnd
        n=newrollCircLen // oLen
        for i in range(n):
            newSeq+=seq
        newSeq=seq[lStart:]+newSeq+seq[0:rEnd]
    else:
        breakpoint=random.randint(0,oLen-1)
        newSeq=seq[breakpoint:]+seq[0:breakpoint]
    return(newSeq)
    

def simuCirc(seq,n):
    if SEPE =='SE':
        read=[]
        for i in range(n):
            seq=getCircRollseq(seq)
            fragment=[i for i in getFragment(seq,fragTargetSize) if len(i) >= fragTargetSize-fragmentdev ]
            for i in fragment:
                read.append(getReadSE(i))
        return(read)
    else:
        read1=[]
        read2=[]
        for i in range(n):
            seq=getCircRollseq(seq)
            fragment=[i for i in getFragment(seq,fragTargetSize) if len(i) >= fragTargetSize-fragmentdev ]
            for i in fragment:
                tmp1,tmp2=getReadPE(i)
                read1.append(tmp1)
                read2.append(tmp2)
        return(read1,read2)
            

def getBSnum(LSnum,psi):
    if psi < 0.005:
        psi = 0.005
    BSnum=LSnum*(1-psi)/psi
    return(int(BSnum))

def getIndex(psi):
    index=[]
    for i in range(len(psi)):
        if psi[i]>1 or psi[i]<0:
            index.append(i)
    return(index)
# run

percentDE=0.05


percentDEnum=int(annoCirc.iloc[:,0].size*percentDE)
DEsample=np.array(random.sample(range(annoCirc.iloc[:,0].size),percentDEnum))
FinalDEsample=DEsample
upDE=int(len(DEsample)/2)
DEsampleUp=DEsample[random.sample(range(len(DEsample)),upDE)]

psi1= np.random.random(annoCirc.iloc[:,0].size)*0.95+0.05
# not DE |psi1-psi2|<percentDE
psi2=psi1+np.random.random(len(psi1))/10 -0.05
psiIndex=getIndex(psi2)
while len(psiIndex)>0:
    psi2[psiIndex]=psi1[psiIndex]+np.random.random(len(psiIndex))/10 -0.05
    psiIndex=getIndex(psi2)
# DE |psi1-psi2|>=0.05

psi2[DEsample]=psi1[DEsample]-np.random.random(len(DEsample))*0.95-0.05
psiIndex=getIndex(psi2[DEsample])
threshold=0.95
while len(psiIndex) > 0:
    psi2[DEsample[psiIndex]]=psi1[DEsample[psiIndex]]-np.random.random(len(psiIndex))*threshold-0.05
    threshold=threshold*0.9
    psiIndex=getIndex(psi2[DEsample])
    

psi2[DEsampleUp]=psi1[DEsampleUp]+np.random.random(len(DEsampleUp))*0.95+0.05
psiIndex=getIndex(psi2[DEsampleUp])
threshold=0
while len(psiIndex) > 0:
    psi2[DEsampleUp[psiIndex]]=psi1[DEsampleUp[psiIndex]]+np.random.random(len(psiIndex))*0.95+0.05+threshold
    psiIndex=getIndex(psi2[DEsampleUp])
    threshold=threshold-0.01
    

#psi matrix
psi1mat=pd.DataFrame(0,index=range(annoCirc.iloc[:,15].size),columns=range(Replicate_Num))
for i in range(Replicate_Num):
    psi1mat.iloc[:,i]=np.random.normal(psi1, psiDelt).tolist()
    psiIndex=getIndex(psi1mat.iloc[:,i].tolist())
    while len(psiIndex) > 0:
        psi1mat.iloc[psiIndex,i]=np.random.normal(psi1[psiIndex], psiDelt).tolist()
        psiIndex=getIndex(psi1mat.iloc[:,i].tolist())
        
psi1mat['psio']=psi1
psi1mat.to_csv(outPath+'/psi1mat.txt',index=False,sep='\t',header=False)
psi2mat=pd.DataFrame(0,index=range(annoCirc.iloc[:,15].size),columns=range(Replicate_Num))
for i in range(Replicate_Num):
    psi2mat.iloc[:,i]=np.random.normal(psi2, psiDelt).tolist()
    psiIndex=getIndex(psi2mat.iloc[:,i].tolist())
    while len(psiIndex) > 0:
        psi2mat.iloc[psiIndex,i]=np.random.normal(psi2[psiIndex], psiDelt).tolist()
        psiIndex=getIndex(psi2mat.iloc[:,i].tolist())
        
psi2mat['psio']=psi2
psi2mat.to_csv(outPath+'/psi2mat.txt',index=False,sep='\t',header=False)
#sys.exit(0)
#transcipt number matrix
transcriptNum=pd.DataFrame(0,index=range(len(annoCirc.iloc[:,15].unique())),columns=range(2*Replicate_Num))
transcriptNum.rename(index=dict(zip(transcriptNum.index,annoCirc.iloc[:,15].unique())),inplace=True)
N=transcriptNum.iloc[:,0].size
rawIsoNum=np.random.randint(50,500,N).tolist()
for i in range(2*Replicate_Num):
    transcriptNum.iloc[:,i]= [i for i in map(int,np.random.normal(rawIsoNum, 1).tolist())]
transcriptNum.to_csv(outPath+'/transcriptNum.txt',index=False,sep='\t',header=False)

####simulate DE gene
randTranscript=np.array(random.sample(range(N),int(0.05*N)))
#DE 1.5< fold change< 5
randFC=(np.random.random(len(randTranscript))*3.5+1.5).tolist()
#up in sample 1
randUp=randTranscript[random.sample(range(len(randTranscript)),int(len(randTranscript)/2))]
randDown=[i for i in randTranscript if i not in randUp]
print(len(randDown))
for i in range(len(randUp)):
    for j in range(Replicate_Num):
        transcriptNum.iat[randUp[i],j]=transcriptNum.iat[randUp[i],j]*randFC[i]
        
for i in range(len(randDown)):
    for j in range(Replicate_Num):
        transcriptNum.iat[randDown[i],Replicate_Num+j]=transcriptNum.iat[randDown[i],Replicate_Num+j]*randFC[i+len(randUp)]

        

#simulation begin
readName='@A00001:001:'
qualityLine=''.join(['J' for i in range(readLength)])
for i in range(Replicate_Num):      #sample_1
    if SEPE == 'SE':
        fastqFile='sample_1_rep'+str(i+1)
        fileName=outPath+'/sample_1_rep'+str(i+1)+'.fastq'
        exec(fastqFile+'=open(\''+fileName+'\',\'w\')')
        fastqFile='sample_2_rep'+str(i+1)
        fileName=outPath+'/sample_2_rep'+str(i+1)+'.fastq'
        exec(fastqFile+'=open(\''+fileName+'\',\'w\')')
    else:
        fastq1File='sample_1_rep'+str(i+1)+'_R1'
        fastq2File='sample_1_rep'+str(i+1)+'_R2'
        file1Name=outPath+'/sample_1_rep'+str(i+1)+'.R1.fastq'
        file2Name=outPath+'/sample_1_rep'+str(i+1)+'.R2.fastq'
        exec(fastq1File+'=open(\''+file1Name+'\',\'w\')')
        exec(fastq2File+'=open(\''+file2Name+'\',\'w\')')
        fastq1File='sample_2_rep'+str(i+1)+'_R1'
        fastq2File='sample_2_rep'+str(i+1)+'_R2'
        file1Name=outPath+'/sample_2_rep'+str(i+1)+'.R1.fastq'
        file2Name=outPath+'/sample_2_rep'+str(i+1)+'.R2.fastq'
        exec(fastq1File+'=open(\''+file1Name+'\',\'w\')')
        exec(fastq2File+'=open(\''+file2Name+'\',\'w\')')
        
readName=('@'+annoCirc.iloc[:,0]+':'+annoCirc.iloc[:,1].map(str)+'|'+annoCirc.iloc[:,2].map(str)+':').tolist()
for i in range(10000):
    print(i)
    circSeq,linearSeq=getIsoform(annoCirc.iloc[i,:].tolist())
    transcriptName=annoCirc.iat[i,15]
    for j in range(Replicate_Num):      #sample_1
        tranN=transcriptNum.ix[transcriptName,j]
        circN=getBSnum(tranN,psi1mat.iloc[i,j])
        if SEPE=='SE':
            fastqFile='sample_1_rep'+str(i+1)
            circRead=simuCirc(circSeq,circN)
            linearRead=simuLinear(linearSeq,tranN)
            
            for p in range(len(circRead)):
                eachreadName=readName[i]+'circular:1:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup=eachreadName+'\\n'+circRead[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastqFile+'.write(\''+eachgroup+'\')')
            for p in range(len(linearRead)):
                eachreadName=readName[i]+'linear:1:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup=eachreadName+'\\n'+circRead[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastqFile+'.write(\''+eachgroup+'\')')
        else:
            fastq1File='sample_1_rep'+str(j+1)+'_R1'
            fastq2File='sample_1_rep'+str(j+1)+'_R2'
            circRead1,circRead2=simuCirc(circSeq,circN)
            linearRead1,linearRead2=simuLinear(linearSeq,tranN)
            for p in range(len(circRead1)):
                eachreadName=readName[i]+'circular:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup1=eachreadName+'/1\\n'+circRead1[p]+'\\n+\\n'+qualityLine+'\\n'
                eachgroup2=eachreadName+'/2\\n'+circRead2[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastq1File+'.write(\''+eachgroup1+'\')')
                eval(fastq2File+'.write(\''+eachgroup2+'\')')
            for p in range(len(linearRead1)):
                eachreadName=readName[i]+'linear:1:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup1=eachreadName+'/1\\n'+linearRead1[p]+'\\n+\\n'+qualityLine+'\\n'
                eachgroup2=eachreadName+'/2\\n'+linearRead2[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastq1File+'.write(\''+eachgroup1+'\')')
                eval(fastq2File+'.write(\''+eachgroup2+'\')')
    for j in range(Replicate_Num):      #sample_1
        tranN=transcriptNum.ix[transcriptName,j+Replicate_Num]
        circN=getBSnum(tranN,psi2mat.iloc[i,j])
        if SEPE=='SE':
            fastqFile='sample_2_rep'+str(j+1)
            circRead=simuCirc(circSeq,circN)
            linearRead=simuLinear(linearSeq,tranN)
            
            for p in range(len(circRead)):
                eachreadName=readName[i]+'circular:1:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup=eachreadName+'\\n'+circRead[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastqFile+'.write(\''+eachgroup+'\')')
            for p in range(len(linearRead)):
                eachreadName=readName[i]+'linear:1:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup=eachreadName+'\\n'+circRead[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastqFile+'.write(\''+eachgroup+'\')')
        else:
            fastq1File='sample_2_rep'+str(j+1)+'_R1'
            fastq2File='sample_2_rep'+str(j+1)+'_R2'
            #print(circSeq)
            #print(circN)
            circRead1,circRead2=simuCirc(circSeq,circN)
            linearRead1,linearRead2=simuLinear(linearSeq,tranN)
            for p in range(len(circRead1)):
                eachreadName=readName[i]+'circular:1:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup1=eachreadName+'/1\\n'+circRead1[p]+'\\n+\\n'+qualityLine+'\\n'
                eachgroup2=eachreadName+'/2\\n'+circRead2[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastq1File+'.write(\''+eachgroup1+'\')')
                eval(fastq2File+'.write(\''+eachgroup2+'\')')
            for p in range(len(linearRead1)):
                eachreadName=readName[i]+'linear:1:'+str(p)+':'+str(j+1)+':'+str(i)
                eachgroup1=eachreadName+'/1\\n'+linearRead1[p]+'\\n+\\n'+qualityLine+'\\n'
                eachgroup2=eachreadName+'/2\\n'+linearRead2[p]+'\\n+\\n'+qualityLine+'\\n'
                eval(fastq1File+'.write(\''+eachgroup1+'\')')
                eval(fastq2File+'.write(\''+eachgroup2+'\')')
