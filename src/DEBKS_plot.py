"""
author: Zelin Liu
email: zlliu@bjmu.edu.cn
license: GPL3
detail: Plot differentially expressed back-splicing base on the output of main script
"""
import sys
if sys.version_info[0] < 3:
    print ("Python Version error: must use python 3 (not supporting python 2)");
    sys.exit(-1);

import re,os,traceback,re,subprocess;
import pandas as pd
import numpy as np

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

    
### software required
#samtools
status,output=commandRun('samtools')
if status !=1:
    print ("Error: samtools is required");
    sys.exit(-1);
    
### input
s1=''; ## sample_1 bam name
s2=''; ## sample_2 bam name
BSfile='';
outDir='';

validArgList=['-s1','-s2','-o','-BS']
for argIndex in range(1,len(sys.argv)): ## going through the all parameters 
    if(sys.argv[argIndex][0]=='-' and sys.argv[argIndex] not in validArgList): ## incorrect argument
        print('Not valid argument: %s' % sys.argv[argIndex]);
        print('Please provide valid arguments.');
        sys.exit();

for paramIndex in range(1,len(sys.argv)): ## going through the all parameters
    if(sys.argv[paramIndex] == '-s1'):  
        paramIndex += 1;  ## increase index
        s1 = sys.argv[paramIndex];
    elif(sys.argv[paramIndex] == '-s2'): 
        paramIndex += 1;  ## increase index
        s2 = sys.argv[paramIndex];
    elif(sys.argv[paramIndex] == '-o'):  
        paramIndex += 1;  ## increase index
        outDir = sys.argv[paramIndex];
    elif (sys.argv[paramIndex] == '-BS'):
        paramIndex += 1;  ## increase index
        BSfile = sys.argv[paramIndex];

if (s1 == '' or  s2 ==''  or BSfile=='' or outDir==''):
    print('Not enough arguments!!');
    print('Usage (with bam files):\n\tpython DEBKS_plot.py -s1 s1File -s2 s2File -o outDir -BS DEBKS_Results');
    sys.exit(-1);

#### check file ####
# sample file
tempSample_1=open(s1,'r')
tempSample_2=open(s2,'r')
sample_1=[i.strip() for i in tempSample_1.readlines()]
sample_2=[i.strip() for i in tempSample_2.readlines()]
for i in sample_1+sample_2:
    if not os.path.exists(i):
        print ("sample file: "+ i +" is not exist");
        sys.exit();    

# DEBKS_Results file
tempBS_file = open(BSfile,'r'); ## open BSfile file
for line in tempBS_file:
  if line.strip()[0]=='#': ## comments, skip this line
    continue;
  Ele = line.strip().split('\t');
  if len(Ele) < 26: ## may be incorrect BSfile format
    print("Incorrect DEBKS_Results file format. File must have 26 tab-delimited columns.");
    sys.exit();
  break;  ## just check the first non-comment column
tempBS_file.close()
        
#### make output file ####
os.system('mkdir -p '+ outDir);

scriptPath = os.path.abspath(os.path.dirname(__file__));  ## absolute script path
outPath = os.path.abspath(outDir); ## absolute output path

s1Path = outPath + '/SAMPLE_1';
os.system('mkdir -p '+ s1Path);
s2Path = outPath + '/SAMPLE_2';
os.system('mkdir -p '+ s2Path);

## making folders for replicates ##
s1rPath = s1Path+'/REP_';
s2rPath = s2Path+'/REP_';
for rr in range(0,len(sample_1)): ## sample_1
    os.system('mkdir -p '+ s1rPath+str(rr+1));
for rr in range(0,len(sample_2)): ## sample_2
    os.system('mkdir -p '+ s2rPath+str(rr+1));

finalPath = outPath+'/DEBKS_plot';
os.system('mkdir -p '+ finalPath);

tempPath = finalPath + '/temp';
os.system('mkdir -p '+ tempPath);

sampleNum=len(sample_1)+len(sample_2)
########## functions here... ############

def doBAMsort():
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1);
        print('sort begin to file: %s' % sample_1[rr])
        if os.path.exists(rTempFolder+"/Aligned.out.sort.bam"):
            continue
        cmd = 'samtools sort -@ 10 -m 2G -o ' + rTempFolder + '/Aligned.out.sort.bam ' + sample_1[rr]
        status,output=commandRun(cmd)
        print('sort down to file: %s' % sample_1[rr])
        if (int(status)!=0): ## it did not go well
            print("error in sort sample_1, rep_%d: %s" % ((rr+1),status));
            print("error detail: %s" % output);
            raise Exception();
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1);
        print('sort begin to file: %s' % sample_2[rr])
        if os.path.exists(rTempFolder+"/Aligned.out.sort.bam"):
            continue
        cmd = 'samtools sort -@ 10 -m 2G -o ' + rTempFolder + '/Aligned.out.sort.bam ' + sample_2[rr]
        status,output=commandRun(cmd)
        print('sort down to file: %s' % sample_2[rr])
        if (int(status)!=0): ## it did not go well
            print("error in sort sample_2, rep_%d: %s" % ((rr+1),status));
            print("error detail: %s" % output);
            raise Exception();
    return;
    
##### end of doBAMsort ####
def buildIndex():
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1);
        print('index begin to file: %s' % sample_1[rr])
        if os.path.exists(rTempFolder+"/Aligned.out.sort.bam.bai"):
            continue
        cmd = 'samtools index ' + rTempFolder + '/Aligned.out.sort.bam'
        status,output=commandRun(cmd)
        print('index down to file: %s' % sample_1[rr])
        if (int(status)!=0): ## it did not go well
            print("error in index sample_1, rep_%d: %s" % ((rr+1),status));
            print("error detail: %s" % output);
            raise Exception();
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1);
        print('index begin to file: %s' % sample_2[rr])
        if os.path.exists(rTempFolder+"/Aligned.out.sort.bam.bai"):
            continue
        cmd = 'samtools index ' + rTempFolder + '/Aligned.out.sort.bam'
        status,output=commandRun(cmd)
        print('index down to file: %s' % sample_2[rr])
        if (int(status)!=0): ## it did not go well
            print("error in index sample_2, rep_%d: %s" % ((rr+1),status));
            print("error detail: %s" % output);
            raise Exception();
    return;
    
##### end of buildIndex ####
def getMapReads():
    mapRead1=[]
    mapRead2=[]
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1);
        print('flagstat begin to file: %s' % sample_1[rr])
        if os.path.exists(rTempFolder+"/flagstat.txt"):
            print('flagstat exists to file: %s' % sample_1[rr])
        else:
            cmd = 'samtools flagstat -@ 10 ' + rTempFolder + '/Aligned.out.sort.bam >' + rTempFolder + '/flagstat.txt'
            status,output=commandRun(cmd)
            print('flagstat down to file: %s' % sample_1[rr])
            if (int(status)!=0): ## it did not go well
                print("error in flagstat sample_1, rep_%d: %s" % ((rr+1),status));
                print("error detail: %s" % output);
                raise Exception();
        tmp=open(rTempFolder + '/flagstat.txt','r')
        mapRead1.append(re.search('(\d+) \+',tmp.readline()).group(1))
        tmp.close()
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1);
        print('flagstat begin to file: %s' % sample_2[rr])
        if os.path.exists(rTempFolder+"/flagstat.txt"):
            print('flagstat exists to file: %s' % sample_1[rr])
        else:
            cmd = 'samtools flagstat -@ 10 ' + rTempFolder + '/Aligned.out.sort.bam >' + rTempFolder + '/flagstat.txt'
            status,output=commandRun(cmd)
            print('flagstat down to file: %s' % sample_1[rr])
            if (int(status)!=0): ## it did not go well
                print("error in flagstat sample_2, rep_%d: %s" % ((rr+1),status));
                print("error detail: %s" % output);
                raise Exception();
        tmp=open(rTempFolder + '/flagstat.txt','r')
        mapRead2.append(re.search('(\d+) \+',tmp.readline()).group(1))
        tmp.close()
    return(','.join(mapRead1),','.join(mapRead2));

#### end of getMapReads ####
def getDepth(chr,start,end):
    targetRegion=str(chr)+':'+str(start)+'-'+str(end)
    matDepth=pd.DataFrame(0,index=range(end-start+1),columns=range(sampleNum))
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1);
        #print('depth begin to file: %s' % sample_1[rr])
        cmd='samtools depth -a -r '+targetRegion+' '+ rTempFolder+'/Aligned.out.sort.bam >'+tempPath+'/tmp.txt'
        status,output=commandRun(cmd)
        #print('depth down to file: %s' % sample_1[rr])
        if (int(status)!=0): ## it did not go well
            print("error in index sample_1, rep_%d: %s" % ((rr+1),status));
            print("error detail: %s" % output);
            raise Exception();
        matDepth.iloc[:,rr]=pd.read_table(tempPath+'/tmp.txt',header=None).iloc[:,2].tolist()
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1);
        #print('depth begin to file: %s' % sample_2[rr])
        cmd='samtools depth -a -r '+targetRegion+' '+ rTempFolder+'/Aligned.out.sort.bam >'+tempPath+'/tmp.txt'
        status,output=commandRun(cmd)
        #print('depth down to file: %s' % sample_1[rr])
        if (int(status)!=0): ## it did not go well
            print("error in index sample_2, rep_%d: %s" % ((rr+1),status));
            print("error detail: %s" % output);
            raise Exception();
        matDepth.iloc[:,(len(sample_2)+rr)]=pd.read_table(tempPath+'/tmp.txt',header=None).iloc[:,2].tolist()
    return(matDepth);
#### end of getDepth ####
######## end of functions ##############


################## actual process ##############
####
#### 1. BAM sort
####
print('BAM sort begin')
try:
    doBAMsort();
except:
    print("There is an exception in BAM sort");
    print("Exception: %s" % sys.exc_info()[0]);
    print("Detail: %s" % sys.exc_info()[1]);
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        print("Stack: %s" % stack);
    sys.exit(-1);
print("done BAM sort'..")

####
#### 2. BAM index
####
print('BAM index begin')
try:
    buildIndex();
except:
    print("There is an exception in BAM index");
    print("Exception: %s" % sys.exc_info()[0]);
    print("Detail: %s" % sys.exc_info()[1]);
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        print("Stack: %s" % stack);
    sys.exit(-1);
print("done BAM index'..")

####
#### 3. Flagstat
####
print('Flagstat begin')
try:
    mapRead1,mapRead2=getMapReads();
except:
    print("There is an exception in flagstat");
    print("Exception: %s" % sys.exc_info()[0]);
    print("Detail: %s" % sys.exc_info()[1]);
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        print("Stack: %s" % stack);
    sys.exit(-1);
print("done flagstat..")

mapReadAll=mapRead1+mapRead2
scriptDir=os.path.dirname(os.path.abspath(__file__))
####
#### 4. Plot BS
####
print('Plot BS begin')
matDEBKS=pd.read_table(BSfile,header=0)
for i in range(matDEBKS.iloc[:,0].size):
    chr=matDEBKS.iat[i,0]
    start=matDEBKS.iat[i,1]
    end=matDEBKS.iat[i,2]
    strand=matDEBKS.iat[i,3]
    exonSize=matDEBKS.iat[i,5]
    exonOff=matDEBKS.iat[i,6]
    geneName=matDEBKS.iat[i,7]
    isoName=matDEBKS.iat[i,8]
    Intron=matDEBKS.iat[i,9]
    baseL=re.match('.*?:(\d+)-.*?|.*',Intron)
    baseR=re.match('.*?\|.*?:.*?-(\d*)',Intron)
    if baseL is not None:
        baseL=int(baseL.group(1))
    else:
        baseL=0
    if baseR is not None:
        baseR=int(baseR.group(1))
    else:
        baseR=0
                
    linearL=matDEBKS.iat[i,10]
    linearR=matDEBKS.iat[i,11]
    SJL1=matDEBKS.iat[i,12]
    SJL2=matDEBKS.iat[i,13]
    SJR1=matDEBKS.iat[i,14]
    SJR2=matDEBKS.iat[i,15]
    BS1=matDEBKS.iat[i,18]
    BS2=matDEBKS.iat[i,19]
    psi1=matDEBKS.iat[i,22]
    psi2=matDEBKS.iat[i,23]
    fdr=matDEBKS.iat[i,25]
    if fdr>0.05:
        continue
    id=str(chr)+'_'+str(start)+'_'+str(end)
    if linearL == 0:
        plotL = start
    else:
        plotL = linearL
    if linearR ==0:
        plotR = end
    else:
        plotR = linearR
    matDepth=getDepth(chr,plotL,plotR)
    fileName=tempPath+'/'+id+'.depth.txt'
    matDepth.to_csv(fileName,sep='\t',index=False,header=False)
    outFile=finalPath+'/'+id+'.pdf'
    status,output=commandRun('Rscript '+scriptDir+'/Rplot.R %s' %' '.join([str(i) for i in [chr,start,end,strand,exonSize,exonOff,geneName,isoName,
    baseL,baseR,linearL,linearR,SJL1,SJL2,SJR1,SJR2,BS1,BS2,
    psi1,psi2,mapRead1,mapRead2,fileName,outFile]]))
    if (int(status)!=0): ## it did not go well
        print("Warrings in Rplot %s: %s" % (id,status))
        print("Warrings detail: %s" % output)
    print('down plot %s' % id)
    
