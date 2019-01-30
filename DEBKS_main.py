#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
author: Zelin Liu
email: zlliu@bjmu.edu.cn
license: GPL3
detail: Analysis differential back-spliced in of ribo-zero RNA-seq
"""
###DEBKS version
DEBKS_ver="1.0.0"
### help info
helpInfo='''
Usage (with FASTQ files):\n\tpython DEBKS_main.py -g genomeFasta -s1 s1File -s2 s2File  -STARindex STARindexDir -gtf gtfFile -o outDir -read readType -len readLength [options]*

Usage (with chimeric and spliced junction files):\n\tpython DEBKS_main.py -g genomeFasta -s1CJ s1CJFile -s2CJ s2CJFile -s1SJ s1SJFile -s2SJ s2SJFile -gtf gtfFile -o outDir -read readType -len readLength [options]*

DEBKS v1.0.0: Analysis of differential percent back-spliced in (C) Zelin Liu, 2019
###Required Parameters:
	-g          <str>       Genome Fasta file

	-STARindex  <str>       STAR alignment index directory

	-s1         <str>       Fastq files of sample 1, replicates in different lines, paired files are separated by semicolon

	-s2         <str>       Fastq files of sample 2, replicates in different lines, paired files are separated by semicolon

	-s1CJ       <str>       Chimeric junction of sample 1 group, replicates in different lines

	-s2CJ       <str>       Chimeric junction of sample 2 group, replicates in different lines

	-s1SJ       <str>       Spliced junction of sample 1 group, replicates in different lines

	-s2SJ       <str>       Spliced junction of sample 2 group, replicates in different lines

	-gtf        <str>       GTF file

	-o          <str>       Output directory of the result files

	-read       <str>       RNA-seq reads are single- or pair-end reads. [single, pair]

	-len        <int>       Read length of RNA-seq reads


###Optional Parameters:
	-h, --help              Show this help message and exit

	-t          <int>       Number of processors [1]

	-c          <float>     Required PBSI difference cutoff between the two samples [0.01]

	-a          <int>       Anchor length for counting chimeric or splicing junctions [6]

	-KeepTemp               Keep the temporary files. Disable by default.
'''
import sys
### python version
if sys.version_info[0] < 3:
    print ("Python version error: python 3 or above is required")
    sys.exit(-1)
    
### import libraries
import subprocess
import re,os,warnings,math,scipy,itertools,logging,time,datetime,pysam,traceback
from scipy import stats
from multiprocessing import Pool
from scipy.optimize import fmin_cobyla
from scipy.optimize import fmin_l_bfgs_b
from math import log
import pandas as pd
import numpy as np
from numpy import *
np.seterr(divide='ignore',invalid='ignore')

# system command
def commandRun(cmd):
    p= subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    status=p.wait()
    status=p.returncode
    outputALL=p.communicate()
    output=''
    for i in outputALL:
        output+=i.decode()
    return(status,output)

### software required
#STAR, gtfToGenePred, CIRCexplorer2
status,output=commandRun('STAR')
if status !=0:
    print ("Error: STAR is required")
    sys.exit(-1)
    
status,output=commandRun('gtfToGenePred')
if status !=255:
    print ("Error: gtfToGenePred is required")
    sys.exit(-1)
    
status,output=commandRun('CIRCexplorer2')
if status !=1:
    print ("Error: CIRCexplorer2 is required")
    sys.exit(-1)
#

#
##### parameter variables ####
### required values
gtf=''     ## gtf file 
outDir=''    ## output directory
readLength = 0 ## read length
readType = ''  ## read type, single or pair
threads = 1
genomeDir= ''
genomeFasta=''
#
### input
s1='' ## sample_1 fastq name
s2='' ## sample_2 fastq name
s1CJ='' ## sample_1 chimeric junction file name
s2CJ='' ## sample_2 chimeric junction file name
s1SJ='' ## sample_1 splicing junction file name
s2SJ='' ## sample_1 splicing junction file name

### with default values
Anchorlength=6 ## Chimeric anchor length
cutoff=0.01
keepTemp = 0 ## keep temp



### checking out the argument names
validArgList=['-g','-s1','-s2','-STARindex','-s1CJ','-s2CJ','-s1SJ','-s2SJ','-read','-t','-gtf','-o','-len','-c','-keepTemp','-a','-h','--help']
for argIndex in range(1,len(sys.argv)): ## going through the all parameters 
    if(sys.argv[argIndex][0]=='-' and sys.argv[argIndex] not in validArgList): ## incorrect argument
        print('Not valid argument: %s' % sys.argv[argIndex])
        print('Please provide valid arguments.')
        print(helpInfo)
        sys.exit()

for paramIndex in range(1,len(sys.argv)): ## going through the all parameters
    if(sys.argv[paramIndex] == '-g'):  ## sample_1
        paramIndex += 1  ## increase index
        genomeFasta = sys.argv[paramIndex]
    elif(sys.argv[paramIndex] == '-s1'):  ## sample_1
        paramIndex += 1  ## increase index
        s1 = sys.argv[paramIndex];    
    elif (sys.argv[paramIndex] == '-s2'):  ## sample_2
        paramIndex += 1  ## increase index
        s2 = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-STARindex'):  ## STAR genome
        paramIndex += 1  ## increase index
        genomeDir = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-s1CJ'):  ## sample_1 chimeric junction
        paramIndex += 1  ## increase index
        s1CJ = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-s2CJ'):  ## sample_2 chimeric junction
        paramIndex += 1  ## increase index
        s2CJ = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-s1SJ'):  ## sample_1 splicing junction
        paramIndex += 1  ## increase index
        s1SJ = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-s2SJ'):  ## sample_2 splicing junction
        paramIndex += 1  ## increase index
        s2SJ = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-gtf'):  ## gtf file
        paramIndex += 1  ## increase index
        gtf = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-o'):  ## output folder
        paramIndex += 1  ## increase index
        outDir = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-read'):  ## readtype, single or pair
        paramIndex += 1  ## increase index
        readType = sys.argv[paramIndex]
        if readType not in ['single','pair']:
            print('-read only accept "single" or "pair" value')
            sys.exit(1)
    elif (sys.argv[paramIndex] == '-t'):  ## threads for mapping
        paramIndex += 1  ## increase index
        threads = int(sys.argv[paramIndex])
    elif (sys.argv[paramIndex] == '-len'):  ## read length
        paramIndex += 1  ## increase index
        readLength = int(sys.argv[paramIndex])
    elif (sys.argv[paramIndex] == '-c'):  ##  c, deltaPBSI
        paramIndex += 1  ## increase index
        cutoff = float(sys.argv[paramIndex])
    elif (sys.argv[paramIndex] == '-a'):  ## anchor length for chimeric junction and linear junction
        paramIndex += 1  ## increase index
        Anchorlength = int(sys.argv[paramIndex])
    elif (sys.argv[paramIndex] == '-keepTemp'):  ## keep temp files, no value needed
        keepTemp = 1  ## keep temp files
    elif (sys.argv[paramIndex] == '-h' or sys.argv[paramIndex] == '--help'):  ## print helpInfo
        print(helpInfo)
        sys.exit()

### checking out the required arguments

if (genomeFasta == '' or ((s1=='' or  s2=='' or genomeDir == '') and (s1CJ=='' or  s2CJ=='' or s1SJ=='' or  s2SJ=='')) or gtf==''  or outDir=='' or readLength==0 or readType==''): ### at least one required param is missing
    print('Not enough arguments!!')
    print(helpInfo)
    sys.exit()

CAlength=Anchorlength
LAlength=Anchorlength
genomeFasta=os.path.abspath(genomeFasta)
gtf=os.path.abspath(gtf)

def listToString(x):
    rVal = ''
    for a in x:
        rVal += a+' '
    return rVal

### process input parameters ####
SEPE = 'SE' ## single-end or pair
if readType=="pair": ## single-end
  SEPE='PE'
##### checking Fasta format ####
tempFasta_file=open(genomeFasta,'r')
for line in tempFasta_file:
  if line.strip()[0]=='#': ## comments, skip this line
    continue
  if line.strip()[0]=='>':
    break
  else:
    print ("Incorrect Fasta file format. Genome Fasta file must be clarified by >.")
    sys.exit()
tempFasta_file.close()

##### checking GTF format #####
#
tempGTF_file = open(gtf,'r') ## open gtf file
for line in tempGTF_file:
  if line.strip()[0]=='#': ## comments, skip this line
    continue
  gtfEle = line.strip().split('\t')
  if len(gtfEle)<9: ## may be incorrect gtf format
    print("Incorrect GTF file format. GTF file must have 9 tab-delimited columns.")
    sys.exit()
  break  ## just check the first non-comment column
tempGTF_file.close()

####### checking readLength ##########
fastq=1
if s1=='' and  s2=='':
    fastq=0

def checkLineNum(file,Num):
    fff=open(file)
    a=fff.readline().strip().split('\t')
    fff.close()
    if len(a) == Num:
        return 1
    else:
        return 0
    
    
if fastq==1: ## fastq file was provided
    genomeDir=os.path.abspath(genomeDir)
    tempSample_1=open(s1,'r')
    tempSample_2=open(s2,'r')
    sample_1=[i.strip() for i in tempSample_1.readlines()]
    sample_2=[i.strip() for i in tempSample_2.readlines()]
    for i in sample_1+sample_2:
        tmp=i.split(';')
        for j in tmp:
            if not os.path.exists(j):
                print ("sample file: "+ j +" is not exist")
                sys.exit()
else: ## CJ and SJ files were provided
    tempSample_1CJ=open(s1CJ,'r')
    tempSample_2CJ=open(s2CJ,'r')
    tempSample_1SJ=open(s1SJ,'r')
    tempSample_2SJ=open(s2SJ,'r')
    sample_1_CJ=[i.strip() for i in tempSample_1CJ.readlines()]
    sample_2_CJ=[i.strip() for i in tempSample_2CJ.readlines()]
    sample_1_SJ=[i.strip() for i in tempSample_1SJ.readlines()]
    sample_2_SJ=[i.strip() for i in tempSample_2SJ.readlines()]
    for i in sample_1_CJ+sample_2_CJ+sample_1_SJ+sample_2_SJ:
        if not os.path.exists(i):
            print("sample file: "+ i +" is not exist")
            sys.exit()
    sample_1=sample_1_CJ
    sample_2=sample_2_CJ
    for CJfile in sample_1_CJ[0:2]: ## examine each file
        if checkLineNum(CJfile,14) == 0:
            print("Incorrect CJ file format. chimeric junction file must have 14 tab-delimited columns.")
            sys.exit()
            break
    for CJfile in sample_2_CJ[0:2]: ## examine each file
        if checkLineNum(CJfile,14) == 0:
            print("Incorrect CJ file format. chimeric junction file must have 14 tab-delimited columns.")
            sys.exit()
            break
    for SJfile in sample_1_SJ[0:2]: ## examine each file
        if checkLineNum(SJfile,9) == 0:
            print("Incorrect SJ file format. chimeric junction file must have 9 tab-delimited columns.")
            sys.exit()
            break
    for SJfile in sample_2_SJ[0:2]: ## examine each file
        if checkLineNum(SJfile,9) == 0:
            print("Incorrect SJ file format. chimeric junction file must have 9 tab-delimited columns.")
            sys.exit()
            break
            
sampleNum=len(sample_1)+len(sample_2)

os.system('mkdir -p '+ outDir)
oFile = open(outDir+'/commands.txt', 'a') ## file that will contain list of commands excuted here

### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outDir+'/DEBKS_main.log' ,
                    filemode='w')

##### Getting Start Time ######
logging.debug('DEBKS version: %s' % DEBKS_ver)
logging.debug('Start the program with [%s]\n', listToString(sys.argv))
startTime = time.time()

scriptPath = os.path.abspath(os.path.dirname(__file__))  ## absolute script path
outPath = os.path.abspath(outDir) ## absolute output path

s1Path = outPath + '/SAMPLE_1'
os.system('mkdir -p '+ s1Path)
s2Path = outPath + '/SAMPLE_2'
os.system('mkdir -p '+ s2Path)

## making folders for replicates ##
s1rPath = s1Path+'/REP_'
s2rPath = s2Path+'/REP_'
for rr in range(0,len(sample_1)): ## sample_1
    os.system('mkdir -p '+ s1rPath+str(rr+1))
for rr in range(0,len(sample_2)): ## sample_2
    os.system('mkdir -p '+ s2rPath+str(rr+1))

finalPath = outPath+'/DEBKS_output'
os.system('mkdir -p '+ finalPath)

tempPath = outPath + '/temp'
os.system('mkdir -p '+ tempPath)

### putting keys in log file
#
logging.debug("################### folder names and associated input files #############")
for fki in range(0,len(sample_1)): ## for each replicate of sample_1
    repTempFolder = "SAMPLE_1\REP_"+str(fki+1)
    associatedFile = sample_1[fki]
    logging.debug(repTempFolder+"\t"+associatedFile)

for fki in range(0,len(sample_2)): ## for each replicate of sample_2
    repTempFolder = "SAMPLE_2\REP_"+str(fki+1)
    associatedFile = sample_2[fki]
    logging.debug(repTempFolder+"\t"+associatedFile)

logging.debug("#########################################################################\n")

#### ln -s chimeric and splicing junction files for no fastq ####
if fastq ==0:
    logging.debug("################# ln -s chimeric junction files for no fastq #############")
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1)
        if os.path.exists(rTempFolder+'/Chimeric.out.junction'):
            continue
        else:
            cmd='ln -s '+os.path.abspath(sample_1_CJ[rr]) +' '+rTempFolder+'/Chimeric.out.junction'
            status,output=commandRun(cmd)
            oFile.write('######  ln -s chimeric junction files for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n')
            oFile.flush()
            if (int(status)!=0): ## it did not go well
                logging.debug("error in ln -s chimeric junction files for sample_1, replicate_%d: %s" % ((rr+1),status))
                logging.debug("error detail: %s" % output)
                raise Exception()
        if os.path.exists(rTempFolder+'/SJ.out.tab'):
            continue
        else:
            cmd='ln -s '+os.path.abspath(sample_1_SJ[rr]) +' '+rTempFolder+'/SJ.out.tab'
            status,output=commandRun(cmd)
            oFile.write('######  ln -s splicing junction files for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n')
            oFile.flush()
            if (int(status)!=0): ## it did not go well
                logging.debug("error in ln -s splicing junction files for sample_1, replicate_%d: %s" % ((rr+1),status))
                logging.debug("error detail: %s" % output)
                raise Exception()

    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1)
        if os.path.exists(rTempFolder+'/Chimeric.out.junction'):
            continue
        else:
            cmd='ln -s '+os.path.abspath(sample_2_CJ[rr]) +' '+rTempFolder+'/Chimeric.out.junction'
            status,output=commandRun(cmd)
            oFile.write('######  ln -s chimeric junction files for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n')
            oFile.flush()
            if (int(status)!=0): ## it did not go well
                logging.debug("error in ln -s chimeric junction files for sample_2, replicate_%d: %s" % ((rr+1),status))
                logging.debug("error detail: %s" % output)
                raise Exception()
        if os.path.exists(rTempFolder+'/SJ.out.tab'):
            continue
        else:
            cmd='ln -s '+os.path.abspath(sample_2_SJ[rr]) +' '+rTempFolder+'/SJ.out.tab'
            status,output=commandRun(cmd)
            oFile.write('######  ln -s splicing junction files for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n')
            oFile.flush()
            if (int(status)!=0): ## it did not go well
                logging.debug("error in ln -s splicing junction files for sample_2, replicate_%d: %s" % ((rr+1),status))
                logging.debug("error detail: %s" % output)
                raise Exception()
    logging.debug("#########################################################################\n")
else:
    pass
########## functions here... ############

def doSTARMapping(): ## do STAR mapping
    logging.debug("mapping the sample1")
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1)
        isgz=0
        if os.path.exists(rTempFolder+"/Chimeric.out.junction") and os.path.exists(rTempFolder+"/SJ.out.tab"):
            continue
        else:
            if os.path.exists(rTempFolder+"/_STARtmp"):
                os.system("rm -rf "+rTempFolder+"/_STARtmp")
                unLoadGenome()
            if os.path.exists(rTempFolder+"/Aligned.out.bam"):
                os.system("rm -rf "+rTempFolder+"/Aligned.out.bam")
        cmd = 'STAR --genomeLoad LoadAndKeep --chimSegmentMin ' + str(CAlength) + ' --outFilterMismatchNmax 3  --runThreadN ' + str(threads) + ' --outSAMstrandField intronMotif --outSAMtype BAM Unsorted'
        cmd += ' --alignSJDBoverhangMin ' + str(LAlength) +' --alignSJoverhangMin '+ str(LAlength) + ' --alignIntronMax 300000 --genomeDir '+genomeDir
        cmd += ' --chimJunctionOverhangMin '+ str(CAlength) +' --outSJfilterOverhangMin -1 '+ str(LAlength) +' -1 -1 --outFileNamePrefix ' + rTempFolder + '/ --readFilesIn '
        if SEPE=='PE': ## pair-end
            cmd += sample_1[rr].split(';')[0]+' '+sample_1[rr].split(';')[1]
            if re.search('.gz$',sample_1[rr].split(';')[0]) is not None:
                isgz=1
        else: ## single-end
            cmd += sample_1[rr]
            if re.search('.gz$',sample_1[rr]) is not None:
                isgz=1
        if isgz == 1:
            cmd += ' --readFilesCommand zcat '
        oFile.write('######  running STAR for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commandRun(cmd)
        logging.debug("mapping sample_1, rep_"+str(rr+1)+" is done with status %s" % status)
        if (int(status)!=0): ## it did not go well
            logging.debug("error in mapping sample_1, rep_%d: %s" % ((rr+1),status))
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)
    logging.debug("mapping the sample2")
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1)
        isgz=0
        if os.path.exists(rTempFolder+"/Chimeric.out.junction") and os.path.exists(rTempFolder+"/SJ.out.tab"):
            continue
        else:
            if os.path.exists(rTempFolder+"/_STARtmp"):
                os.system("rm -rf "+rTempFolder+"/_STARtmp")
                unLoadGenome()
            if os.path.exists(rTempFolder+"/Aligned.out.bam"):
                os.system("rm -rf "+rTempFolder+"/Aligned.out.bam")
        cmd = 'STAR --genomeLoad LoadAndKeep --chimSegmentMin ' + str(CAlength) + '  --runThreadN ' + str(threads) + ' --outSAMstrandField intronMotif --outSAMtype BAM Unsorted'
        cmd += ' --alignSJDBoverhangMin ' + str(LAlength) +' --alignSJoverhangMin '+ str(LAlength) + ' --alignIntronMax 300000 --genomeDir '+genomeDir
        cmd += ' --chimJunctionOverhangMin '+ str(CAlength) +' --outSJfilterOverhangMin -1 '+ str(LAlength) +' -1 -1 --outFileNamePrefix ' + rTempFolder + '/ --readFilesIn '
        if SEPE=='PE': ## pair-end
            cmd += sample_2[rr].split(';')[0]+' '+sample_2[rr].split(';')[1]
            if re.search('.gz$',sample_1[rr].split(';')[0]) is not None:
                isgz=1
        else: ## single-end
            cmd += sample_2[rr]
            if re.search('.gz$',sample_1[rr]) is not None:
                isgz=1
        if isgz == 1:
            cmd += ' --readFilesCommand zcat '
        oFile.write('######  running STAR for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commandRun(cmd)
        logging.debug("mapping sample_2, rep_"+str(rr+1)+" is done with status %s" % status)
        if (int(status)!=0): ## it did not go well
            logging.debug("error in mapping sample_2, rep_%d: %s" % ((rr+1),status))
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)
        unLoadGenome()
    os.chdir(tempPath)
    
    return
##### end of doSTARMapping ####
def unLoadGenome():
    cmd='STAR --genomeLoad Remove  --genomeDir '+genomeDir
    oFile.write('######  running STAR for remove load genome  #####\n'+cmd+'\n#\n')
    oFile.flush()
    status,output=commandRun(cmd)
    logging.debug("Remove STAR load genome %s" % status)
    if (int(status)!=0): ## it did not go well
        logging.debug("Warring in remove load genome %s" % status)
        logging.debug("Warring detail: %s" % output)
        return
    logging.debug(output)
##### end of unLoadGenome ####

def doCIRCexplorer2():
    if os.path.exists(tempPath+'/All.Chimeric.out.junction'):
        logging.debug(tempPath+'/All.Chimeric.out.junction has already exist')
        return()
    logging.debug("merge all chimeric junctions")
    cmd='cat'
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1)
        cmd+=' '+rTempFolder+'/Chimeric.out.junction'
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1)
        cmd+=' '+rTempFolder+'/Chimeric.out.junction'
    cmd+=' >'+tempPath+'/All.Chimeric.out.junction'
    oFile.write('######  cat all samples chimeric junction files  #####\n'+cmd+'\n#\n')
    oFile.flush()
    status,output=commandRun(cmd)
    logging.debug("cat all Chimeric.out junction is down %s" % status)
    if (int(status)!=0): ## it did not go well
        logging.debug("error in cat all Chimeric.out junction %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)
    logging.debug("run CIRCexplorer2 parse to all chimeric juncions")
    cmd='CIRCexplorer2 parse -t STAR '+ tempPath+'/All.Chimeric.out.junction'
    oFile.write('######  Detect circRNA  #####\n'+cmd+'\n#\n')
    os.chdir(tempPath)
    oFile.flush()
    status,output=commandRun(cmd)
    logging.debug("CIRCexplorer2 parse is down %s" % status)
    if (int(status)!=0): ## it did not go well
        logging.debug("error in CIRCexplorer2 parse %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

###### end of doCIRCexplorer2 ####


def annotateCirc():
    if os.path.exists(tempPath+'/circularRNA_known.txt'):
        logging.debug(tempPath+'/circularRNA_known.txt has already exist')
        return()
    logging.debug("Annotate circRNA")
    logging.debug("gtfToGenePred begin")
    if os.path.exists(tempPath+'/GenePred.txt'):
        logging.debug("gtfToGenePred has already existed")
    else:
        cmd = 'gtfToGenePred -genePredExt -ignoreGroupsWithoutExons '+ gtf + ' ' +tempPath+'/GenePred.txt'
        oFile.write('######  gtfToGenePred  #####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commandRun(cmd)
        logging.debug("gtfToGenePred is down %s" % status)
        if (int(status)!=0): ## it did not go well
            logging.debug("error in  gtfToGenePred%s" % status)
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)
    logging.debug("Make genePred format fit CIRCexplorer2 format")
    cmd=' awk -F \'\\t\' \'{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}\' '+ tempPath+'/GenePred.txt >'+ tempPath+'/GenePred.CIRC.txt'
    oFile.write('######  GenePredtoCIRCGP  #####\n'+cmd+'\n#\n')
    oFile.flush()
    status,output=commandRun(cmd)
    logging.debug("GenePredtoCIRCGP is down %s" % status)
    if (int(status)!=0): ## it did not go well
        logging.debug("error in  GenePredtoCIRCGP%s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)
    logging.debug("CIRCexplorer2 annotate begin")
    cmd='CIRCexplorer2 annotate -r '+tempPath+'/GenePred.CIRC.txt'+ ' -g '+ genomeFasta + ' -b '+tempPath+'/back_spliced_junction.bed -o '+tempPath+'/circularRNA_known.txt' 
    oFile.write('######  CIRCexplorer2 annotate   #####\n'+cmd+'\n#\n')
    oFile.flush()
    status,output=commandRun(cmd)
    logging.debug("CIRCexplorer2 annotate is down %s" % status)
    if (int(status)!=0): ## it did not go well
        logging.debug("error in  CIRCexplorer2 annotate%s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)
    
###### end of annotateCirc####

def doCIRCexplorer2foreach():
    logging.debug("run CIRCexplorer2 parse to each chimeric juncions")
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1)
        if os.path.exists(rTempFolder+'/circularRNA_known.txt'):
            logging.debug(rTempFolder+'/circularRNA_known.txt has already existed')
            continue
        cmd='CIRCexplorer2 parse -t STAR '+ rTempFolder+'/Chimeric.out.junction'
        os.chdir(rTempFolder)
        oFile.write('######  Parse circRNA  #####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commandRun(cmd)
        logging.debug("CIRCexplorer2 parse is down %s" % status)
        if (int(status)!=0): ## it did not go well
            logging.debug("error in CIRCexplorer2 parse %s" % status)
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)
        cmd='CIRCexplorer2 annotate -r '+tempPath+'/GenePred.CIRC.txt'+ ' -g '+ genomeFasta + ' -b '+rTempFolder+'/back_spliced_junction.bed -o '+rTempFolder+'/circularRNA_known.txt'
        oFile.write('######  Annotate circRNA  #####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commandRun(cmd)
        logging.debug("CIRCexplorer2 annotate is down %s" % status)
        if (int(status)!=0): ## it did not go well
            logging.debug("error in CIRCexplorer2 parse %s" % status)
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)
        status,output=commandRun('sed -i \'s/circular_RNA\///g\' ' + rTempFolder+'/circularRNA_known.txt')
    for rr in range(0,len(sample_2)): ## for each replicate of sample_1
        rTempFolder = s2rPath+str(rr+1)
        if os.path.exists(rTempFolder+'/circularRNA_known.txt'):
            logging.debug(rTempFolder+'/circularRNA_known.txt has already existed')
            continue
        cmd='CIRCexplorer2 parse -t STAR '+ rTempFolder+'/Chimeric.out.junction'
        os.chdir(rTempFolder)
        oFile.write('######  Detect circRNA  #####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commandRun(cmd)
        logging.debug("CIRCexplorer2 parse is down %s" % status)
        if (int(status)!=0): ## it did not go well
            logging.debug("error in CIRCexplorer2 parse %s" % status)
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)
        cmd='CIRCexplorer2 annotate -r '+tempPath+'/GenePred.CIRC.txt'+ ' -g '+ genomeFasta + ' -b '+rTempFolder+'/back_spliced_junction.bed -o '+rTempFolder+'/circularRNA_known.txt'
        oFile.write('######  Annotate circRNA  #####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commandRun(cmd)
        logging.debug("CIRCexplorer2 annotate is down %s" % status)
        if (int(status)!=0): ## it did not go well
            logging.debug("error in CIRCexplorer2 parse %s" % status)
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)
        status,output=commandRun('sed -i \'s/circular_RNA\///g\' ' + rTempFolder+'/circularRNA_known.txt')
###### end of doCIRCexplorer2foreach####

def getCircRNAdetail():
    logging.debug("Get circRNA detail")
    status,output=commandRun('sed -i \'s/circular_RNA\///g\' ' + tempPath+'/circularRNA_known.txt')
    os.system('awk -F "\t" \'{if($4>='+str(2*(len(sample_1)+len(sample_2)))+'){print $0}}\' '+ tempPath+'/circularRNA_known.txt > '+tempPath+'/circularRNA_filter.txt')
    annoBS=pd.read_table(tempPath+'/circularRNA_filter.txt',header=None)
    annoBS=annoBS.assign(exonLen=0)
    title=['chr','start','end','name','score','strand','thickStart','thickEnd','itemRgb','exonCount','exonSizes','exonOffsets','readNumber','circType','geneName','isoformName','index','flankIntron']
    titledict=dict()
    for i in range(annoBS.columns.size-1):
        titledict[i]=title[i]
    rowIndex= dict(zip([ i  for i in range(annoBS.iloc[:,0].size)],annoBS[0]+'|'+annoBS[1].map(str)+'|'+annoBS[2].map(str)))
    annoBS.rename(index=rowIndex,columns=titledict,inplace=True)
    annoBS=annoBS.iloc[[not i for i in annoBS.index.duplicated()],:]
    exonLen=[]
    annoBS.iloc[:,18]=annoBS.iloc[:,10].map(lambda x: sum([int(i) for i in x.split(',')]))
    annoBS.to_csv(tempPath+'/annoBS.txt',index=True,sep='\t',header=True)
    matBS=pd.DataFrame(0,index=range(len(annoBS.index)),columns=range(sampleNum))
    matBS.rename(index=dict(zip([ i  for i in range(annoBS.iloc[:,0].size)],annoBS.index.tolist())),inplace=True)
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1)
        BSJbed=pd.read_table(rTempFolder+'/circularRNA_known.txt',header=None)
        BSrowIndex= dict(zip([ i  for i in range(BSJbed.iloc[:,0].size)],BSJbed[0]+'|'+BSJbed[1].map(str)+'|'+BSJbed[2].map(str)))
        BSJbed.rename(index=BSrowIndex,inplace=True)
        BSJbed=BSJbed.iloc[[not i for i in BSJbed.index.duplicated()],:]
        BSJbed=BSJbed.reindex(annoBS.index,axis=0)
        matBS.iloc[:,rr]=BSJbed.iloc[:,3]
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
        rTempFolder = s2rPath+str(rr+1)
        BSJbed=pd.read_table(rTempFolder+'/circularRNA_known.txt',header=None)
        BSrowIndex= dict(zip([ i  for i in range(BSJbed.iloc[:,0].size)],BSJbed[0]+'|'+BSJbed[1].map(str)+'|'+BSJbed[2].map(str)))
        BSJbed.rename(index=BSrowIndex,inplace=True)
        BSJbed=BSJbed.iloc[[not i for i in BSJbed.index.duplicated()],:]
        BSJbed=BSJbed.reindex(annoBS.index,axis=0)
        matBS.iloc[:,rr+len(sample_1)]=BSJbed.iloc[:,3]
    matBS.fillna(0,inplace=True)
    matBS.to_csv(tempPath+'/matBS.txt',index=True,sep='\t',header=False)
    return(matBS, annoBS)
##### end of getCircRNAdetail ####


def combineLinearJunction():
    logging.debug("Get splicing junction")
    BSrowIndex=matBS.index
    matSJ=pd.DataFrame('',index=range(len(BSrowIndex)),columns=range(sampleNum))
    matSJ_L=pd.DataFrame('',index=range(len(BSrowIndex)),columns=range(sampleNum))
    matSJ_R=pd.DataFrame('',index=range(len(BSrowIndex)),columns=range(sampleNum))
    matSJ.rename(index=dict(zip([ i  for i in range(annoBS.iloc[:,0].size)],annoBS.iloc[:,0]+'|'+annoBS.iloc[:,1].map(str)+'|'+annoBS.iloc[:,2].map(str))),inplace=True)
    annoSJ=pd.DataFrame('',index=range(len(BSrowIndex)),columns=range(sampleNum*4))
    annoSJ.rename(index=dict(zip([ i  for i in range(annoBS.iloc[:,0].size)],annoBS.iloc[:,0]+'|'+annoBS.iloc[:,1].map(str)+'|'+annoBS.iloc[:,2].map(str))),inplace=True)
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
        rTempFolder = s1rPath+str(rr+1)
        SJtab=pd.read_table(rTempFolder+'/SJ.out.tab',header=None)
        SJtab['startKey']=SJtab[0]+"|"+(SJtab[1]-1).map(str)
        SJtab['endKey']=SJtab[0]+"|"+SJtab[2].map(str)
        tempMergeL=matBS.merge(SJtab,left_on="startKey",right_on="endKey")
        tempMergeL.iloc[:,sampleNum+9]=(tempMergeL.iloc[:,sampleNum+9]+tempMergeL.iloc[:,sampleNum+10]).map(str)  #also count multiple mapping reads on junctions
        tempMergeR=matBS.merge(SJtab,left_on="endKey",right_on="startKey")
        tempMergeR.iloc[:,sampleNum+9]=(tempMergeR.iloc[:,sampleNum+9]+tempMergeR.iloc[:,sampleNum+10]).map(str)  #also count multiple mapping reads on junctions
        preKeyL=tempMergeL.at[0,'key']
        preL=str(tempMergeL.at[0,'1_y'])
        preCountL=tempMergeL.iat[0,sampleNum+9]
        for i in range(1,tempMergeL.iloc[:,0].size):
            nowKeyL=tempMergeL.at[i,'key']
            nowL=str(tempMergeL.at[i,'1_y'])
            nowCountL=tempMergeL.iat[i,sampleNum+9]
            if nowKeyL == preKeyL:
                preL+=','+nowL
                preCountL+=','+nowCountL
            else:
                annoSJ.at[preKeyL,4*rr]=preL
                annoSJ.at[preKeyL,4*rr+1]=preCountL
                preKeyL=nowKeyL
                preL=nowL
                preCountL=nowCountL
        annoSJ.at[preKeyL,4*rr]=preL
        annoSJ.at[preKeyL,4*rr+1]=preCountL
        preKeyR=tempMergeR.at[0,'key']
        preR=str(tempMergeR.at[0,'2_y'])
        preCountR=tempMergeR.iat[0,sampleNum+9]
        for i in range(1,tempMergeR.iloc[:,0].size):
            nowKeyR=tempMergeR.at[i,'key']
            nowR=str(tempMergeR.at[i,'2_y'])
            nowCountR=tempMergeR.iat[i,sampleNum+9]
            if nowKeyR == preKeyR:
                preR+=','+nowR
                preCountR+=','+nowCountR
            else:
                annoSJ.at[preKeyR,4*rr+2]=preR
                annoSJ.at[preKeyR,4*rr+3]=preCountR
                preKeyR=nowKeyR
                preR=nowR
                preCountR=nowCountR
        annoSJ.at[preKeyR,4*rr+2]=preR
        annoSJ.at[preKeyR,4*rr+3]=preCountR
    for rr in range(len(sample_1),sampleNum): ## for each replicate of sample_1
        rTempFolder = s2rPath+str(rr-len(sample_1)+1)
        SJtab=pd.read_table(rTempFolder+'/SJ.out.tab',header=None)
        SJtab['startKey']=SJtab[0]+"|"+(SJtab[1]-1).map(str)
        SJtab['endKey']=SJtab[0]+"|"+SJtab[2].map(str)
        tempMergeL=matBS.merge(SJtab,left_on="startKey",right_on="endKey")
        tempMergeL.iloc[:,sampleNum+9]=(tempMergeL.iloc[:,sampleNum+9]+tempMergeL.iloc[:,sampleNum+10]).map(str)   #also count multiple mapping reads on junctions
        tempMergeR=matBS.merge(SJtab,left_on="endKey",right_on="startKey")
        tempMergeR.iloc[:,sampleNum+9]=(tempMergeR.iloc[:,sampleNum+9]+tempMergeR.iloc[:,sampleNum+10]).map(str)   #also count multiple mapping reads on junctions
        preKeyL=tempMergeL.at[0,'key']
        preL=str(tempMergeL.at[0,'1_y'])
        preCountL=tempMergeL.iat[0,sampleNum+9]
        for i in range(1,tempMergeL.iloc[:,0].size):
            nowKeyL=tempMergeL.at[i,'key']
            nowL=str(tempMergeL.at[i,'1_y'])
            nowCountL=tempMergeL.iat[i,sampleNum+9]
            if nowKeyL == preKeyL:
                preL+=','+nowL
                preCountL+=','+nowCountL
            else:
                annoSJ.at[preKeyL,4*rr]=preL
                annoSJ.at[preKeyL,4*rr+1]=preCountL
                preKeyL=nowKeyL
                preL=nowL
                preCountL=nowCountL
        annoSJ.at[preKeyL,4*rr]=preL
        annoSJ.at[preKeyL,4*rr+1]=preCountL
        preKeyR=tempMergeR.at[0,'key']
        preR=str(tempMergeR.at[0,'2_y'])
        preCountR=tempMergeR.iat[0,sampleNum+9]
        for i in range(1,tempMergeR.iloc[:,0].size):
            nowKeyR=tempMergeR.at[i,'key']
            nowR=str(tempMergeR.at[i,'2_y'])
            nowCountR=tempMergeR.iat[i,sampleNum+9]
            if nowKeyR == preKeyR:
                preR+=','+nowR
                preCountR+=','+nowCountR
            else:
                annoSJ.at[preKeyR,4*rr+2]=preR
                annoSJ.at[preKeyR,4*rr+3]=preCountR
                preKeyR=nowKeyR
                preR=nowR
                preCountR=nowCountR
        annoSJ.at[preKeyR,4*rr+2]=preR
        annoSJ.at[preKeyR,4*rr+3]=preCountR
    annoSJ.to_csv(tempPath+'/annoSJ.txt',index=True,sep='\t',header=False)
    annoBS.iloc[:,1]=annoBS.iloc[:,1].map(str)
    annoBS.iloc[:,2]=annoBS.iloc[:,2].map(str)
    for i in range(annoSJ.iloc[:,0].size):
        Intron=annoBS.iat[i,17]
        exonS=annoBS.iat[i,1]
        exonE=annoBS.iat[i,2]
        baseL=re.match('.*?:(\d+)-'+exonS+'.*',Intron)
        baseR=re.match('.*?:'+exonE+'-(\d*)',Intron)
        for j in range(sampleNum):
            tmpl=0
            tmpr=0
            
            if baseL is not None:
                tmp1=annoSJ.iat[i,4*j]
                if len(tmp1)>0:
                    vecL=[int(i) for i in tmp1.split(',')]
                    if int(baseL.group(1))+1 in vecL:
                        idx=vecL.index(int(baseL.group(1))+1)
                        tmpl=int(annoSJ.iat[i,4*j+1].split(',')[idx])
            if baseR is not None:
                tmp2=annoSJ.iat[i,4*j+2]
                if len(tmp2) > 0:
                    vecR=[int(i) for i in tmp2.split(',')]
                    if int(baseR.group(1)) in vecR:
                        idx=vecR.index(int(baseR.group(1)))
                        tmpr=int(annoSJ.iat[i,4*j+3].split(',')[idx])
            matSJ.iat[i,j]=tmpl+tmpr
            matSJ_L.iat[i,j]=tmpl
            matSJ_R.iat[i,j]=tmpr
    list1_L=matSJ_L.iloc[:,0].map(str)
    list2_L=matSJ_L.iloc[:,len(sample_1)].map(str)
    list1_R=matSJ_R.iloc[:,0].map(str)
    list2_R=matSJ_R.iloc[:,len(sample_1)].map(str)
    for i in range(1,len(sample_1)):
        list1_L=list1_L+','+matSJ_L.iloc[:,i].map(str)
        list1_R=list1_R+','+matSJ_R.iloc[:,i].map(str)
    for i in range(len(sample_1)+1,sampleNum):
        list2_L=list2_L+','+matSJ_L.iloc[:,i].map(str)
        list2_R=list2_R+','+matSJ_R.iloc[:,i].map(str)
    matSJ.to_csv(tempPath+'/matSJ.txt',index=True,sep='\t',header=False)
    return (matSJ,annoSJ,list1_L.tolist(),list2_L.tolist(),list1_R.tolist(),list2_R.tolist())

###### end of combineLinearJunction####


def getannoTrans():
    annoFile=pd.read_table(tempPath+'/GenePred.CIRC.txt',header=None)
    annoFile.rename(index=dict(zip(annoFile.index,annoFile.iloc[:,1])),columns=dict(zip([i for i in range(11)],['gene','transcript','chr','strand','start','end','cdsS','cdsE','exonNum','exonS','exonE'])),inplace=True)
    annoTrans=pd.DataFrame(0,index=range(annoBS.iloc[:,0].size),columns=range(4))
    for i in range(annoBS.iloc[:,0].size):
        isoform=annoBS.iat[i,15]
        isoExonS=annoFile.at[isoform,'exonS']
        isoExonE=annoFile.at[isoform,'exonE']
        if len(isoExonS)>1:
            isoExonS=[int(i) for i in isoExonS[0:(len(isoExonS)-1)].split(',')]
            isoExonE=[int(i) for i in isoExonE[0:(len(isoExonE)-1)].split(',')]
        else:
            continue
        Intron=annoBS.iat[i,17]
        exonS=annoBS.iat[i,1]
        exonE=annoBS.iat[i,2]
        baseL=re.match('.*?:(\d+)-'+exonS+'.*',Intron)
        baseR=re.match('.*?:'+exonE+'-(\d*)',Intron)
        if baseL is not None:
            baseL=int(baseL.group(1))
            annoTrans.iat[i,1]=baseL
            if baseL in isoExonE:
                annoTrans.iat[i,0]=isoExonS[isoExonE.index(baseL)]
        if baseR is not None:
            baseR=int(baseR.group(1))
            annoTrans.iat[i,2]=baseR
            if baseR in isoExonS:
                annoTrans.iat[i,3]=isoExonE[isoExonS.index(baseR)]
    annoTrans.to_csv(tempPath+'/annoTrans.txt',index=False,sep='\t',header=False)
    return(annoTrans)



##### end of getannoTrans #####
#binomial MLE optimization functions
def logit(x):
    if x<0.01:
        x=0.01
    if x>0.99:
        x=0.99
    return(log(x/(1-x)))

def myfunc_multivar(x,*args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3]
    #print('psi1');print(psi1);print('psi2');print(psi2)
    sum1=0;sum2=0
    for i in psi1:
        sum1+=pow(logit(i)-logit(x[0]),2)
    sum1=sum1/var1/2
    for i in psi2:
        sum2+=pow(logit(i)-logit(x[1]),2)
    sum2=sum2/var2/2
    return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x[0]),2)+pow(stats.norm.ppf(x[1]),2)-2*rho*stats.norm.ppf(x[0])*stats.norm.ppf(x[1])))

def myfunc_multivar_der(x,*args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3]
    sum1=0;sum2=0
    for i in psi1:
        sum1+=-2*(logit(i)-logit(x[0]))/x[0]/(1-x[0])
    sum1=sum1/var1/2
    res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[0])-2*rho*stats.norm.ppf(x[1]))/stats.norm.pdf(stats.norm.ppf(x[0]))
    #print('1');print(x[1]);print(res1);print(stats.norm.pdf(stats.norm.ppf(x[0])))
    for i in psi2:
        sum2+=-2*(logit(i)-logit(x[1]))/x[1]/(1-x[1])
    sum2=sum2/var2/2
    res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[1])-2*rho*stats.norm.ppf(x[0]))/stats.norm.pdf(stats.norm.ppf(x[1]))
    return(np.array([res1,res2]))

def myfunc_1(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3]
    sum1=0;sum2=0
    for i in psi1:
        sum1+=pow(logit(i)-logit(x+cutoff),2)
    sum1=sum1/var1/2
    for i in psi2:
        sum2+=pow(logit(i)-logit(x),2)
    sum2=sum2/var2/2
    return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x+cutoff),2)+pow(stats.norm.ppf(x),2)-2*rho*stats.norm.ppf(x+cutoff)*stats.norm.ppf(x)))

def myfunc_der_1(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3]
    sum1=0;sum2=0
    for i in psi1:
        sum1+=-2*(logit(i)-logit(x+cutoff))/(x+cutoff)/(1-x-cutoff)
    sum1=sum1/var1/2
    res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x+cutoff)-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf(x+cutoff))
    for i in psi2:
        sum2+=-2*(logit(i)-logit(x))/x/(1-x)
    sum2=sum2/var2/2
    res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf(x+cutoff))/stats.norm.pdf(stats.norm.ppf(x))
    return(res1+res2)

def myfunc_2(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3]
    sum1=0;sum2=0
    for i in psi1:
        sum1+=pow(logit(i)-logit(x),2)
    sum1=sum1/var1/2
    for i in psi2:
        sum2+=pow(logit(i)-logit(x+cutoff),2)
    sum2=sum2/var2/2
    return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x+cutoff),2)+pow(stats.norm.ppf(x),2)-2*rho*stats.norm.ppf(x+cutoff)*stats.norm.ppf(x)))

def myfunc_der_2(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3]
    sum1=0;sum2=0
    for i in psi1:
        sum1+=-2*(logit(i)-logit(x))/(x)/(1-x)
    sum1=sum1/var1/2
    res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x+cutoff)-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf(x+cutoff))
    for i in psi2:
        sum2+=-2*(logit(i)-logit(x+cutoff))/(x+cutoff)/(1-x-cutoff)
    sum2=sum2/var2/2
    res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf(x+cutoff))/stats.norm.pdf(stats.norm.ppf(x))
    return(res1+res2)

def myfunc_marginal_2der(x, args):
    I=args[0];S=args[1];beta=args[2];var=args[3]
    inclusion_length=args[4]
    skipping_length=args[5]
    #print('test x');print(x);print(var)
    temp1=1/pow(x,2)/pow(1-x,2)*(((2*x-1)*(logit(beta)-logit(x))-1)/var-1);#temp1=(2*x-1)/pow(x,2)/pow(1-x,2)*((logit(beta)-logit(x)-1/(2*x-1))/var+1)
    temp2=I*skipping_length*((2*inclusion_length+skipping_length)*x+skipping_length*(1-x))/pow(x,2)/pow(inclusion_length*x+skipping_length*(1-x),2)
    temp3=S*inclusion_length*((inclusion_length+2*skipping_length)*(1-x)+inclusion_length*x)/pow(1-x,2)/pow(inclusion_length*x+skipping_length*(1-x),2)
    #print('test');print(beta);print(x);print(var);print(temp1);print(temp2);print(temp3)
    return(temp1-temp2-temp3)
    
def myfunc_marginal(x, *args):
    beta=x
    I=args[0];S=args[1];psi=args[2];var=args[3]
    inclusion_length=args[4];skipping_length=args[5]
    sum=0
    for i in range(len(psi)):
        new_psi=inclusion_length*psi[i]/(inclusion_length*psi[i]+skipping_length*(1-psi[i]))
        f1=I[i]*log(new_psi)+S[i]*log(1-new_psi)-pow(logit(psi[i])-logit(beta),2)/2/var-log(psi[i])-log(1-psi[i])
        f1_2der=abs(myfunc_marginal_2der(psi[i],[I[i],S[i],beta,var,inclusion_length,skipping_length]))
        sum+=(0.5*log(abs(f1_2der)+0.00001)-f1)
    return(sum)

def myfunc_marginal_der(x, *args):
    beta=x
    I=args[0];S=args[1];psi=args[2];var=args[3]
    inclusion_length=args[4];skipping_length=args[5]
    sum=0
    for i in range(len(psi)):
        new_psi=inclusion_length*psi[i]/(inclusion_length*psi[i]+skipping_length*(1-psi[i]))
        f1_3der=1*(2*psi[i]-1)/pow(psi[i],2)/pow(1-psi[i],2)/beta/(1-beta)/var
        f1_2der=myfunc_marginal_2der(psi[i],[I[i],S[i],beta,var,inclusion_length,skipping_length])
        #print('test2');print(f1_2der)
        f1_der=(logit(psi[i])-logit(beta))/beta/(1-beta)/var
        f1=I[i]*log(new_psi)+S[i]*log(1-new_psi)-pow(logit(psi[i])-logit(beta),2)/2/var+log(psi[i])+log(1-psi[i])
        sum+=(0.5*f1_3der/(f1_2der)-f1_der)
    return(sum)

def myfunc_marginal_1(x, *args):
    beta2=x;beta1=x+cutoff
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3]
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7]
    inclusion_length=args[8];skipping_length=args[9]
    return(myfunc_marginal(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length))

def myfunc_marginal_1_der(x, *args):
    beta2=x;beta1=x+cutoff
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3]
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7]
    inclusion_length=args[8];skipping_length=args[9]
    return(myfunc_marginal_der(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal_der(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length))

def myfunc_marginal_2(x, *args):
    beta2=x+cutoff;beta1=x
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3]
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7]
    inclusion_length=args[8];skipping_length=args[9]
    return(myfunc_marginal(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length))

def myfunc_marginal_2_der(x, *args):
    beta2=x+cutoff;beta1=x
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3]
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7]
    inclusion_length=args[8];skipping_length=args[9]
    return(myfunc_marginal_der(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal_der(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length))
    
def myfunc_individual(x,*args):
    I=args[0];S=args[1];beta=args[2];var=args[3]
    inclusion_length=args[4]
    skipping_length=args[5]
    new_psi=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x))
    return(-1*(I*log(new_psi)+S*log(1-new_psi)-(logit(x)-logit(beta))*(logit(x)-logit(beta))/2/var-log(x)-log(1-x)))

def myfunc_individual_der(x,*args):
    I=args[0];S=args[1];beta=args[2];var=args[3]
    inclusion_length=args[4]
    skipping_length=args[5]
    new_psi=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x))
    new_psi_der=inclusion_length*skipping_length/pow(inclusion_length*x+skipping_length*(1-x),2)
    return(-1*(I/new_psi*new_psi_der-S/(1-new_psi)*new_psi_der-(logit(x)-logit(beta))/var/x/(1-x)-1/x+1/(1-x) ))

def myfunc_likelihood(x, args):
    I=args[0];S=args[1];beta=args[2];var=args[3];sum=0;N=I+S
    #return(-1*(-log(sqrt((I+S)*x*(1-x)))-(I-(I+S)*x)*(I-(I+S)*x)/2/((I+S)*x*(1-x))-log(sqrt(var))-(x-beta)*(x-beta)/2/var))
    #print('debug');print(N);print(var);print(x);print(beta)
    if N==0:
        return(0)
    return(-0.5*((I-N*x)*(I-N*x)/(N*x)+(S-N*(1-x))*(S-N*(1-x))/(N*(1-x)))-log(sqrt(var))-(x-beta)*(x-beta)/2/var)

def MLE_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length)
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0
    beta_0=sum(psi1)/len(psi1)
    beta_1=sum(psi2)/len(psi2)
    var1=10*scipy.var(np.array(psi1)-beta_0)
    var2=10*scipy.var(np.array(psi2)-beta_1)
    if var1<=0.01:
        var1=0.01
    if var2<=0.01:
        var2=0.01
    #print('var1');print(var1);print('var2');print(var2)
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1
        #iteration of beta
        beta_0=sum(psi1)/len(psi1)
        beta_1=sum(psi2)/len(psi2)
        #print('var1');print(var1);print('var2');print(var2)
        #if abs(sum(psi1)/len(psi1)-sum(psi2)/len(psi2))>cutoff:
        if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff
        else:#minize psi1 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta1 = max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff
        #print('constrain_1xopt');print('theta');print(theta1);print(theta2);print(xopt)
        #else:
        #    theta1=sum(psi1)/len(psi1);theta2=sum(psi2)/len(psi2)
        beta_0=theta1;beta_1=theta2
        #iteration of psi
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0
        #print('constrain_2xopt')
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
            #likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length])
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
            #likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length])
        #print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0])
        psi1=new_psi1;psi2=new_psi2
        #print('count');print(count);print('previous_sum');print(previous_sum);print('current_sum');print(current_sum)
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum)
        previous_sum=current_sum
    #print('constrain');print(theta1);print(theta2);print(psi1);print(psi2);print(current_sum);print(likelihood_sum)
    #print(xopt)
    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]])
    #return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]])

def MLE_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length)
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0
    beta_0=sum(psi1)/len(psi1)
    beta_1=sum(psi2)/len(psi2)
    var1=10*scipy.var(np.array(psi1)-beta_0)
    var2=10*scipy.var(np.array(psi2)-beta_1)
    if var1<=0.01:
        var1=0.01
    if var2<=0.01:
        var2=0.01
    #print('var1');print(var1);print('var2');print(var2)
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1
        #iteration of beta
        beta_0=sum(psi1)/len(psi1)
        beta_1=sum(psi2)/len(psi2)
        xopt=fmin_l_bfgs_b(myfunc_multivar,[beta_0,beta_1],myfunc_multivar_der,args=[psi1,psi2,var1,var2],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1)
        beta_0=float(xopt[0][0])
        beta_1=float(xopt[0][1])
        #print('unconstrain_1xopt');print(xopt)
        #print('theta');print(beta_0);print(beta_1);print('theta_end')
        #iteration of psi
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
            #likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length])
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
            #likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length])
        #print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0])
        psi1=new_psi1;psi2=new_psi2
        #print('count');print(count);('previous_sum');print(previous_sum);print('current_sum');print(current_sum)
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum)
        previous_sum=current_sum
    #print('unconstrain');print(beta_0);print(beta_0+beta_1);print(psi1);print(psi2);print(current_sum);print(likelihood_sum)
    #print(xopt)
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2]])
    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]])
    #return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]])

def MLE_marginal_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    #initial value
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length)
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0
    beta_0=sum(psi1)/len(psi1)
    beta_1=sum(psi2)/len(psi2)
    var1=10*scipy.var(np.array(psi1)-beta_0)
    var2=10*scipy.var(np.array(psi2)-beta_1)
    if var1<=0.01:
        var1=0.01
    if var2<=0.01:
        var2=0.01
    #print('var1');print(var1);print('var2');print(var2)
    
    #MLE of the full likelihood
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1
        #iteration of beta
        beta_0=sum(psi1)/len(psi1)
        beta_1=sum(psi2)/len(psi2)
        xopt=fmin_l_bfgs_b(myfunc_multivar,[beta_0,beta_1],myfunc_multivar_der,args=[psi1,psi2,var1,var2],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1)
        beta_0=float(xopt[0][0])
        beta_1=float(xopt[0][1])
        #print('unconstrain_MLE_xopt');print(xopt)
        #print('theta');print(beta_0);print(beta_1);print('theta_end')
        #iteration of psi
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1])
            #likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length])
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1])
            #likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length])
        #print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0])
        psi1=new_psi1;psi2=new_psi2
        #print('count');print(count);print('previous_sum');print(previous_sum);print('current_sum');print(current_sum)
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum)
        previous_sum=current_sum
    #print('unconstrain');print(beta_0);print(beta_0+beta_1);print(psi1);print(psi2);print(current_sum);print(likelihood_sum)
    #print('unconstrain_psi_MLE');print(psi1);print(psi2)
    #print('unconstrain_beta_MLE');print(beta_0);print(beta_1)
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2]])
    
    #MLE of the marginal likelihood
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1
        #print('unconstrain_MLE_marginal_value_1_der');print(myfunc_marginal_der(beta_0,i1,s1,psi1,var1,effective_inclusion_length,effective_skipping_length))
        #print('unconstrain_MLE_marginal_value_2_der');print(myfunc_marginal_der(beta_1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length))
        xopt1=fmin_l_bfgs_b(myfunc_marginal,[beta_0],myfunc_marginal_der,args=[i1,s1,psi1,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
        xopt2=fmin_l_bfgs_b(myfunc_marginal,[beta_1],myfunc_marginal_der,args=[i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
        beta_0=float(xopt1[0][0])
        beta_1=float(xopt2[0][0])
        #print('unconstrain_1xopt');print(xopt1)
        #print('unconstrain_2xopt');print(xopt2)
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1])
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1])
        psi1=new_psi1;psi2=new_psi2
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum)
        previous_sum=current_sum
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2]])
    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]])

def MLE_marginal_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    #initial value
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length)
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0
    beta_0=sum(psi1)/len(psi1)
    beta_1=sum(psi2)/len(psi2)
    var1=10*scipy.var(np.array(psi1)-beta_0)
    var2=10*scipy.var(np.array(psi2)-beta_1)
    if var1<=0.01:
        var1=0.01
    if var2<=0.01:
        var2=0.01
    #print('var1');print(var1);print('var2');print(var2)
    
    #MLE of the full likelihood
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1
        #iteration of beta
        beta_0=sum(psi1)/len(psi1)
        beta_1=sum(psi2)/len(psi2)
        #print('var1');print(var1);print('var2');print(var2)
        #if abs(sum(psi1)/len(psi1)-sum(psi2)/len(psi2))>cutoff:
        if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff
        else:#minize psi1 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta1 = max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff
        #print('constrain_1xopt');print('theta');print(theta1);print(theta2);print(xopt)
        #else:
        #    theta1=sum(psi1)/len(psi1);theta2=sum(psi2)/len(psi2)
        beta_0=theta1;beta_1=theta2
        #iteration of psi
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0
        #print('constrain_2xopt')
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
            #likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length])
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
            #likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length])
        #print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0])
        psi1=new_psi1;psi2=new_psi2
        #print('count');print(count);print('previous_sum');print(previous_sum);print('current_sum');print(current_sum)
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum)
        previous_sum=current_sum
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2]])
    #print('constrain');print(theta1);print(theta2);print(psi1);print(psi2);print(current_sum);print(likelihood_sum)
    #print(xopt)
    
    #MLE of the marginal likelihood
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1
        if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_marginal_1,[beta_1],myfunc_marginal_1_der,args=[i1,s1,psi1,var1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            beta_1 = max(min(float(xopt[0]),1-cutoff),0);beta_0=beta_1+cutoff
        else:#minize psi1 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_marginal_2,[beta_0],myfunc_marginal_2_der,args=[i1,s1,psi1,var1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            beta_0 = max(min(float(xopt[0]),1-cutoff),0);beta_1=beta_0+cutoff
        #print('constrain_xopt');print(xopt)
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1)
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt)
        psi1=new_psi1;psi2=new_psi2
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum)
        previous_sum=current_sum
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2]])

    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]])
    #return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]])

    
#Random Sampling Function
def likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag,id):
    #print('testing'+str(id))
    if flag==0:
        return([1,1])
    else:
        res=MLE_marginal_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length)
        if abs(res[1][2]-res[1][3])<=cutoff:
            #print('1<=cutoff');print(res);print((res[1][2]-res[1][3]))
            return([1,res[0]])
        else:
            res_constrain=MLE_marginal_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length)
            #print('2>cutoff');print(res);print(res_constrain)
            #print(abs(res_constrain[0]-res[0]));print('2end')
            #print(abs(res_constrain[2]-res[2]));print('2end_marginal');            
            if len(i1)<=3:
                return([1-scipy.stats.chi2.cdf(2*(abs(res_constrain[0]-res[0])),1)])
            else:
                return([1-scipy.stats.chi2.cdf(2*(abs(res_constrain[0]-res[0])),1)])

#MultiProcessorFunction
def MultiProcessorPool(n_original_diff):
    i1=n_original_diff[0];i2=n_original_diff[1];s1=n_original_diff[2];s2=n_original_diff[3]
    effective_inclusion_length=n_original_diff[4];effective_skipping_length=n_original_diff[5]
    flag=n_original_diff[6];id=n_original_diff[7]
    P=likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag,id)
    return(P)

#Function for vector handling
def vec2float(vec):
    res=[]
    for i in vec:
        res.append(float(i))
    return(res)

#add 1 in both inclusion and skipping counts for robustness in small sample size
def vecAddOne(vec):
    res=[]
    for i in vec:
        res.append(i+1)
    return(res)

def vecprod(vec):
    res=1
    for i in vec:
        res=res*i
    return(res)

def vecadd(vec1,vec2):
    res=[]
    for i in range(len(vec1)):
        res.append(vec1[i]+vec2[i])
    return(res)
    
def vec2remove0psi(inc,skp):
    res1=[];res2=[]
    for i in range(len(inc)):
        if (inc[i]!=0) | (skp[i]!=0):
            res1.append(inc[i]);res2.append(skp[i])
    return([res1,res2])

def vec2psi(inc,skp,effective_inclusion_length,effective_skipping_length):
    psi=[]
    inclusion_length=effective_inclusion_length
    skipping_length=effective_skipping_length
    for i in range(len(inc)):
        if (float(inc[i])+float(skp[i]))==0:
            psi.append(0.5)
        else:
            psi.append(float(inc[i])/inclusion_length/(float(inc[i])/inclusion_length+float(skp[i])/skipping_length))
    return(psi)

def vec210(vec):
    res=[]
    for i in vec:
        if i>0:
            res.append(1)
        else:
            res.append(-1)
    return(res)

def myttest(vec1,vec2):
    if (len(vec1)==1) & (len(vec2)==1):
        res=stats.ttest_ind([vec1[0],vec1[0]],[vec2[0],vec2[0]])
    else:
        res=stats.ttest_ind(vec1,vec2)
    return(res)
def myorder(p,reverse):
    p_dict={}
    #print('len(p)');print(len(p))
    for index in range(len(p)):
        if not(p[index] in p_dict):
            p_dict[p[index]]=[]
        p_dict[p[index]].append(index+1)
    #print('p_dict');print(p_dict)
    o_index=sorted(p_dict,reverse=reverse);o=[]
    #print('o_index');print(o_index)
    for index in o_index:
        for this_p in p_dict[index]:
            o.append(this_p)
    return(o)

def mycummin(p):
    res=[]
    for i in range(len(p)):
        res.append(min(p[:(i+1)]))
    return(res)

def myFDR(p):
    lp=len(p)
    i=range(lp,0,-1)
    o=myorder(p,True)
    ro=myorder(o,False)
    p_new=[]
    for index in range(len(o)):
        p_new.append(p[o[index]-1]*lp/i[index])
    p_mycummin=mycummin(p_new)
    p_mycummin_new=[]
    for index in range(len(p_mycummin)):
        p_mycummin_new.append(min(p_mycummin[index],1))
    res=[]
    for index in range(len(p_mycummin_new)):
        res.append(p_mycummin_new[ro[index]-1])
    return(res)
    
def adjustPsi(inc,skp,effective_inclusion_length,effective_skipping_length):
    if effective_inclusion_length == 0:
        print ([inc,skp,effective_inclusion_length,effective_skipping_length])
    inc_adj=np.array(inc)/effective_inclusion_length
    skp_adj=np.array(skp)/effective_skipping_length
    return(','.join([str(i) for i in (inc_adj/(inc_adj+skp_adj)).tolist()]))

def myand(l1,l2):
    return([i for i in range(len(l1)) if l1[i] & l2[i] ==True])
    
######## end of functions ##############


################## actual process ##############
####
#### 1. STAR mapping
####
logging.debug("start mapping..")
if fastq==1: ## no junction files, start mapping
    try:
        doSTARMapping()
    except:
        logging.debug("There is an exception in mapping")
        logging.debug("Exception: %s" % sys.exc_info()[0])
        logging.debug("Detail: %s" % sys.exc_info()[1])
        for stack in traceback.extract_tb(sys.exc_info()[2]):
            logging.debug("Stack: %s" % stack)
        sys.exit(-1)
else: ## junction files are provided
    logging.debug("junction files are provided. skip mapping ..")

logging.debug("done mapping..")


####
#### 2. CIRCexplorer2 detection
####
logging.debug("Detect circRNA by CIRCexplorer2 parse")
try:
    doCIRCexplorer2()
except:
    logging.debug("There is an exception in detection")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    for stack in traceback.extract_tb(sys.exc_info()[2]):
            logging.debug("Stack: %s" % stack)
    sys.exit(-1)
logging.debug("done circRNA detection..")

####
#### 3. annotation circRNA
####
logging.debug('Annotate circRNA by CIRCexplorer2 annotate')
try:
    annotateCirc()
except:
    logging.debug("There is an exception in annoation")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        logging.debug("Stack: %s" % stack)
    sys.exit(-1)
logging.debug("done circRNA annotation..")

####
#### 4. CIRCexplorer2 detection for each sample
####
logging.debug('Detect circRNA by CIRCexplorer2 parse for each sample')
try:
    doCIRCexplorer2foreach()
except:
    logging.debug("There is an exception in detection")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        logging.debug("Stack: %s" % stack)
    sys.exit(-1)
logging.debug("done circRNA detection..")

####
#### 5. Combine back-splicing
####
#
logging.debug('Combine back-splicing')
try:
    matBS,annoBS=getCircRNAdetail()
except:
    logging.debug("There is an exception in combine back-splicing junction")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        logging.debug("Stack: %s" % stack)
    sys.exit(-1)
logging.debug("done combine back-splicing junction..")
#
matBS['key']=annoBS.iloc[:,0]+'|'+annoBS.iloc[:,1].map(str)+'|'+annoBS.iloc[:,2].map(str)
matBS['startKey']=annoBS.iloc[:,0]+'|'+annoBS.iloc[:,1].map(str)
matBS['endKey']=annoBS.iloc[:,0]+'|'+annoBS.iloc[:,2].map(str)

####
#### 6. Combine linear splicing
####
#
logging.debug('Combine linear splicing')
try:
    matSJ,annoSJ,list1_L,list2_L,list1_R,list2_R=combineLinearJunction()
except:
    logging.debug("There is an exception in combine linear splicing junction")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        logging.debug("Stack: %s" % stack)
    sys.exit(-1)
logging.debug("done combine linear splicing junction..")
'''
#### 7. Collapse back-splicing
####
#

logging.debug('Collapse back-splicing junction')
try:
    matCBJ,annoCBJ=collapseBJ()
except:
    logging.debug("There is an exception in collapse back-splicing junction")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        logging.debug("Stack: %s" % stack)
    sys.exit(-1)
logging.debug("done collapse back-splicing junction..")
'''
####
#### 7. Get linear transcript detail
#
logging.debug('Get linear transcript detail')
try:
    annoTrans=getannoTrans()
    annoTrans.rename(index=dict(zip([ i  for i in range(annoBS.iloc[:,0].size)],matBS.index)),inplace=True)
except:
    logging.debug("There is an exception in get linear transcript detail")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    for stack in traceback.extract_tb(sys.exc_info()[2]):
        logging.debug("Stack: %s" % stack)
    sys.exit(-1)
logging.debug("done get linear transcript detail..")

#### 8. Cal effective length
####
#
logging.debug('Calculate effective length')
matBS=matBS.iloc[:,0:sampleNum]
matSBJ=matSJ

matStat=pd.DataFrame('',index=range(len(annoBS.index)),columns=range(7))
matStat.rename(index=dict(zip([ i  for i in range(annoBS.iloc[:,0].size)],matBS.index)),columns={0:'id',1:'inc1',2:'inc2',3:'bs1',4:'bs2',5:'effective_inclusion_length',6:'effective_bs_length'},inplace=True)
matStat['id']=matBS.index
matStat['effective_inclusion_length']=4*(readLength-Anchorlength)
matStat['effective_bs_length']=2*(readLength-Anchorlength)
annoBS.iloc[:,18]=annoBS.iloc[:,18].map(int)

    
for i in range(matStat.iloc[:,0].size):
    #print(i)
    if matStat.iat[i,6] > annoBS.iat[i,18]:
        matStat.iat[i,6]=annoBS.iat[i,18]
        matStat.iat[i,5]=matStat.iat[i,5]-2*(readLength-Anchorlength)+annoBS.iat[i,18]

        
for i in matStat.index[myand((annoTrans.iloc[:,0]==0).tolist() , (annoTrans.iloc[:,2]==0).tolist())]:
    matStat.at[i,'effective_inclusion_length']=0


for i in matStat.index[myand((annoTrans.iloc[:,0]==0).tolist(),(annoTrans.iloc[:,2]!=0).tolist())]:
    if (readLength-Anchorlength)>annoBS.at[i,'exonLen']:
        matStat.at[i,'effective_inclusion_length']=annoBS.at[i,'exonLen']+readLength-Anchorlength
    else:
        matStat.at[i,'effective_inclusion_length']=2*(readLength-Anchorlength)


for i in matStat.index[myand((annoTrans.iloc[:,2]==0).tolist() ,(annoTrans.iloc[:,0]!=0).tolist())]:
    if (readLength-Anchorlength)>annoBS.at[i,'exonLen']:
        matStat.at[i,'effective_inclusion_length']=annoBS.at[i,'exonLen']+readLength-Anchorlength
    else:
        matStat.at[i,'effective_inclusion_length']=2*(readLength-Anchorlength)

logging.debug('done calculate effective length')

#### 9. Statistic analysis
####
#
logging.debug('Statistic analysis')

matStat.loc[:,'inc1']=matSBJ.iloc[:,0].map(str).tolist()
matStat.loc[:,'inc2']=matSBJ.iloc[:,len(sample_1)].map(str).tolist()
matStat.loc[:,'bs1']=matBS.iloc[:,0].map(str).tolist()
matStat.loc[:,'bs2']=matBS.iloc[:,len(sample_1)].map(str).tolist()

for i in range(1,len(sample_1)):
    tmp=matStat.loc[:,'inc1']+','+matSBJ.iloc[:,i].map(str)
    matStat.loc[:,'inc1']=tmp.tolist()
for i in range(len(sample_1)+1,sampleNum):
    tmp=matStat.loc[:,'inc2']+','+matSBJ.iloc[:,i].map(str)
    matStat.loc[:,'inc2']=tmp.tolist()

for i in range(1,len(sample_1)):
    tmp=matStat.loc[:,'bs1']+','+matBS.iloc[:,i].map(str)
    matStat.loc[:,'bs1']=tmp.tolist()
for i in range(len(sample_1)+1,sampleNum):
    tmp=matStat.loc[:,'bs2']+','+matBS.iloc[:,i].map(str)
    matStat.loc[:,'bs2']=tmp.tolist()

    
    

  
list_n_original_diff=[];probability=[];psi_list_1=[];psi_list_2=[]
for i in range(matStat.iloc[:,0].size):
    inc1=matSBJ.iloc[i,0:len(sample_1)].tolist()
    inc2=matSBJ.iloc[i,len(sample_1):sampleNum].tolist()
    skp1=matBS.iloc[i,0:len(sample_1)].tolist()
    skp2=matBS.iloc[i,len(sample_1):sampleNum].tolist()
    effective_inclusion_length=int(matStat.iat[i,5])
    effective_skipping_length=int(matStat.iat[i,6])
    inc1=vec2float(inc1);skp1=vec2float(skp1);inc2=vec2float(inc2);skp2=vec2float(skp2)

    if (vecprod(vecadd(inc1,skp1))==0) | (vecprod(vecadd(inc2,skp2))==0) | (effective_inclusion_length==0):
        psi_list_1.append('')
        psi_list_2.append('')
        if effective_inclusion_length ==0:
            effective_inclusion_length=readLength
        list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,i])
    else:
        psi_list_1.append(adjustPsi(inc1,skp1,effective_inclusion_length,effective_skipping_length))
        psi_list_2.append(adjustPsi(inc2,skp2,effective_inclusion_length,effective_skipping_length))
        inc1=vecAddOne(inc1);skp1=vecAddOne(skp1);inc2=vecAddOne(inc2);skp2=vecAddOne(skp2)
        list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,1,i])
rho=0.9

if threads>1:
    pool=Pool(processes=threads)
    probability=pool.map(MultiProcessorPool,list_n_original_diff)
else:
    for i in range(len(list_n_original_diff)):
        probability.append(MultiProcessorPool(list_n_original_diff[i]))
        
P=pd.DataFrame(probability).iloc[:,0].tolist()

logging.debug("done Statistic analysis..")
logging.debug("Prepare DEBKS_results output..")
matStat['chr']=annoBS.loc[:,'chr'].tolist()
matStat['start']=annoBS.loc[:,'start'].tolist()
matStat['end']=annoBS.loc[:,'end'].tolist()
matStat['strand']=annoBS.loc[:,'strand'].tolist()
matStat['exonCount']=annoBS.loc[:,'exonCount'].tolist()
matStat['exonSizes']=annoBS.loc[:,'exonSizes'].tolist()
matStat['exonOffsets']=annoBS.loc[:,'exonOffsets'].tolist()
matStat['geneID']=annoBS.loc[:,'geneName'].tolist()
matStat['isoformID']=annoBS.loc[:,'isoformName'].tolist()
matStat['flankIntron']=annoBS.loc[:,'flankIntron'].tolist()
matStat['linearExonL']=annoTrans.iloc[:,0].tolist()
matStat['linearExonR']=annoTrans.iloc[:,3].tolist()
matStat['PBSI1']=psi_list_1
matStat['PBSI2']=psi_list_2
matStat['SJL1']=list1_L
matStat['SJL2']=list2_L
matStat['SJR1']=list1_R
matStat['SJR2']=list2_R
matStat['P']=P
matStat=matStat[matStat.loc[:,'PBSI1'].map(len)>0]
FDR=myFDR(matStat['P'].tolist())
matStat['FDR']=FDR
title=['chr','start','end','strand','exonCount','exonSizes',
'exonOffsets','geneID','isoformID','flankIntron',
'linearExonL','linearExonR','SJL1','SJL2','SJR1','SJR2','inc1','inc2','bs1','bs2',
'effective_inclusion_length','effective_bs_length','PBSI1','PBSI2','P','FDR']
matStat=matStat.ix[:,title]
matStat.to_csv(outPath+'/DEBKS_output/DEBKS_results.txt',index=False,sep='\t',header=True)
logging.debug("done DEBKS_results output..")
### clean up temp folder if necessary ###
if keepTemp==0: ## delete temp folder, by default
  os.system('rm -rf '+ tempPath)
  logging.debug("Temp folder is deleted..");  
#
#############
## calculate total running time
#############
logging.debug("Program ended")
currentTime = time.time()
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))
sys.exit(0)
