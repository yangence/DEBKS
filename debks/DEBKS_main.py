#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
author: Zelin Liu
email: zlliu@bjmu.edu.cn
license: GPL3
detail: Analysis differential back-spliced in of rRNA-depleted RNA-seq
"""
from .version import __version__
# DEBKS version
DEBKS_ver = __version__
# help info
helpInfo = '''
Usage (with FASTQ files):\n\tDEBKS -g genomeFasta -s1 s1File -s2 s2File  -STARindex STARindexDir -gtf gtfFile -o outDir -read readType -len readLength [options]*

Usage (with junction files):\n\tDEBKS -g genomeFasta -s1CJ s1CJFile -s2CJ s2CJFile -s1SJ s1SJFile -s2SJ s2SJFile -gtf gtfFile -o outDir -read readType -len readLength [options]*

DEBKS '''+DEBKS_ver+''' : Analysis of differential percent back-spliced in (C) Zelin Liu, 2020
###Required Parameters:

	-g          <str>       Genome Fasta file

	-gtf        <str>       GTF file

	-o          <str>       Output directory of the result files

	-read       <str>       RNA-seq reads are single- or pair-end reads. [single, pair]

	-len        <int>       Read length of RNA-seq reads

###Only required with FASTQ files

	-STARindex  <str>       STAR alignment index directory

	-s1         <str>       Fastq files of sample 1, replicates in different lines, paired files are separated by semicolon

	-s2         <str>       Fastq files of sample 2, replicates in different lines, paired files are separated by semicolon

###Only required with junction files

	-s1CJ       <str>       Chimeric junction of sample 1 group, replicates in different lines

	-s2CJ       <str>       Chimeric junction of sample 2 group, replicates in different lines

	-s1SJ       <str>       Spliced junction of sample 1 group, replicates in different lines

	-s2SJ       <str>       Spliced junction of sample 2 group, replicates in different lines

###Optional Parameters:

	-h, --help              Show this help message and exit

	-p                      Sample 1 group and sample 2 group is paired

	-n          <int>       Required total juction reads in all samples to filter out low expressed circRNAs [2*samples]

	-t          <int>       Number of processors [1]

	-c          <float>     Required PBSI difference cutoff between the two samples [0.05]

	-a          <int>       Minimum overhang length for counting chimeric or splicing junctions [6]

	-keepTemp               Keep the temporary files. Disable by default.
'''
import sys
# python version
'''
if sys.version_info[0] < 3:
    print ("Python version error: python 3 or above is required")
    sys.exit(-1)
    
'''
# import libraries
import subprocess
import re,os,warnings,math,scipy,itertools,logging,time,datetime,pysam,traceback,gzip
from multiprocessing import Pool
from math import log
import pandas as pd
import numpy as np
from numpy import *
np.seterr(divide='ignore',invalid='ignore')
########## functions ############
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

def listToString(x):
    rVal = ''
    for a in x:
        rVal += a+' '
    return rVal
####################

# software required
# STAR, gtfToGenePred, CIRCexplorer2
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
cutoff=0.05
nThreshold=-1
groupPair=False
keepTemp = 0 ## keep temp
### checking out the argument names
validArgList=['-g','-s1','-s2','-STARindex','-s1CJ','-s2CJ','-s1SJ','-s2SJ','-read','-t','-gtf','-o','-len','-c','-keepTemp','-a','-h','--help','-p','-n']
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
    elif (sys.argv[paramIndex] == '-n'):  ## anchor length for chimeric junction and linear junction
        paramIndex += 1  ## increase index
        nThreshold = int(sys.argv[paramIndex])
    elif (sys.argv[paramIndex] == '-keepTemp'):  ## keep temp files, no value needed
        keepTemp = 1  ## keep temp files
    elif (sys.argv[paramIndex] == '-p'):  ## group is paired
        groupPair = True
    elif (sys.argv[paramIndex] == '-h' or sys.argv[paramIndex] == '--help'):  ## print helpInfo
        print(helpInfo)
        sys.exit()

### checking out the required arguments

if (genomeFasta == '' or ((s1=='' or  s2=='' or genomeDir == '') and (s1CJ=='' or  s2CJ=='' or s1SJ=='' or  s2SJ=='')) or gtf==''  or outDir=='' or readLength==0 or readType==''): ### at least one required param is missing
    print('Not enough arguments!!')
    print(helpInfo)
    sys.exit()

if groupPair:
    from .rMATs_paired import *
else:
    from .rMATs import *

CAlength=Anchorlength
LAlength=Anchorlength
genomeFasta=os.path.abspath(genomeFasta)
gtf=os.path.abspath(gtf)

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
if gtf[(len(gtf)-3):]=='.gz':
    tempGTF_file = gzip.open(gtf,'rb') ## open gtf file
else:
    tempGTF_file = open(gtf,'r')
    
for line in tempGTF_file:
  if gtf[(len(gtf)-3):]=='.gz':
    line=line.decode()
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
            print("Error: Incorrect CJ file format. chimeric junction file must have 14 tab-delimited columns.")
            sys.exit()
    for CJfile in sample_2_CJ[0:2]: ## examine each file
        if checkLineNum(CJfile,14) == 0:
            print("Error: Incorrect CJ file format. chimeric junction file must have 14 tab-delimited columns.")
            sys.exit()
    for SJfile in sample_1_SJ[0:2]: ## examine each file
        if checkLineNum(SJfile,9) == 0:
            print("Error: Incorrect SJ file format. chimeric junction file must have 9 tab-delimited columns.")
            sys.exit()
    for SJfile in sample_2_SJ[0:2]: ## examine each file
        if checkLineNum(SJfile,9) == 0:
            print("Error: Incorrect SJ file format. chimeric junction file must have 9 tab-delimited columns.")
            sys.exit()
            
sampleNum=len(sample_1)+len(sample_2)
if len(sample_1)==1 or len(sample_2)==1:
    print('Error: Each group need at least two samples!')
    sys.exit()
if nThreshold==-1:
    nThreshold=2*sampleNum
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

def main():
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
            cmd = 'STAR --genomeLoad LoadAndKeep --chimSegmentMin ' + str(CAlength) + ' --chimOutType WithinBAM  --runThreadN ' + str(threads) + '  --outSAMtype BAM Unsorted'
            cmd += ' --alignSJDBoverhangMin ' + str(LAlength) +' --alignSJoverhangMin '+ str(LAlength) + '  --genomeDir '+genomeDir
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
            cmd = 'STAR --genomeLoad LoadAndKeep --chimSegmentMin ' + str(CAlength) + ' --chimOutType WithinBAM  --runThreadN ' + str(threads) + ' --outSAMtype BAM Unsorted'
            cmd += ' --alignSJDBoverhangMin ' + str(LAlength) +' --alignSJoverhangMin '+ str(LAlength) + ' --genomeDir '+genomeDir
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
        cmd='STAR --genomeLoad Remove --outFileNamePrefix '+tempPath+'/ --genomeDir '+genomeDir
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
        os.system('awk -F "\t" \'{if($4>='+str(nThreshold)+'){print $0}}\' '+ tempPath+'/circularRNA_known.txt > '+tempPath+'/circularRNA_filter.txt')
        annoBS=pd.read_csv(tempPath+'/circularRNA_filter.txt',header=None,sep='\t')
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
            BSJbed=pd.read_csv(rTempFolder+'/circularRNA_known.txt',header=None,sep='\t')
            BSrowIndex= dict(zip([ i  for i in range(BSJbed.iloc[:,0].size)],BSJbed[0]+'|'+BSJbed[1].map(str)+'|'+BSJbed[2].map(str)))
            BSJbed.rename(index=BSrowIndex,inplace=True)
            BSJbed=BSJbed.iloc[[not i for i in BSJbed.index.duplicated()],:]
            BSJbed=BSJbed.reindex(annoBS.index,axis=0)
            matBS.iloc[:,rr]=BSJbed.iloc[:,3]
        for rr in range(0,len(sample_2)): ## for each replicate of sample_2
            rTempFolder = s2rPath+str(rr+1)
            BSJbed=pd.read_csv(rTempFolder+'/circularRNA_known.txt',header=None,sep='\t')
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
            SJtab=pd.read_csv(rTempFolder+'/SJ.out.tab',header=None,sep='\t')
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
            SJtab=pd.read_csv(rTempFolder+'/SJ.out.tab',header=None,sep='\t')
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
        annoFile=pd.read_csv(tempPath+'/GenePred.CIRC.txt',header=None,sep='\t')
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
    ##### PBSI #############
    def psi2PBSI(psi):
        pbsi=[]
        for i in psi:
            eachPbsi=[]
            eachPsi=i.split(',')
            for j in eachPsi:
                if j!='':
                    eachPbsi.append(str(1-float(j)))
                else:
                    eachPbsi.append('')
            pbsi.append(','.join(eachPbsi))
        return(pbsi)
    
    ######## end of functions ##############

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
    matStat.loc[:,'bs1']=matBS.iloc[:,0].map(int).map(str).tolist()
    matStat.loc[:,'bs2']=matBS.iloc[:,len(sample_1)].map(int).map(str).tolist()

    for i in range(1,len(sample_1)):
        tmp=matStat.loc[:,'inc1']+','+matSBJ.iloc[:,i].map(str)
        matStat.loc[:,'inc1']=tmp.tolist()
    for i in range(len(sample_1)+1,sampleNum):
        tmp=matStat.loc[:,'inc2']+','+matSBJ.iloc[:,i].map(str)
        matStat.loc[:,'inc2']=tmp.tolist()

    for i in range(1,len(sample_1)):
        tmp=matStat.loc[:,'bs1']+','+matBS.iloc[:,i].map(int).map(str)
        matStat.loc[:,'bs1']=tmp.tolist()
    for i in range(len(sample_1)+1,sampleNum):
        tmp=matStat.loc[:,'bs2']+','+matBS.iloc[:,i].map(int).map(str)
        matStat.loc[:,'bs2']=tmp.tolist()

        
        
    def vec_remove_na_pair(inc1,skp1,inc2,skp2):
        res_inc1=[];res_skp1=[];res_inc2=[];res_skp2=[];
        for i in range(len(inc1)):
            if ((float(inc1[i])+float(skp1[i]))>0) & ((float(inc2[i])+float(skp2[i]))>0):
                res_inc1.append(inc1[i]);
                res_skp1.append(skp1[i]);
                res_inc2.append(inc2[i]);
                res_skp2.append(skp2[i]);
        return([res_inc1,res_skp1,res_inc2,res_skp2])
    def vec_remove_na(inc,skp):
        res_inc=[];res_skp=[];
        for i in range(len(inc)):
            if (float(inc[i])+float(skp[i]))>0:
                res_inc.append(inc[i]);
                res_skp.append(skp[i]);
        return([res_inc,res_skp]);
        
    list_n_original_diff=[];psi_list_1=[];psi_list_2=[]
    for i in range(matStat.iloc[:,0].size):
        inc1=matSBJ.iloc[i,0:len(sample_1)].tolist()
        inc2=matSBJ.iloc[i,len(sample_1):sampleNum].tolist()
        skp1=matBS.iloc[i,0:len(sample_1)].tolist()
        skp2=matBS.iloc[i,len(sample_1):sampleNum].tolist()
        effective_inclusion_length=int(matStat.iat[i,5])
        effective_skipping_length=int(matStat.iat[i,6])
        inc1=vec2float(inc1);skp1=vec2float(skp1);inc2=vec2float(inc2);skp2=vec2float(skp2)
        if effective_inclusion_length <=0:
            effective_inclusion_length=readLength
        if groupPair:
            temp=vec_remove_na_pair(inc1,skp1,inc2,skp2);inc1_nona=temp[0];skp1_nona=temp[1];inc2_nona=temp[2];skp2_nona=temp[3];
            if (len(inc1_nona)==0) | (len(inc2_nona)==0):
                psi_list_1.append('')
                psi_list_2.append('')
                list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,i]);
            else:
                psi_list_1.append(adjustPsi(inc1,skp1,effective_inclusion_length,effective_skipping_length))
                psi_list_2.append(adjustPsi(inc2,skp2,effective_inclusion_length,effective_skipping_length))
                inc1=inc1_nona;skp1=skp1_nona;inc2=inc2_nona;skp2=skp2_nona;
                inc1=vecAddOne(inc1);skp1=vecAddOne(skp1);inc2=vecAddOne(inc2);skp2=vecAddOne(skp2);
                #psi_avg_1=sum(array(inc1)/effective_inclusion_length/(array(inc1)/effective_inclusion_length+array(skp1)/effective_skipping_length))/len(inc1);
                #psi_avg_2=sum(array(inc2)/effective_inclusion_length/(array(inc2)/effective_inclusion_length+array(skp2)/effective_skipping_length))/len(inc2);
                if (len(inc1)==1) |(len(inc2)==1):
                    list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,i]);
                else:
                    list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,1,i]);
        else:
            temp1=vec_remove_na(inc1,skp1);temp2=vec_remove_na(inc2,skp2);
            inc1_nona=temp1[0];skp1_nona=temp1[1];inc2_nona=temp2[0];skp2_nona=temp2[1];
            if (len(inc1_nona)==0) | (len(inc2_nona)==0) | (effective_inclusion_length==0):
                psi_list_1.append('')
                psi_list_2.append('')
                list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,i])
            else:
                psi_list_1.append(adjustPsi(inc1,skp1,effective_inclusion_length,effective_skipping_length))
                psi_list_2.append(adjustPsi(inc2,skp2,effective_inclusion_length,effective_skipping_length))
                inc1=inc1_nona;skp1=skp1_nona;inc2=inc2_nona;skp2=skp2_nona;
                inc1=vecAddOne(inc1);skp1=vecAddOne(skp1);inc2=vecAddOne(inc2);skp2=vecAddOne(skp2)
                list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,1,i])
    P=getrMATS_P(list_n_original_diff,cutoff,threads)
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
    matStat['PBSI1']=psi2PBSI(psi_list_1)
    matStat['PBSI2']=psi2PBSI(psi_list_2)
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
    matStat=matStat.loc[:,title]
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

if __name__ == '__main__':
    main()
