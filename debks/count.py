'''
Usage: DEBKS count -c circ_pos (-n name | -f file) [-l library_type] [-t threads]  [-a hangover_len]  [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -c circ_pos                 File of circRNA position in tab format (the first three columns are chr, start, end).
    -n name                     Names of sorted bam file, separated by comma (e.g. bam1,bam2,bam3).
    -f file                     File includes names of sorted bam file, separated by line breaks.
    -l library_type             0 represent unstrandard, 1 represent fr-firststrand (e.g., dUTP protocol), 2 represent fr-secondstrand [default: 0].
    -t threads                  Number of threads to detect DE circRNA [default: 4].
    -a hangover_len             The minimum hangover length for junction quantification [default: 6].
    -o output                   Output prefix [default: linear].
'''
import pysam,pandas as pd, numpy as np,time,os,sys,traceback
from multiprocessing import Pool
def fileCheck(file):
    for i in file:
        if not os.path.exists(i):
            sys.exit('ERROR: %s is not exist!!!' % i)

def readCircFile(fileName):
    if not os.path.exists(fileName):
        sys.exit('ERROR: %s is not exist!!!' % fileName)
    fff=open(fileName)
    fl=fff.readline().strip().split('\t')
    fff.close()
    if len(fl) <3:
        sys.exit('ERROR: The number of columns in %s is not enough!!!' % fileName)
    isHeader=None
    for i in fl[1]:
        if i not in ['0','1','2','3','4','5','6','7','8','9']:
            isHeader=0
            break
    mat=pd.read_csv(fileName,sep='\t',header=isHeader)
    return(mat)

def count(options):
    threads=int(options['-t'])
    hangOver=int(options['-a'])
    if options['-n']:
        fileArr=options['-n'].split(',')
    else:
        if os.path.exists(options['-f']):
            fileArr=[i.strip() for i in open(options['-f']).readlines()]
        else:
            sys.exit('ERROR: %s is not exist!!!' % options['-f'])
    fileCheck(fileArr)
    outJun=options['-o']+'_jun.txt'
    outDetail=options['-o']+'_detail.txt'
    global circPos
    circPos=readCircFile(options['-c'])
    threads=min(threads,len(fileArr))
    stime = time.time()
    print('Begin counting %d circRNAs from %d files' % (circPos.shape[0],len(fileArr)))
    junDfList=[]
    taskArr=[]
    for i in fileArr:
        taskArr.append([i,options['-l'],hangOver])
    if threads>1:
        pool=Pool(processes=threads)
        junDfList=pool.map(coutSample,taskArr)
        pool.close()
        pool.join()
    else:
        for i in range(len(taskArr)):
            junDfList.append(coutSample(taskArr[i]))
    jun_num=[]
    jun_detail=[]
    for i in junDfList:
        jun_num.append(i['readN'].tolist())
        jun_detail.append(i['jun'].tolist())
    jun_num=pd.DataFrame(jun_num).T
    jun_detail=pd.DataFrame(jun_detail).T
    jun_num=pd.concat([circPos.iloc[:,0:3],jun_num],axis=1)
    jun_detail=pd.concat([circPos.iloc[:,0:3],jun_detail],axis=1)
    etime=time.time()-stime
    jun_num.to_csv(outJun,sep='\t',index=False,header=False)
    jun_detail.to_csv(outDetail,sep='\t',index=False,header=False)
    print ("Finished in %d seconds" % (etime))

#### function part ####
def getReadInfo(read):
    rS=read.reference_start
    rE=read.reference_end
    qS=read.query_alignment_start
    qE=read.query_alignment_end
    readcigar=read.cigar
    exonS=[rS]
    exonE=[]
    currentPos=rS-1
    queryS=[qS]
    queryE=[]
    queryPos=qS-1
    i=0
    for each in readcigar:
        i=i+1
        if i == len(readcigar):
            if each[0] ==0:
                currentPos+=each[1]
                queryPos+=each[1]
            elif each[0] == 1:
                queryPos+=each[1]
            elif each[0] == 2:
                currentPos+=each[1]
            exonE.append(currentPos)
            queryE.append(queryPos)
        else:
            if each[0] ==0:
                currentPos+=each[1]
                queryPos+=each[1]
            elif each[0] == 1:
                queryPos+=each[1]
            elif each[0] == 2:
                currentPos+=each[1]
            elif each[0] == 3:
                exonE.append(currentPos)
                queryE.append(queryPos)
                currentPos+=each[1]
                exonS.append(currentPos+1)
                queryS.append(queryPos+1)
    return([read.query_name,read.reference_name,rS,rE,qS,qE,list(zip(exonS,exonE)),list(zip(queryS,queryE))],read.flag)

def judgePos(read,pos,start=1,hangOver=6):
    isMap=0
    count=0
    candi=0
    tmpKey=''
    for i in read[6]:
        if i[0]<=pos and i[1]>=pos:
            isMap=1
            disDiff=pos-i[0]
            break
        count+=1
    if isMap:
        readPos=read[7][count][0]+disDiff
    else:
        return([candi,tmpKey])
    if readPos>=read[4]+hangOver and readPos<=read[5]-hangOver:
        if start:
            if read[6][0][0]+disDiff<=pos:
                candi=1
                if count>0:
                    if read[6][count][0]==pos:
                        tmpKey=str(read[1])+':'+str(read[6][count-1][1])+'-'+str(pos)
                    else:
                        tmpKey='noJun'
                else:
                    tmpKey='noJun'

        else:
            if read[6][len(read[6])-1][1]-disDiff>=pos:
                candi=1
                if count<len(read[6])-1:
                    if read[6][count][1]==pos:
                        tmpKey=str(read[1])+':'+str(pos)+'-'+str(read[6][count+1][0])
                    else:
                        tmpKey='noJun'
                else:
                    tmpKey='noJun'
    return([candi,tmpKey])
    
def countJun(chr,pos,bamfile,start=1,strandType='0',hangOver=6):
    pos=int(pos)
    try:
        iter=bamfile.fetch(chr,pos-1,pos)
    except:
        print("ERROR: %s" % sys.exc_info()[0])
        print("Detail: %s" % sys.exc_info()[1])
        for stack in traceback.extract_tb(sys.exc_info()[2]):
            print("Stack: %s" % stack)
        sys.exit('ERROR: make sure %s is suitable!!!' % bamfile.filename)
    countJunc={}
    readName=[]
    for x in iter:
        isCount=[0]
        eachRead,readFlag=getReadInfo(x)
        if strandType=='0':
            isCount=judgePos(eachRead,pos-1,start,hangOver)
        elif strandType=='1':
            if readFlag & 1==1:
                if readFlag & 80 ==80 or readFlag & 160 == 160:
                    isCount=judgePos(eachRead,pos-1,start,hangOver)
            else:
                if readFlag & 16==16:
                    isCount=judgePos(eachRead,pos-1,start,hangOver)
        elif strandType=='2':
            if readFlag & 1==1:
                if readFlag & 96 == 96 or readFlag & 144 == 144:
                    isCount=judgePos(eachRead,pos-1,start,hangOver)
            else:
                if readFlag & 32 == 32:
                    isCount=judgePos(eachRead,pos-1,start,hangOver)
        else:
             sys.exit('Unknown RNA-seq library type.')
        if isCount[0]:
            if countJunc.__contains__(isCount[1]):
                countJunc[isCount[1]]+=1
            else:
                countJunc[isCount[1]]=1
            readName.append(eachRead[0])
    countJuncList=[]
    for i in countJunc.keys():
        countJuncList.append(i+':'+str(countJunc[i]))
    return('|'.join(countJuncList),readName)

def countTwoJunc(x,bamfile,strandType='0',hangOver=6):  #1st: chr; 2nd: start; 3rd: end; # 1-based position
    junS,readS=countJun(x[0],x[1],bamfile,start=1,strandType=strandType,hangOver=hangOver)
    junE,readE=countJun(x[0],x[2],bamfile,start=0,strandType=strandType,hangOver=hangOver)
    readCount=len(set(readS+readE))
    return(junS+';'+junE,readCount)

def coutSample(x):
    fileName=x[0]
    strandType=x[1]
    hangOver=x[2]
    try:
        bamfile=pysam.AlignmentFile(fileName,'rb')
        
    except:
        print("ERROR: %s" % sys.exc_info()[0])
        print("Detail: %s" % sys.exc_info()[1])
        for stack in traceback.extract_tb(sys.exc_info()[2]):
            print("Stack: %s" % stack)
        sys.exit('ERROR: make sure %s is suitable!!!' % fileName)
    
    junDf=circPos.apply(lambda x: countTwoJunc(x,bamfile,strandType,hangOver),axis=1).tolist()
    junDf=pd.DataFrame(junDf,columns=['jun','readN'])
    return(junDf)