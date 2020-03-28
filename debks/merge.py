'''
Usage: DEBKS  merge (-s software | -d designate) (-n name | -f file)  [-b pos_based]  [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -s software                 Name of software used to detect circRNA (e.g. ciri2, circexplorer2, find_circ).
    -d designate                Designate the column postion of chr,start,end,circ,linear in the output file of detection software (e.g. 1,2,3,4,5).
                                The first four postion is required.
    -b pos_based                0- or 1-based position for designated start and end position, separated by comma [default: 1,1].
    -n name                     Names of output file of detection software, separated by comma (e.g. s1,s2,s3).
    -f file                     File includes names of output file of detection software, separated by line breaks.
    -o output                   Output prefix [default: merge].
'''
import pysam,pandas as pd, numpy as np,time,os,sys,traceback

def fileCheck(file):
    for i in file:
        if not os.path.exists(i):
            sys.exit('ERROR: %s is not exist!!!' % i)

def checkD(x):
    x=x.split(',')
    try:
        arr=[int(i) for i in x]
    except:
        sys.exit('ERROR: make sure "-d" options is suitable!!!')
    if len(arr) not in [4,5]:
        sys.exit('ERROR: "-d" options required four or five positions!!!')
    if 0 in arr:
        sys.exit('ERROR: "-d" options is 1-based, 0 position is not allowed!!!')
    return(arr)

def adjPos(pos_mat,based):
    based=based.split(',')
    if len(based)!=2:
        sys.exit('ERROR: "-b" only accept two parameter separated by comma!!!')
    if based[0]=='0':
        pos_mat['start']=pos_mat['start'].map(lambda x:int(x)+1).tolist()
    if based[1]=='0':
        pos_mat['start']=pos_mat['start'].map(lambda x:int(x)+1).tolist()
    return(pos_mat)
def merge(options):
    output=options['-o']
    soft_all=['ciri2','circexplorer2','find_circ']
    pos_ciri2=[2,3,4,5,7] # 1-based
    pos_circexplorer2=[1,2,3,4] # start: 0-based  end: 1-based
    pos_find_circ=[1,2,3,7] # Both start and end are 0-based
    pos_dict={'ciri2':pos_ciri2,'circexplorer2':pos_circexplorer2,'find_circ':pos_find_circ}
    if options['-s']:
        software=options['-s'].lower()
        if software not in soft_all:
            sys.exit('ERROR: %s is not included in our default software list, please try "-d" to designate the position' % options['-s'])
        circ_pos=pos_dict[software]
    else:
        circ_pos=checkD(options['-d'])
    circ_pos=[i-1 for i in circ_pos]
    # read file
    fileArr=[]
    if options['-n']:
        fileArr=options['-n'].split(',')
    else:
        if os.path.exists(options['-f']):
            fileArr=[i.strip() for i in open(options['-f']).readlines()]
        else:
            sys.exit('ERROR: %s is not exist!!!' % options['-f'])
    fileArr_2=list(set(fileArr))
    if len(fileArr_2)<len(fileArr):
        sys.exit('ERROR: file names is replicated!!!')
    fileCheck(fileArr)
    stime=time.time()
    print('Merge %d files.' % len(fileArr))
    if len(circ_pos)<5:
        pos_mat,circ_mat=comMat(fileArr,circ_pos)
    else:
        pos_mat,circ_mat,linear_mat=comMat(fileArr,circ_pos)
    if options['-s']:
        if software=='circexplorer2':
            pos_mat['start']=pos_mat['start'].map(lambda x:int(x)+1).tolist()
        elif software=='find_circ':
            pos_mat['start']=pos_mat['start'].map(lambda x:int(x)+1).tolist()
            pos_mat['end']=pos_mat['end'].map(lambda x:int(x)+1).tolist()
    if options['-d'] and options['-b']:
        pos_mat=adjPos(pos_mat,options['-b'])
    pos_mat.to_csv(output+'_pos.txt',sep='\t',header=True,index=None)
    circ_mat=pd.concat([pos_mat,circ_mat],axis=1)
    circ_mat.to_csv(output+'_circ.txt',sep='\t',header=True,index=None)
    if len(circ_pos)==5:
        linear_mat=pd.concat([pos_mat,linear_mat],axis=1)
        linear_mat.to_csv(output+'_linear.txt',sep='\t',header=True,index=None)
    etime=time.time()-stime
    print('Finished in %d second.' % etime)

#### function part ####
def readMatFile(fileName,pos):
    fff=open(fileName)
    fl=fff.readline().strip().split('\t')
    f2=fff.readline().strip().split('\t')
    fff.close()
    if len(f2) <(max(pos)+1):
        sys.exit('ERROR: The number of columns in %s is not enough!!!' % fileName)
    isHeader=None
    for i in fl[pos[1]]:
        if i not in ['0','1','2','3','4','5','6','7','8','9']:
            isHeader=1
            break
    if isHeader:
        mat=pd.read_csv(fileName,sep='\t',skiprows=[0],header=None).iloc[:,pos]
    else:
        mat=pd.read_csv(fileName,sep='\t',header=None).iloc[:,pos]
    return(mat)

def comMat(fileArr,pos):
    circ_dict={}
    linear_dict={}
    for i in fileArr:
        tmpMat=readMatFile(i,pos)
        for j in range(tmpMat.shape[0]):
            tmp=tmpMat.iloc[j,:].tolist()
            circID=str(tmp[0])+'|'+str(tmp[1])+'|'+str(tmp[2])
            if not circ_dict.__contains__(circID):
                circ_dict[circID]={}
            circ_dict[circID][i]=tmp[3]
            if len(pos)==5:
                if not linear_dict.__contains__(circID):
                    linear_dict[circID]={}
                linear_dict[circID][i]=tmp[4]
    circ_mat={}
    linear_mat={}
    pos_mat=[]
    for circID in circ_dict.keys():
        pos_mat.append(circID.split('|'))
        circ_mat[circID]=[]
        linear_mat[circID]=[]
        for i in fileArr:
            if circ_dict[circID].__contains__(i):
                circ_mat[circID].append(circ_dict[circID][i])
            else:
                circ_mat[circID].append(0)
            if len(pos)==5:
                if linear_dict[circID].__contains__(i):
                    linear_mat[circID].append(linear_dict[circID][i])
                else:
                    linear_mat[circID].append(0)
    circ_mat=list(circ_mat.values())
    circ_mat=pd.DataFrame(circ_mat,columns=fileArr)
    pos_mat=pd.DataFrame(pos_mat,columns=['chr','start','end'])
    #circ_mat=pd.concat([pos_mat,circ_mat],axis=1)
    if len(pos)<5:
        return(pos_mat,circ_mat)
    else:
        linear_mat=list(linear_mat.values())
        linear_mat=pd.DataFrame(linear_mat,columns=fileArr)
        #linear_mat=pd.concat([pos_mat,linear_mat],axis=1)
        return(pos_mat,circ_mat,linear_mat)