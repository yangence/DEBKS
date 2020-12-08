'''
Usage: DEBKS dec -c circ (-l linear |--c2 circ2) [-p] [-n sample_num] [-f circRNA_filter] [-s sample_filter] [-t threads] [-d cutoff] [-r read_len] [-a hangover_len] [-e circ_len] [--e2 circ2_len] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -c circ                     Circular junction counts file in tab format, the first three columns is circRNA position (chr,start,end).
    -l linear                   Linear junction counts file, the format is same with -c.
    --c2 circ2                  Circular junction counts file in tab format, the format is same with -c.
    -p                          Sample group is paired.
    -n sample_num               Number of samples in each group, separated by comma (e.g. 2,3).
    -f circRNA_filter           Required total circular juction counts in all samples to filter out low expressed circRNAs [default: 0].
    -s sample_filter            Required lowest CJC of each circRNA to filter out sample [default: 0].
    -t threads                  Number of threads to detect DE circRNA [default: 4].
    -d cutoff                   Cutoff of DE circRNA in two groups [default: 0.05].
    -r read_len                 Length of RNA-seq reads [default: 150].
    -a hangover_len             The minimum hangover length for junction quantification [default: 6].
    -e circ_len                 File of circRNA length in tab format. The first four columns is circRNA position (chr, start, end) and circRNA length.
    --e2 circ2_len              File of circRNA2 length in tab format, only work with --c2. The format is same with -e.
    -o output                   Output file [default: dec_circRNA.txt].
'''
import sys,time,pandas as pd,numpy as np,docopt,os,time

def readJuncFile(fileName,reqLen=7):
    if not os.path.exists(fileName):
        sys.exit('ERROR: %s is not exist!!!' % fileName)
    fff=open(fileName)
    fl=fff.readline().strip().split('\t')
    fff.close()
    if len(fl) <reqLen:
        sys.exit('ERROR: The number of columns in %s is not enough!!!' % fileName)
    isHeader=None
    for i in fl[1]:
        if i not in ['0','1','2','3','4','5','6','7','8','9']:
            isHeader=0
            break
    juncMat=pd.read_csv(fileName,sep='\t',header=isHeader)
    return(juncMat)

def filterCirc(mat,filterNum):
    rowSum=mat.apply(lambda x: sum(x[3:]),axis=1).tolist()
    passID=[]
    for i in range(mat.shape[0]):
        if rowSum[i]>=filterNum:
            passID.append(i)
    return(passID)
    
def getCircLen(lenFile,circID,adjustReadLength):
    circLenFile=readJuncFile(lenFile,reqLen=4)
    circID3=(circLenFile.iloc[:,0]+':'+circLenFile.iloc[:,1].map(str)+'|'+circLenFile.iloc[:,2].map(str)).tolist()
    circLen=circLenFile.iloc[:,3].map(int).tolist()
    if len(circLen) != len(circID):
        sys.exit('ERROR: The number of lines in %s file (%d) is different with junction file (%d)!!!' % (lenFile,len(circLen),len(circID)))
    if circID !=circID3:
        sys.exit('ERROR: circRNA position is different in junction file and %s, DEBKS do not adjust order automatically!!!' % lenFile)
    if min(circLen)==0:
        sys.exit('ERROR: some circRNAs are designeted to 0 nt length in %s!!!' % lenFile)
    circLen=[int(i) for i in circLen]
    adjust_inc_len=[]
    adjust_skp_len=[]
    for i in circLen:
        if i<2*adjustReadLength:
            adjust_inc_len.append(2*adjustReadLength+i)
            adjust_skp_len.append(i)
        else:
            adjust_inc_len.append(4*adjustReadLength)
            adjust_skp_len.append(2*adjustReadLength)
    return(adjust_inc_len,adjust_skp_len)

def dec(options):
    #print(options)
    # number parameter
    threads=int(options['-t'])
    cutoff=min(float(options['-d']),1.00)
    readLength=int(options['-r'])
    hangOver=int(options['-a'])
    filterNum=max(int(options['-f']),0)
    filterSample=max(int(options['-s']),0)
    # is paired
    if options['-p']:
        from .rMATs_paired import getrMATS_P,vec2float,vecAddOne,adjustPsi,myFDR
        groupPair=True
    else:
        from .rMATs import getrMATS_P,vec2float,vecAddOne,adjustPsi,myFDR
        groupPair=False
    global getrMATS_P,vec2float,vecAddOne,adjustPsi,myFDR
    # junction file
    circFile=options['-c']
    circMat=readJuncFile(circFile)
    if options['-l']:
        linearFile=options['-l']
        linearMat=readJuncFile(linearFile)
    if options['--c2']:
        circ2File=options['--c2']
        linearMat=readJuncFile(circ2File)

    if circMat.shape[0]==linearMat.shape[0]:
        circID=(circMat.iloc[:,0]+':'+circMat.iloc[:,1].map(str)+'|'+circMat.iloc[:,2].map(str)).tolist()
        circID2=(linearMat.iloc[:,0]+':'+linearMat.iloc[:,1].map(str)+'|'+linearMat.iloc[:,2].map(str)).tolist()
        if not circID==circID2:
            sys.exit('ERROR: circRNA position is different in two junction file!!!')
    else:
        sys.exit('ERROR: The number of circRNA is different in two junction file!!!')
    # sample number
    sampleNum=circMat.shape[1]-3
    if options['-n']:
        group_num=options['-n'].split(',')
        group1_num=int(group_num[0])
        group2_num=int(group_num[1])
        if group1_num+group2_num !=sampleNum:
            sys.exit('ERROR: Sample number mismatch the junction file!!!')
    else:
        group1_num=int((circMat.shape[1]-3)/2)
        group2_num=sampleNum-group1_num
    if group1_num<2 or group2_num<2:
        sys.exit('ERROR: Require replicates!!!')
    # adjust inclusion and exclusion length
    adjustReadLength=readLength-hangOver
    if options['-e']:
        adjust_inc_len,adjust_skp_len=getCircLen(options['-e'],circID,adjustReadLength)
    else:
        adjust_inc_len=[4*adjustReadLength for i in range(circMat.shape[0])]
        adjust_skp_len=[2*adjustReadLength for i in range(circMat.shape[0])]
    if options['--c2']:
        if options['--e2']:
            tmp,adjust_inc_len=getCircLen(options['--e2'],circID2,adjustReadLength)
        else:
            adjust_inc_len=[2*adjustReadLength for i in range(circMat.shape[0])]
    
    if filterNum>0:
        passID=filterCirc(circMat,filterNum)
        if len(passID)>0:
            print('Keep %d circRNAs for dec.' % len(passID))
            circMat=circMat.iloc[passID,:].copy()
            linearMat=linearMat.iloc[passID,:].copy()
            adjust_inc_len=[adjust_inc_len[i] for i in passID]
            adjust_skp_len=[adjust_skp_len[i] for i in passID]
        else:
            sys.exit('WARNING: Filter keep zero circRNA for dec.')
    print('Begin detecting DE from %d circRNAs in %d vs %d samples.' % (circMat.shape[0],group1_num,group2_num))
    stime=time.time()
    decNum=stat(circMat,linearMat,group1_num,sampleNum,adjust_inc_len,adjust_skp_len,cutoff,groupPair,threads,options['-o'],options['--c2'],filterSample)
    etime=time.time()-stime
    print('Identify %d DEC with %f cutoff in %d seconds.' %(decNum,cutoff, etime))
##### functions part #############
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
def getMean(x):
    num=len(x)
    xSum=0
    for i in x:
        if i == 'nan':
            num=num-1
            continue
        try:
            tmp=float(i)
            xSum=xSum+tmp
        except:
            num=num-1
    if num ==0:
        return(0)
    else:
        return(xSum/num)
def get_delta(PBSI1,PBSI2):
    delta_pbsi=[]
    for i in range(len(PBSI1)):
        eachPbsi1_mean=getMean(PBSI1[i].split(','))
        eachPbsi2_mean=getMean(PBSI2[i].split(','))
        delta_pbsi.append(eachPbsi1_mean-eachPbsi2_mean)
    return(delta_pbsi)
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

def vec_remove_low_pair(inc1,skp1,inc2,skp2,num):
    res_inc1=[];res_skp1=[];res_inc2=[];res_skp2=[];
    for i in range(len(inc1)):
        if (float(skp1[i])>=num) & (float(skp2[i])>=num):
            res_inc1.append(inc1[i]);
            res_skp1.append(skp1[i]);
            res_inc2.append(inc2[i]);
            res_skp2.append(skp2[i]);
    return([res_inc1,res_skp1,res_inc2,res_skp2])
def vec_remove_low(inc,skp,num):
    res_inc=[];res_skp=[];
    for i in range(len(inc)):
        if float(skp[i])>=num:
            res_inc.append(inc[i]);
            res_skp.append(skp[i]);
    return([res_inc,res_skp]);

def checkInf(x):
    xt=[]
    for i in x:
        if i==float('inf') or i==float('-inf'):
            xt.append(float('nan'))
        else:
            xt.append(i)
    return(xt)
def stat(circMat,linearMat,group1_num,sampleNum,adjust_inc_len,adjust_skp_len,cutoff,groupPair,threads,output,isc2,filterSample):
    list_n_original_diff=[];psi_list_1=[];psi_list_2=[]
    for i in range(circMat.shape[0]):
        inc1=linearMat.iloc[i,3:(group1_num+3)].tolist()
        inc2=linearMat.iloc[i,(group1_num+3):(sampleNum+3)].tolist()
        skp1=circMat.iloc[i,3:(group1_num+3)].tolist()
        skp2=circMat.iloc[i,(group1_num+3):(sampleNum+3)].tolist()
        effective_inclusion_length=adjust_inc_len[i]
        effective_skipping_length=adjust_skp_len[i]
        inc1=vec2float(inc1);skp1=vec2float(skp1);inc2=vec2float(inc2);skp2=vec2float(skp2)
        inc1=checkInf(inc1);skp1=checkInf(skp1);inc2=checkInf(inc2);skp2=checkInf(skp2)
        if effective_inclusion_length <=0:
            effective_inclusion_length=readLength
        if groupPair:
            temp=vec_remove_na_pair(inc1,skp1,inc2,skp2);inc1_nona=temp[0];skp1_nona=temp[1];inc2_nona=temp[2];skp2_nona=temp[3];
            temp=vec_remove_low_pair(inc1_nona,skp1_nona,inc2_nona,skp2_nona,filterSample);inc1_nona=temp[0];skp1_nona=temp[1];inc2_nona=temp[2];skp2_nona=temp[3];
            if (len(inc1_nona)==0) | (len(inc2_nona)==0):
                psi_list_1.append('')
                psi_list_2.append('')
                list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,i]);
            else:
                psi_list_1.append(adjustPsi(inc1,skp1,effective_inclusion_length,effective_skipping_length))
                psi_list_2.append(adjustPsi(inc2,skp2,effective_inclusion_length,effective_skipping_length))
                inc1=inc1_nona;skp1=skp1_nona;inc2=inc2_nona;skp2=skp2_nona;
                inc1=vecAddOne(inc1);skp1=vecAddOne(skp1);inc2=vecAddOne(inc2);skp2=vecAddOne(skp2);
                if (len(inc1)==1) |(len(inc2)==1):
                    list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,i]);
                else:
                    list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,1,i]);
        else:
            temp1=vec_remove_na(inc1,skp1);temp2=vec_remove_na(inc2,skp2);inc1_nona=temp1[0];skp1_nona=temp1[1];inc2_nona=temp2[0];skp2_nona=temp2[1];
            temp1=vec_remove_low(inc1_nona,skp1_nona,filterSample);temp2=vec_remove_low(inc2_nona,skp2_nona,filterSample);inc1_nona=temp1[0];skp1_nona=temp1[1];inc2_nona=temp2[0];skp2_nona=temp2[1];
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
    FDR=myFDR(P)
    pbsi1=psi2PBSI(psi_list_1)
    pbsi2=psi2PBSI(psi_list_2)
    deltPBSI=get_delta(pbsi1,pbsi2)
    cjc1=circMat.iloc[:,3].map(str)
    cjc2=circMat.iloc[:,(group1_num+3)].map(str)
    ljc1=linearMat.iloc[:,3].map(str)
    ljc2=linearMat.iloc[:,(group1_num+3)].map(str)
    for i in range(4,group1_num+3):
        cjc1=cjc1+','+circMat.iloc[:,i].map(str)
    for i in range(group1_num+4,sampleNum+3):
        cjc2=cjc2+','+circMat.iloc[:,i].map(str)
    for i in range(4,group1_num+3):
        ljc1=ljc1+','+linearMat.iloc[:,i].map(str)
    for i in range(group1_num+4,sampleNum+3):
        ljc2=ljc2+','+linearMat.iloc[:,i].map(str)
    if isc2:
        matStat=pd.DataFrame({'chr':circMat.iloc[:,0],'start':circMat.iloc[:,1],'end':circMat.iloc[:,2],
        'cjc_1':cjc1.tolist(),'cjc_2':cjc2.tolist(),'cjc2_1':ljc1.tolist(),
        'cjc2_2':ljc2.tolist(),'adj_cjc_len':adjust_skp_len,'adj_cjc2_len':adjust_inc_len,
        'pbsi1':pbsi1,'pbsi2':pbsi2,'delta_pbsi':deltPBSI,'P':P,'FDR':FDR})
        matStat.to_csv(output,index=False,sep='\t',header=True)
    else:
        matStat=pd.DataFrame({'chr':circMat.iloc[:,0],'start':circMat.iloc[:,1],'end':circMat.iloc[:,2],
        'cjc_1':cjc1.tolist(),'cjc_2':cjc2.tolist(),'ljc_1':ljc1.tolist(),
        'ljc_2':ljc2.tolist(),'adj_cjc_len':adjust_skp_len,'adj_ljc_len':adjust_inc_len,
        'pbsi1':pbsi1,'pbsi2':pbsi2,'delta_pbsi':deltPBSI,'P':P,'FDR':FDR})
        matStat.to_csv(output,index=False,sep='\t',header=True)
    return(sum(np.array(FDR)<0.05))
