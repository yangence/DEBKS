import pandas as pd,random, numpy as np,math,sys,os
np.random.seed(10)

sigma=float(sys.argv[1])
outFile=sys.argv[2]
if outFile[-1]!='/':
    outFile=outFile+'/'
#gtfFile='/media/data4/lzl/annotation/GENCODE/humanV19/gencode.v19.annotation.gene.gtf'
gtfFile=sys.argv[3]
gtfDf=pd.read_csv(gtfFile,sep='\t',header=None)
geneID=list(set(gtfDf.iloc[:,8].map(lambda x: x.split('; ')[0]).tolist()))

#geneCovDf=pd.read_csv('/media/data4/lzl/DEBKS/script/outJunction_mouse.txt',sep='\t',names=['c1','c2','s1','s2'])
geneCovDf=pd.read_csv(sys.argv[4],sep='\t',names=['c1','c2','s1','s2'])
choizeGeneCov= np.random.choice(geneCovDf.shape[0], size=len(geneID), replace=True)

percentDiff=0.5
cutoff=0.05
groupNum1=3
groupNum2=3

def getCommonPsi_beta(geneID,cutoff):
    majorNum=int(len(geneID)*0.95)
    minorNum=len(geneID)-majorNum
    psi_value1=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
    psi_value2=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
    psi_g1=[]
    psi_g2=[]
    diffNum=0
    normalNum=len(geneID)
    while len(psi_g1)!=len(geneID):
        for i in range(len(geneID)):
            p1=psi_value1[i]
            p2=psi_value2[i]
            if abs(p1-p2)<=cutoff:
                if normalNum!=0:
                    psi_g1.append(p1)
                    psi_g2.append(p2)
                    normalNum+=-1
            else:
                if diffNum!=0:
                    psi_g1.append(p1)
                    psi_g2.append(p2)
                    diffNum+=-1
        psi_value1=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
        psi_value2=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
    return(psi_g1,psi_g2)

def getChangePsi_beta(geneID,cutoff,percentDiff):
    majorNum=int(len(geneID)*0.95)
    minorNum=len(geneID)-majorNum
    psi_value1=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
    psi_value2=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
    psi_g1=[]
    psi_g2=[]
    diffNum=len(geneID)
    normalNum=0
    while len(psi_g1)!=len(geneID):
        for i in range(len(geneID)):
            p1=psi_value1[i]
            p2=psi_value2[i]
            if abs(p1-p2)<=cutoff:
                if normalNum!=0:
                    psi_g1.append(p1)
                    psi_g2.append(p2)
                    normalNum+=-1
            else:
                if diffNum!=0:
                    psi_g1.append(p1)
                    psi_g2.append(p2)
                    diffNum+=-1
        psi_value1=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
        psi_value2=np.random.permutation(list(np.random.beta(1,12,majorNum))+list(np.random.beta(7,10,minorNum)))
    return(psi_g1,psi_g2)


def getCommonPsi(geneID,cutoff):
    psi_value=np.random.rand(len(geneID))
    psi_g1=[]
    # in cutoff
    for i in psi_value:
        tmp=np.random.rand()*cutoff *pow(-1,np.random.randint(2))
        while not (tmp+i>=0 and tmp+i<=1):
            tmp=np.random.rand()/(1/cutoff) *pow(-1,np.random.randint(2))
        psi_g1.append(tmp+i)
    return(psi_value,psi_g1)

def getChangePsi(geneID,cutoff,percentDiff):
    psi_value=np.random.rand(len(geneID))
    psi_g1=[]
    # in cutoff
    for i in psi_value:
        tmp=np.random.rand()/(1/cutoff) *pow(-1,np.random.randint(2))
        while not (tmp+i>0 and tmp+i<1):
            tmp=np.random.rand()/(1/cutoff) *pow(-1,np.random.randint(2))
        psi_g1.append(tmp+i)

    #choiceMore=np.random.choice(len(geneID),int(len(geneID)*percentDiff), replace=False)
    choiceMore=range(len(geneID))
    # out cutoff
    for i in choiceMore:
        tmp=(np.random.rand()*(1-cutoff)+cutoff)*pow(-1,np.random.randint(2))
        while not (tmp+psi_value[i]>0 and tmp+psi_value[i]<1):
            tmp=(np.random.rand()*(1-cutoff)+cutoff)*pow(-1,np.random.randint(2))
        psi_g1[i]=psi_value[i]+tmp
    psi_g2=psi_value.tolist()
    return(psi_g1,psi_g2)

def getRepPsi(psi,sigma):
    rep=[]
    for i in psi:
        tmp=np.random.normal(0,sigma)
        while not (tmp+i>=0 and tmp+i<=1):
            tmp=np.random.normal(0,sigma)
        rep.append(tmp+i)
    return(rep)
    
def getJun(junM,junS):
    tmp=np.random.normal(junM,junS)
    while tmp<=0:
        tmp=np.random.normal(junM,junS)
    return(tmp)

def getCircLinear(psi,gene):
    gene=np.array(gene)
    circ=psi*gene
    linear=gene-circ
    return([int(round(i)) for i in circ],[int(round(i)) for i in linear])

def getEachJunction(psi1,psi2,jun):
    if psi1 == 0 or psi2==0 or psi1 ==1 or psi2==1:
        psi3=psi1/(psi1+psi2)
        return(psi3*jun,(1-psi3)*jun,0)
    linear=jun/(1+(psi1/(1-psi1))+(psi2/(1-psi2)))
    circ1=linear*(psi1/(1-psi1))
    circ2=linear*(psi2/(1-psi2))
    return(int(circ1),int(circ2),int(linear))
def getJunction(psi1,psi2,jun):
    l=[]
    c1=[]
    c2=[]
    for i in range(len(psi1)):
        tmp1,tmp2,tmp3=getEachJunction(psi1[i],psi2[i],jun[i])
        c1.append(tmp1);c2.append(tmp2);l.append(tmp3)
    return(c1,c2,l)
    
# total counts
geneGroup1_rep1=[]
geneGroup1_rep2=[]
geneGroup1_rep3=[]
geneGroup2_rep1=[]
geneGroup2_rep2=[]
geneGroup2_rep3=[]
for i in choizeGeneCov:
    geneGroup1_rep1.append(int(getJun(geneCovDf.iloc[i,0],geneCovDf.iloc[i,2])))
    geneGroup1_rep2.append(int(getJun(geneCovDf.iloc[i,0],geneCovDf.iloc[i,2])))
    geneGroup1_rep3.append(int(getJun(geneCovDf.iloc[i,0],geneCovDf.iloc[i,2])))
    geneGroup2_rep1.append(int(getJun(geneCovDf.iloc[i,1],geneCovDf.iloc[i,3])))
    geneGroup2_rep2.append(int(getJun(geneCovDf.iloc[i,1],geneCovDf.iloc[i,3])))
    geneGroup2_rep3.append(int(getJun(geneCovDf.iloc[i,1],geneCovDf.iloc[i,3])))

# circRNA1 without changes
c1_psi_g1,c1_psi_g2=getCommonPsi_beta(geneID,cutoff)
c1_psi_g1_rep1=getRepPsi(c1_psi_g1,sigma)
c1_psi_g1_rep2=getRepPsi(c1_psi_g1,sigma)
c1_psi_g1_rep3=getRepPsi(c1_psi_g1,sigma)
c1_psi_g2_rep1=getRepPsi(c1_psi_g2,sigma)
c1_psi_g2_rep2=getRepPsi(c1_psi_g2,sigma)
c1_psi_g2_rep3=getRepPsi(c1_psi_g2,sigma)

# circRNA2 with changes 
c2_psi_g1,c2_psi_g2=getChangePsi_beta(geneID,cutoff,percentDiff)
c2_psi_g1_rep1=getRepPsi(c2_psi_g1,sigma)
c2_psi_g1_rep2=getRepPsi(c2_psi_g1,sigma)
c2_psi_g1_rep3=getRepPsi(c2_psi_g1,sigma)
c2_psi_g2_rep1=getRepPsi(c2_psi_g2,sigma)
c2_psi_g2_rep2=getRepPsi(c2_psi_g2,sigma)
c2_psi_g2_rep3=getRepPsi(c2_psi_g2,sigma)

# linear junction
c1_g1_rep1,c2_g1_rep1,l_g1_rep1=getJunction(c1_psi_g1_rep1,c2_psi_g1_rep1,geneGroup1_rep1)
c1_g1_rep2,c2_g1_rep2,l_g1_rep2=getJunction(c1_psi_g1_rep2,c2_psi_g1_rep2,geneGroup1_rep2)
c1_g1_rep3,c2_g1_rep3,l_g1_rep3=getJunction(c1_psi_g1_rep3,c2_psi_g1_rep3,geneGroup1_rep3)

c1_g2_rep1,c2_g2_rep1,l_g2_rep1=getJunction(c1_psi_g2_rep1,c2_psi_g2_rep1,geneGroup2_rep1)
c1_g2_rep2,c2_g2_rep2,l_g2_rep2=getJunction(c1_psi_g2_rep2,c2_psi_g2_rep2,geneGroup2_rep2)
c1_g2_rep3,c2_g2_rep3,l_g2_rep3=getJunction(c1_psi_g2_rep3,c2_psi_g2_rep3,geneGroup2_rep3)


def outputDf(geneID,circCov1,circCov2,linearCov,c1_psi,c2_psi,c1_g_psi,c2_g_psi,group='g1',rep='rep1',sigma=0):
    covDf=pd.DataFrame({'geneID':geneID,'circ1':circCov1,'circ2':circCov2,
                        'linear':linearCov,'c1_psi':c1_psi,'c2_psi':c2_psi,'c1_g_psi':c1_g_psi,'c2_g_psi':c2_g_psi})
    covDf.iloc[:,0]=covDf.iloc[:,0].map(lambda x:x[9:-1])
    covDf.to_csv(outFile+"coverage_%s_%s_%s.txt" % (str(sigma),group,rep),sep='\t',header=None,index=None)
    
outputDf(geneID,c1_g1_rep1,c2_g1_rep1,l_g1_rep1,c1_psi_g1_rep1,c2_psi_g1_rep1,c1_psi_g1,c2_psi_g1,group='g1',rep='rep1',sigma=sigma)
outputDf(geneID,c1_g1_rep2,c2_g1_rep2,l_g1_rep2,c1_psi_g1_rep2,c2_psi_g1_rep2,c1_psi_g1,c2_psi_g1,group='g1',rep='rep2',sigma=sigma)
outputDf(geneID,c1_g1_rep3,c2_g1_rep3,l_g1_rep3,c1_psi_g1_rep3,c2_psi_g1_rep3,c1_psi_g1,c2_psi_g1,group='g1',rep='rep3',sigma=sigma)
outputDf(geneID,c1_g2_rep1,c2_g2_rep1,l_g2_rep1,c1_psi_g2_rep1,c2_psi_g2_rep1,c1_psi_g2,c2_psi_g2,group='g2',rep='rep1',sigma=sigma)
outputDf(geneID,c1_g2_rep2,c2_g2_rep2,l_g2_rep2,c1_psi_g2_rep2,c2_psi_g2_rep2,c1_psi_g2,c2_psi_g2,group='g2',rep='rep2',sigma=sigma)
outputDf(geneID,c1_g2_rep3,c2_g2_rep3,l_g2_rep3,c1_psi_g2_rep3,c2_psi_g2_rep3,c1_psi_g2,c2_psi_g2,group='g2',rep='rep3',sigma=sigma)

