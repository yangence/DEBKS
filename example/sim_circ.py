import pandas as pd,random, numpy as np,math,sys
np.random.seed(10)

sigma=float(sys.argv[1])
outFile=sys.argv[2]
if outFile[-1]!='/':
    outFile=outFile+'/'
gtfFile=sys.argv[3]
gtfDf=pd.read_csv(gtfFile,sep='\t',header=None)
geneID=list(set(gtfDf.iloc[:,8].map(lambda x: x.split('; ')[0]).tolist()))

psi_value=np.random.rand(len(geneID))
psi_g1=[]
percentDiff=0.5
for i in psi_value:
    tmp=np.random.rand()/20 *pow(-1,np.random.randint(2))
    while not (tmp+i>=0 and tmp+i<=1):
        tmp=np.random.rand()/20
    psi_g1.append(tmp+i)

choiceMore=np.random.choice(len(geneID),int(len(geneID)*percentDiff), replace=False)

for i in choiceMore:
    tmp=(np.random.rand()*95+5)/100*pow(-1,np.random.randint(2))
    while not (tmp+psi_value[i]>=0 and tmp+psi_value[i]<=1):
        tmp=(np.random.rand()*95+5)/100*pow(-1,np.random.randint(2))
    psi_g1[i]=psi_value[i]+tmp

psi_g2=psi_value

def getRepPsi(psi,sigma):
    rep=[]
    for i in psi:
        tmp=np.random.normal(0,sigma)
        while not (tmp+i>=0 and tmp+i<=1):
            tmp=np.random.normal(0,sigma)
        rep.append(tmp+i)
    return(rep)
    
def getRepGene(gene,sigma):
    rep=abs(gene+np.random.normal(0,sigma,len(gene)))
    return(rep)
    

psi_g1_rep1=getRepPsi(psi_g1,sigma)
psi_g1_rep2=getRepPsi(psi_g1,sigma)
psi_g1_rep3=getRepPsi(psi_g1,sigma)
psi_g2_rep1=getRepPsi(psi_g2,sigma)
psi_g2_rep2=getRepPsi(psi_g2,sigma)
psi_g2_rep3=getRepPsi(psi_g2,sigma)


geneCovDf=pd.read_csv('outJunction_mouse.txt',sep='\t',names=['c1','c2','s1','s2'])
choizeGeneCov= np.random.choice(geneCovDf.shape[0], size=len(geneID), replace=True)

def getJun(junM,junS):
    tmp=np.random.normal(junM,junS)
    while tmp<=0:
        tmp=np.random.normal(junM,junS)
    return(tmp)

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


def getCircLinear(psi,gene):
    gene=np.array(gene)
    circ=psi*gene
    linear=gene-circ
    return([[int(round(i)) for i in circ],[int(round(i)) for i in linear]])

[circ_g1_rep1,linear_g1_rep1]=getCircLinear(psi_g1_rep1,geneGroup1_rep1)
[circ_g1_rep2,linear_g1_rep2]=getCircLinear(psi_g1_rep2,geneGroup1_rep2)
[circ_g1_rep3,linear_g1_rep3]=getCircLinear(psi_g1_rep3,geneGroup1_rep3)

[circ_g2_rep1,linear_g2_rep1]=getCircLinear(psi_g2_rep1,geneGroup2_rep1)
[circ_g2_rep2,linear_g2_rep2]=getCircLinear(psi_g2_rep2,geneGroup2_rep2)
[circ_g2_rep3,linear_g2_rep3]=getCircLinear(psi_g2_rep3,geneGroup2_rep3)


def outputDf(geneID,circCov,linearCov,psi,epsi,group='g1',rep='rep1',sigma=0):
    covDf=pd.DataFrame({'geneID':geneID,'circ':circCov,'gene':linearCov,'groupPsi':psi,'eachPsi':epsi})
    covDf.iloc[:,0]=covDf.iloc[:,0].map(lambda x:x[9:-1])
    covDf.to_csv(outFile+"coverage_%s_%s_%s.txt" % (str(sigma),group,rep),sep='\t',header=None,index=None)


outputDf(geneID,circ_g1_rep1,linear_g1_rep1,psi=psi_g1,epsi=psi_g1_rep1,group='g1',rep='rep1',sigma=sigma)
outputDf(geneID,circ_g1_rep2,linear_g1_rep2,psi=psi_g1,epsi=psi_g1_rep2,group='g1',rep='rep2',sigma=sigma)
outputDf(geneID,circ_g1_rep3,linear_g1_rep3,psi=psi_g1,epsi=psi_g1_rep3,group='g1',rep='rep3',sigma=sigma)
outputDf(geneID,circ_g2_rep1,linear_g2_rep1,psi=psi_g2,epsi=psi_g2_rep1,group='g2',rep='rep1',sigma=sigma)
outputDf(geneID,circ_g2_rep2,linear_g2_rep2,psi=psi_g2,epsi=psi_g2_rep2,group='g2',rep='rep2',sigma=sigma)
outputDf(geneID,circ_g2_rep3,linear_g2_rep3,psi=psi_g2,epsi=psi_g2_rep3,group='g2',rep='rep3',sigma=sigma)

