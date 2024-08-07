'''
Usage: DEBKS anno -c circ_pos -g gene_anno -m genome [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -c circ_pos                 File of circRNA position in tab format (the first three columns are chr, start, end).
    -g gene_anno                Tabix indexed gtf File of gene annotation.
    -m genome                   Fasta file of genome.
    -o output                   Output prefix [default: anno].
'''

import pysam,pandas as pd, numpy as np,time,os,sys,traceback

def fileCheck(file):
    for i in file:
        if not os.path.exists(i):
            sys.exit('ERROR: %s is not exist!!!' % i)

def readGTFfile(fileName):
    if not os.path.exists(fileName):
        sys.exit('ERROR: %s is not exist!!!' % fileName)
    try:
        tabixfile = pysam.TabixFile(fileName)
        return(tabixfile)
    except:
        sys.exit('ERROR: make sure %s is sorted and tabix indexed!!!' % fileName)

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

def readFaFile(fileName):
    if not os.path.exists(fileName):
        sys.exit('ERROR: %s is not exist!!!' % fileName)
    fff=open(fileName)
    fl=fff.readline()
    fff.close()
    if fl[0]!='>':
        sys.exit('ERROR:  %s need to be Fasta format!!!' % fileName)
    try:
        faFile=pysam.FastaFile(fileName)
        return(faFile)
    except:
        sys.exit('ERROR: make sure %s is fai indexed!!!' % fileName)

def anno(options):
    gtfFile=options['-g']
    faFile=options['-m']
    outAnno=options['-o']+'_detail.txt'
    outLen=options['-o']+'_len.txt'
    gtfFile=readGTFfile(gtfFile)
    faFile=readFaFile(faFile)
    circPos=readCircFile(options['-c']).iloc[:,0:3]
    circPos.columns=['chr','start','end']
    stime=time.time()
    print('Begin %s annotation for %d circRNAs.' % (options['-c'],circPos.shape[0]))
    annoPos=pd.DataFrame(circPos.apply(lambda x: pos2ano(x,gtfFile,faFile),axis=1).tolist(),columns=['type','strand','gene_id','gene_name','subtype','exon_start','exon_end','len'])
    annoPos_com=pd.concat([circPos,annoPos],axis=1)
    annoPos_com.to_csv(outAnno,index=False,sep='\t')
    annoPos['len']=annoPos['len'].map(str).tolist()
    lenDf0=annoPos['len'].map(lambda x: x.split(';')[0])
    lenDf=pd.concat([circPos,lenDf0],axis=1)
    lenDf.to_csv(outLen,index=False,sep='\t')
    etime=time.time()-stime
    print('Finished in %d seconds.' % etime)

##### function part ####
def getSeq(chr,start,end,faFile):
    try:
        return(faFile.fetch(chr,start,end))
    except:
        sys.exit('ERROR: make sure circ_pos and genome is match!!!' )

def getGene(chr,start,end,gtfFile):
    geneDict={}
    geneID_list=[]
    geneName_list=[]
    geneID_Dict={}
    geneStrand_list=[]
    if chr not in gtfFile.contigs:
        return(['intergenic'])
    for gtf in gtfFile.fetch(chr, start, end+1,parser=pysam.asGTF()):
        if gtf.feature=='gene':
            geneID_list.append(gtf.gene_id)
            geneName_list.append(gtf.gene_name)
            geneStrand_list.append(gtf.strand)
            geneID_Dict[gtf.gene_id]=[gtf.start,gtf.end]
        if gtf.feature=='exon':
            if geneDict.__contains__(gtf.transcript_id):
                geneDict[gtf.transcript_id].append([gtf.start,gtf.end,gtf.gene_id,gtf.gene_name,gtf.strand])
            else:
                geneDict[gtf.transcript_id]=[]
                geneDict[gtf.transcript_id].append([gtf.start,gtf.end,gtf.gene_id,gtf.gene_name,gtf.strand])
    if len(geneID_list)==0:
        return(['intergenic'])
    elif len(geneDict.keys())==0:
        return(['intronic',geneID_list,geneName_list,geneStrand_list])
    else:
        annoArr={}
        annoScore={'intergenic': 1, 'intronic': 3, 'inExon': 4, 'isExon': 7}
        for i in geneDict.keys():
            annoArr[i]=''
            tmpTranArr=geneDict[i]
            tmpExoArr_1=tmpTranArr[0]
            tmpGene=geneID_Dict[tmpExoArr_1[2]]
            if tmpExoArr_1[0]<start:
                annoArr[i]='inExon'
            elif tmpExoArr_1[0]==start:
                annoArr[i]='isExon'
            else:
                if tmpGene[0]>start:
                    annoArr[i]='intergenic'
                else:
                    annoArr[i]='intronic'
            tmpExoArr_L=tmpTranArr[-1]
            if tmpExoArr_L[1]<end:
                if tmpGene[1]<end:
                    annoArr[i]+='-intergenic'
                else:
                    annoArr[i]+='-intronic'
            elif tmpExoArr_L[1]==end:
                annoArr[i]+='-isExon'
            else:
                annoArr[i]+='-inExon'
        proTran=[]
        maxScore=0
        for i in annoArr.keys():
            tmpType=annoArr[i].split('-')
            tmpScore=annoScore[tmpType[0]]+annoScore[tmpType[1]]
            if tmpScore>maxScore:
                maxScore=tmpScore
                proTran=[i]
            elif tmpScore==maxScore:
                proTran.append(i)
            else:
                continue
        detailStart_arr=[]
        detailEnd_arr=[]
        detailLen_arr=[]
        detailGeneID_arr=[]
        detailGeneName_arr=[]
        detailAnno_arr=[]
        detailStrand_arr=[]
        for i in proTran:
            tmpS=[]
            tmpE=[]
            tmpTranArr=geneDict[i]
            lenExon=len(tmpTranArr)
            if lenExon==1:
                tmpS.append(start)
                tmpE.append(end)
            for j in range(lenExon):
                if j==0:
                    tmpS.append(start)
                    tmpE.append(tmpTranArr[j][1])
                elif j==lenExon-1:
                    tmpS.append(tmpTranArr[j][0])
                    tmpE.append(end)
                else:
                    tmpS.append(tmpTranArr[j][0])
                    tmpE.append(tmpTranArr[j][1])
            detailStart_arr.append(','.join([str(i+1) for i in tmpS]))
            detailEnd_arr.append(','.join([str(i) for i in tmpE]))
            detailLen_arr.append(sum(np.array(tmpE)-np.array(tmpS)))
            detailGeneID_arr.append(tmpTranArr[0][2])
            detailGeneName_arr.append(tmpTranArr[0][3])
            detailAnno_arr.append(annoArr[i])
            detailStrand_arr.append(tmpTranArr[0][4])
        return(['gene',
                list(set(zip(detailAnno_arr,detailGeneID_arr,detailGeneName_arr, detailStrand_arr,detailStart_arr,detailEnd_arr,detailLen_arr)))])
                
def getGenePair(geneID,geneName):
    geneZip=list(set(zip(geneID,geneName)))
    geneID=[]
    geneName=[]
    for i in geneZip:
        geneID.append(i[0])
        geneName.append(i[1])
    return(geneID,geneName)

def getGeneStrandPair(strand,geneID,geneName):
    geneZip=list(set(zip(strand,geneID,geneName)))
    strand=[]
    geneID=[]
    geneName=[]
    for i in geneZip:
        strand.append(i[0])
        geneID.append(i[1])
        geneName.append(i[2])
    return(strand,geneID,geneName)

def splitGenePair(x):
    geneID=[]
    geneName=[]
    for i in x:
        geneID.append(i[1])
        geneName.append(i[2])
    return(geneID,geneName)

def splitGeneAll(x,id):
    subType=[]
    geneID=[]
    geneName=[]
    start=[]
    end=[]
    circLen=[]
    strandList=[]
    for i in id:
        subType.append(x[i][0])
        geneID.append(x[i][1])
        geneName.append(x[i][2])
        strandList.append(x[i][3])
        start.append(x[i][4])
        end.append(x[i][5])
        circLen.append(str(x[i][6]))
    return(subType,geneID,geneName,strandList,start,end,circLen)

def subType2Type(x):
    topType=[]
    exonList=['isExon-inExon','inExon-isExon']
    mixList=['isExon-intronic','intronic-isExon','inExon-intronic','intronic-inExon',
             'isExon-intergenic','intergenic-isExon','inExon-intergenic','intergenic-inExon']
    for i in x:
        arr=list(set(i.split('-')))
        if len(arr)==1:
            if arr[0]=='intronic':
                topType.append('intronic')
            elif arr[0][2:]=='Exon':
                topType.append('exonic')
            else:
                topType.append('unknown')
        else:
            if i in exonList:
                topType.append('exon')
            elif i in mixList:
                topType.append('mix')
            else:
                topType.append('unknown')
    return(topType)

def adjGene(chr,start,end,gtfFile,strand):
    anno=getGene(chr,start,end,gtfFile)
    #type; strand; geneID; geneName; subType; start; end; len
    if anno[0]=='intergenic':
        return('intergenic',strand,'.','.','intergenic',start+1,end,end-start)
    elif anno[0]=='intronic':# type == intronic or intergenic; subType: antisense
        strand_list=anno[3]
        if strand in ['+','-']:
            keepID=[]
            for i in range(len(strand_list)):
                if strand_list[i]==strand:
                    keepID.append(i)
            if len(keepID)==0:
                geneID,geneName=getGenePair(anno[1],anno[2])
                return('intergenic',strand,';'.join(geneID),';'.join(geneName),'antisense',start+1,end,end-start)
            else:
                geneID=[anno[1][i] for i in keepID]
                geneName=[anno[2][i] for i in keepID]
                geneID,geneName=getGenePair(geneID,geneName)
                if (len(set(geneID))==1):
                    return('intronic',strand,geneID[0],geneName[0],'intronic',start+1,end,end-start)
                return('intronic',strand,';'.join(geneID),';'.join(geneName),'intronic',start+1,end,end-start)
        else:
            strand_list,geneID,geneName=getGeneStrandPair(strand_list,anno[1],anno[2])
            if (len(set(geneID))==1):
                return('intronic',strand_list[0],geneID[0],geneName[0],'intronic',start+1,end,end-start)
            return('intronic',';'.join(strand_list),';'.join(geneID),';'.join(geneName),'intronic',start+1,end,end-start)
    else:
        if strand in ['+','-']:
            keepID=[]
            for i in range(len(anno[1])):
                if anno[1][i][3]==strand:
                    keepID.append(i)
            if len(keepID)==0:
                geneID,geneName=splitGenePair(anno[1])
                if(len(set(geneID))==1):
                    return('intergenic',strand,geneID[0],geneName[0],'antisense',start+1,end,end-start)
                return('intergenic',strand,';'.join(geneID),';'.join(geneName),'antisense',start+1,end,end-start)
            else:
                subType,geneID,geneName,strand_list,start,end,circLen=splitGeneAll(anno[1],keepID)
                topType=subType2Type(subType)
                if(len(set(geneID))==1):
                    return(topType[0],strand,geneID[0],geneName[0],subType[0],';'.join(start),';'.join(end),';'.join(circLen))
                return(';'.join(topType),strand,';'.join(geneID),';'.join(geneName),';'.join(subType),';'.join(start),';'.join(end),';'.join(circLen))
        else:
            keepID=range(len(anno[1]))
            subType,geneID,geneName,strand_list,start,end,circLen=splitGeneAll(anno[1],keepID)
            topType=subType2Type(subType)
            if(len(set(geneID))==1):
                return(topType[0],strand_list[0],geneID[0],geneName[0],subType[0],';'.join(start),';'.join(end),';'.join(circLen))
            return(';'.join(topType),';'.join(strand_list),';'.join(geneID),';'.join(geneName),';'.join(subType),';'.join(start),';'.join(end),';'.join(circLen))

def pos2ano(pos,gtfFile,faFile):
    chr=pos[0]
    start=pos[1]-1
    end=pos[2]-1
    sSeq=getSeq(chr,start-2,start,faFile)
    eSeq=getSeq(chr,end+1,end+3,faFile)
    motifSeq=sSeq+eSeq
    strand='.'
    if motifSeq=='ACCT' or motifSeq=='GCCT' or motifSeq=='ATGT':
        strand='-'
    elif motifSeq=='AGGT' or motifSeq=='AGGC' or motifSeq=='ACAT':
        strand='+'
    return(list(adjGene(chr,start,end+1,gtfFile,strand)))
