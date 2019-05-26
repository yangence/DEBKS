############################
#author: Zelin Liu
#email: zlliu@bjmu.edu.cn
#license: GPL3
#detail: Plot differentially expressed back-splicing based on the call of DEBKS_plot.py
#############################
options(stringsAsFactors = F)
"%&%"=function(a,b)paste0(a,b)

# required args
args<-commandArgs(T)

#temp
#args=c('chr11', 581491, 582081, '+', '115,120', '0,470', 'ENSG00000070047.7', 'ENST00000534320.1', 576592, 587258, 576485, 587464, '41,24,49,49,4,21', '37,60,66,36,30,21', '81,60,93,103,17,46', '70,118,104,73,30,43', '4.0,0.0,3.0,2.0,0.0,1.0', '0.0,3.0,1.0,2.0,0.0,0.0', '0.9101123595505618,1.0,0.9393939393939394,0.9626168224299065,1.0,0.9583333333333333', '1.0,0.9516129032258065,0.9811320754716982,0.9480519480519481,1.0,1.0', '98503656,78159750,114271402,115552918,148935318,152134018', '106029948,128533350,126735728,112039666,67850850,223838280', '/media/data4/lzl/DEBKS/data/BCresult/DEBKS_plot/temp/chr11_581491_582081.depth.txt', '/media/data4/lzl/DEBKS/data/BCresult/DEBKS_plot/chr11_581491_582081.pdf')
#
chr=args[1]
start=as.numeric(args[2])
end=as.numeric(args[3])
strand=args[4]
exonSize=as.numeric(unlist(strsplit(args[5],',')))
exonOff=as.numeric(unlist(strsplit(args[6],',')))
geneName=args[7]
isoName=args[8]
intronL=as.numeric(args[9])
intronR=as.numeric(args[10])
linearL=as.numeric(args[11])
linearR=as.numeric(args[12])
SJL1=as.numeric(unlist(strsplit(args[13],',')))
SJL2=as.numeric(unlist(strsplit(args[14],',')))
SJL=c(SJL1,SJL2)
SJR1=as.numeric(unlist(strsplit(args[15],',')))
SJR2=as.numeric(unlist(strsplit(args[16],',')))
SJR=c(SJR1,SJR2)
BS1=as.numeric(unlist(strsplit(args[17],',')))
BS2=as.numeric(unlist(strsplit(args[18],',')))
BS=c(BS1,BS2)
psi1=as.numeric(unlist(strsplit(args[19],',')))
psi2=as.numeric(unlist(strsplit(args[20],',')))
psi=round(c(psi1,psi2),2)
mapRead1=as.numeric(unlist(strsplit(args[21],',')))
mapRead2=as.numeric(unlist(strsplit(args[22],',')))
mapRead=c(mapRead1,mapRead2)
fileName=args[23]
outFile=args[24]
matDepth=read.delim(fileName,header=F)
matDepth=t(t(matDepth)/mapRead)*10^9
# 
id=chr%&%":"%&%start%&%"|"%&%end
sampleNum1=length(psi1)
sampleNum2=length(psi2)
sampleNum=sampleNum1+sampleNum2

#functions
plotBox=function(xstart,size,maxvalue,height,strand){

  for(i in 1:length(xstart)){
    rect(xleft = xstart[i],xright = xstart[i]+size[i]-1,ybottom = 0,ytop = height)
  }
  if(length(xstart)==1){
    return(0)
  }
  for(i in 1:(length(xstart)-1)){
    lines(c(xstart[i]+size[i]-1,xstart[i+1]),c(height/2,height/2),lwd=2)
    plotArrow(xstart[i]+size[i]-1,xstart[i+1],maxvalue,height/2,strand)
  }
}
plotArrow=function(start,end,maxvalue,height,strand){
  stepSize=maxvalue/100
  arrowSize=height/10
  if(start+2*stepSize>end){
      return(0)
  }
  if(strand=='+'){
     for(i in seq(start+stepSize,end-stepSize,stepSize)){
      lines(c(i,i-stepSize/2),c(height,height+arrowSize))
      lines(c(i,i-stepSize/2),c(height,height-arrowSize))
    }
  }else{
    for(i in seq(start+stepSize,end-stepSize,stepSize)){
      lines(c(i,i+stepSize/2),c(height,height+arrowSize))
      lines(c(i,i+stepSize/2),c(height,height-arrowSize))
    }
  }
 
  
}
# plot parameter
{
pdf(outFile,width = 10, height = (sampleNum+1)*2)
par(mar=c(2,4,4,2),mfrow=c(sampleNum+1,1))
for(i in 1:sampleNum){
  if(i<=sampleNum1){
    col=rgb(230,15,15,maxColorValue = 255)
    title_sample='sample_1_rep_'%&%i%&%' PBSI:'%&%' '%&%psi[i]
  }else{
    col=rgb(15,15,230,maxColorValue = 255)
    title_sample='sample_2_rep_'%&%(i-sampleNum1)%&%' PBSI:'%&%' '%&%psi[i]
  }
  plot(matDepth[,i],type='h',xaxt="n",bty='n',xaxs='i',col=col,ylab='RPkM',las=1,xlab='',main=title_sample)
    
  
  if(linearL == 0){
    indexBS=8
    indexBE=end-start-9
    heightS=matDepth[indexBS,i]
    if(linearR ==0){
      break()
    }else{
      indexR=intronR-start+8
      heightR=matDepth[indexR,i]
      heightE=matDepth[indexBE,i]
      midR=(indexBE+indexR)/2
      maxheight=max(matDepth[indexR:nrow(matDepth),i])
      lines(c(indexBE,midR),c(heightE,maxheight),col="#BEBADA",lwd=2)
      lines(c(midR,indexR),c(maxheight,heightR),col="#BEBADA",lwd=2)
      text(x=midR,y=maxheight*0.9,labels = SJR[i])
    }
  }else{
    indexBS=start-linearL+8
    indexBE=end-linearL-9
    indexL=intronL-linearL-9
    heightL=matDepth[indexL,i]
    heightS=matDepth[indexBS,i]
    midL=(indexL+indexBS)/2
    maxheight=max(matDepth[1:indexL,i])
    lines(c(indexL,midL),c(heightL,maxheight),col="#BEBADA",lwd=2)
    lines(c(midL,indexBS),c(maxheight,heightS),col="#BEBADA",lwd=2)
    text(x=midL,y=maxheight*0.9,labels = SJL[i])
      if(linearR ==0){

      }else{
        indexR=intronR-linearL+8
        heightR=matDepth[indexR,i]
        heightE=matDepth[indexBE,i]
        midR=(indexBE+indexR)/2
        maxheight=max(matDepth[indexR:nrow(matDepth),i])
        lines(c(indexBE,midR),c(heightE,maxheight),col="#BEBADA",lwd=2)
        lines(c(midR,indexR),c(maxheight,heightR),col="#BEBADA",lwd=2)
        text(x=midR,y=maxheight*0.9,labels = SJR[i])
      }
    
  }
  midSE=(indexBS+indexBE)/2

  maxheight=max(matDepth[indexBS:indexBE,i])
  text(x=midSE,y=maxheight*0.9,labels = BS[i])
  r=midSE-indexBS
  x=r*cos(seq(0,50,0.5)/100*2*pi)+midSE
  y=c((maxheight-heightE)*sin(seq(0,24,0.5)/100*2*pi)+heightE,(maxheight-heightS)*sin(seq(24.5,50,0.5)/100*2*pi)+heightS)
  lines(x,y,lwd=2,col="#8DD3C7")
}
index=c(1,floor(nrow(matDepth)/4),floor(nrow(matDepth)/2),floor(nrow(matDepth)*3/4),nrow(matDepth))
axis(1,at = (1:nrow(matDepth))[index],labels = (linearL:(linearL+nrow(matDepth)-1))[index],line=1)
#plot(matDepth[,i],ylim=c(0,100),type='n',xaxs='i',xaxt='n',xlab=geneName%&%" ("%&%strand%&%") "%&%id,bty='n',yaxt='n')
plot(matDepth[,i],ylim=c(0,100),type='n',xaxt="n",yaxt='n',bty='n',xaxs='i',col=col,ylab='',las=1,main=geneName%&%" ("%&%strand%&%") "%&%id,xlab='')
maxheight=100
plotBox(exonOff+start+1-linearL,exonSize,maxvalue =  nrow(matDepth),maxheight,strand)
if(linearL != 0 & linearR !=0){
  rect(xleft = 1,xright = intronL-linearL+1,ybottom = 0,ytop = maxheight)
  rect(xleft = intronR-linearL+1,xright = linearR-linearL+1,ybottom = 0,ytop = maxheight)
  lines(c(intronL-linearL+2,start-linearL+1),c(maxheight/2,maxheight/2),lwd=2)
  plotArrow(start=intronL-linearL+1,end=start-linearL+1,maxvalue =  nrow(matDepth),height=maxheight/2,strand)
  lines(c(end-linearL+1,intronR-linearL+1),c(maxheight/2,maxheight/2),lwd=2)
  plotArrow(start=end-linearL+1,end=intronR-linearL+1,maxvalue =  nrow(matDepth),height=maxheight/2,strand)
}else if(linearL == 0){
  rect(xleft = intronR-start+1,xright = linearR-start+1,ybottom = 0,ytop = maxheight)
  
  lines(c(end-start+1,intronR-start+1),c(maxheight/2,maxheight/2),lwd=2)
  plotArrow(end-start+1,intronR-start+1,maxvalue =  nrow(matDepth),height=maxheight/2,strand)
  
}else{
  rect(xleft = 1,xright = intronL-linearL+1,ybottom = 0,ytop = maxheight)
  lines(c(intronL-linearL+2,start-linearL+1),c(maxheight/2,maxheight/2),lwd=2)
  plotArrow(intronL-linearL+2,start-linearL+1,maxvalue =  nrow(matDepth),height=maxheight/2,strand)
}

dev.off()
}
