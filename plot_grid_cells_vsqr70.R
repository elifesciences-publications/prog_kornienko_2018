t<-read.table(paste(ep@resultsDirectory,"cells.class.table",sep="/"),sep="\t")
t2<-read.table(paste(ep@resultsDirectory,"cells.sig.changes",sep="/"),sep="\t")

load(paste(ep@resultsDirectory,"grid.pairs.id",sep="/"))
load(paste(ep@resultsDirectory,"grid.pairs.map",sep="/"))
load(paste(ep@resultsDirectory,"grid.pairs.ifra",sep="/"))

#########################################################################################

cellid = cells1$cell.id[t$grid==1]
source(paste(indir,'/plot_example_cells_grid.R',sep=""))

x = (round(sum(t$grid==1 & (t2$mapcor==1 | t2$fr==1 | t2$hd.dir==1))/sum(t$grid==1)*100))
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in map similarity, rate or direction: ",x,"%",sep=""))
x = round(sum(t$grid==1 & t2$mapcor==1 & t2$fr!=1 & t2$hd.dir!=1)/sum(t$grid==1)*100)
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in map similarity: ",x,"%",sep=""))
x = round(sum(t$grid==1 & t2$mapcor==1 & t2$fr==1 & t2$hd.dir!=1)/sum(t$grid==1)*100)
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in map similarity and rate: ",x,"%",sep=""))
x = round(sum(t$grid==1 & t2$mapcor==1 & t2$fr!=1 & t2$hd.dir==1)/sum(t$grid==1)*100)
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in map similarity and direction: ",x,"%",sep=""))
x = round(sum(t$grid==1 & t2$mapcor==1 & t2$fr==1 & t2$hd.dir==1)/sum(t$grid==1)*100)
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in map similarity, rate and direction: ",x,"%",sep=""))
x = round(sum(t$grid==1 & t2$mapcor!=1 & t2$fr==1 & t2$hd.dir!=1)/sum(t$grid==1)*100)
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in rate: ",x,"%",sep=""))
x = round(sum(t$grid==1 & t2$mapcor!=1 & t2$fr==1 & t2$hd.dir==1)/sum(t$grid==1)*100)
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in rate and direction: ",x,"%",sep=""))
x = round(sum(t$grid==1 & t2$mapcor!=1 & t2$fr!=1 & t2$hd.dir==1)/sum(t$grid==1)*100)
print(paste("Percentage of N = ",sum(t$grid==1), " grid cells with significant changes in direction: ",x,"%",sep=""))


# grid cells changing significantly
# ids=cells1$cell.id[(t$grid==1 & (t2$mapcor==1 | t2$fr==1 | t2$hd.dir==1))]
# ids = id[,1]%in%ids | id[,2]%in%ids 

cexn=0.6
par(mfrow=c(2,3))
for (i in 1:2){
  # if(i==1){x=grid.pairs.ifra[,ids==F]}
  # if(i==2){x=grid.pairs.map[,ids==F]}
  if(i==1){x=grid.pairs.ifra;ylab=c(0,.05,.1)}
  if(i==2){x=grid.pairs.map;ylab=c(0,.05,.1)}
  if(i==1){
    plot(x[1,],x[3,],cex=cexn,pch=20,axes=F,xlim=c(-.6,.6),ylim=c(-.6,.6),
    main = paste("Between conditions, r=",round(cor(x[1,],x[3,]),2),sep = ""),
    xlab = "IFR association vp1",ylab = "IFR association vp2")
    axis(side=1,tck=-.05,pos=-.6, at=c(-.6,0,.6))
    axis(side=2,tck=-.05,pos=-.6, at=c(-.6,0,.6),las=1)
  }
  
  if(i==2){
    plot(x[1,],x[3,],cex=cexn,pch=20,axes=F,xlim=c(-1,1),ylim=c(-1,1),
    main = paste("Between conditions, r=",round(cor(x[1,],x[3,]),2),sep = ""),
    xlab = "Pairwise map sim. vp1",ylab = "Pairwise map sim. vp2")
    axis(side=1,tck=-.05,pos=-1, at=c(-1,0,1))
    axis(side=2,tck=-.05,pos=-1, at=c(-1,0,1),las=1)
  }
  
  if(i==1){
    plot(x[1,],x[2,],cex=cexn,pch=20,axes=F,xlim=c(-.6,.6),ylim=c(-.6,.6),
    main = paste("Within condition, r=",round(cor(x[1,],x[2,]),2),sep = ""),
    xlab = "IFR association vp1",ylab = "IFR association vp1")
    axis(side=1,tck=-.05,pos=-.6, at=c(-.6,0,.6))
    axis(side=2,tck=-.05,pos=-.6, at=c(-.6,0,.6),las=1)
    
    x1=c(abs(x[1,]-x[2,]));  x2=c(abs(x[1,]-x[3,]))
    n=dim(x)[2]
    ci1=2*sd(x1)/sqrt(n)
    ci2=2*sd(x2)/sqrt(n)
    m1=mean(x1)
    m2=mean(x2)
    y=cbind(rep(.75,n),rep(1.88,n))
    
    barplot(c(m2,m1),ylim = c(0,.1),col = c("darkolivegreen2","grey47"),axes=F,ylab = "IFR association",main = "Reorganization IFR association")

    lines(c(y[1,2]-.125,y[1,2]+.125),c(m1+ci1,m1+ci1),lwd=1.0)
    lines(c(y[1,2]-.125,y[1,2]+.125),c(m1-ci1,m1-ci1),lwd=1.0)
    lines(c(y[1,1]-.125,y[1,1]+.125),c(m2+ci2,m2+ci2),lwd=1.0)
    lines(c(y[1,1]-.125,y[1,1]+.125),c(m2-ci2,m2-ci2),lwd=1.0)
    lines(c(y[1,2],y[1,2]),c(m1-ci1,m1+ci1),lwd=1.0)
    lines(c(y[1,1],y[1,1]),c(m2-ci2,m2+ci2),lwd=1.0)
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,0.1,0),at=c(y[1,1],y[1,2]),labels = c("Between","Within"))
    axis(side = 2, las=2, pos=0,tck=-0.05,at=ylab,mgp=c(0.15,0.7,0))
  }
  
  if(i==2){
    plot(x[1,],x[2,],cex=cexn,pch=20,axes=F,xlim=c(-1,1),ylim=c(-1,1),
    main = paste("Within condition, r=",round(cor(x[1,],x[2,]),2),sep = ""),
    xlab = "Pairwise map sim. vp1",ylab = "Pairwise map sim. vp1")
    axis(side=1,tck=-.05,pos=-1, at=c(-1,0,1))
    axis(side=2,tck=-.05,pos=-1, at=c(-1,0,1),las=1)
    
    x1=c(abs(x[1,]-x[2,]));  x2=c(abs(x[1,]-x[3,]))
    n=dim(x)[2]
    
    ci1=2*sd(x1)/sqrt(n)
    ci2=2*sd(x2)/sqrt(n)
    m1=mean(x1)
    m2=mean(x2)
    y=cbind(rep(.75,n),rep(1.88,n))

    barplot(c(m2,m1),ylim = c(0,.1),col = c("darkolivegreen2","grey47"),axes=F,ylab = "Pairwise map sim.",main = "Reorganization map sim.")

    lines(c(y[1,2]-.125,y[1,2]+.125),c(m1+ci1,m1+ci1),lwd=1.0)
    lines(c(y[1,2]-.125,y[1,2]+.125),c(m1-ci1,m1-ci1),lwd=1.0)
    lines(c(y[1,1]-.125,y[1,1]+.125),c(m2+ci2,m2+ci2),lwd=1.0)
    lines(c(y[1,1]-.125,y[1,1]+.125),c(m2-ci2,m2-ci2),lwd=1.0)
    lines(c(y[1,2],y[1,2]),c(m1-ci1,m1+ci1),lwd=1.0)
    lines(c(y[1,1],y[1,1]),c(m2-ci2,m2+ci2),lwd=1.0)
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,0.1,0),at=c(y[1,1],y[1,2]),labels = c("Between","Within"))
    axis(side = 2, las=2, pos=0,tck=-0.05,at=ylab,mgp=c(0.15,0.7,0))
  }
}
