t<-read.table(paste(ep@resultsDirectory,"cells.class.table",sep="/"),sep="\t")
load(paste(ep@resultsDirectory,"cells1",sep="/"))
load(paste(ep@resultsDirectory,"spectrum",sep="/"))
load(paste(ep@resultsDirectory,"spectrum.id",sep="/"))
load(paste(ep@resultsDirectory,"spectrum.freq",sep="/"))
load(paste(ep@resultsDirectory,"light.changes",sep="/"))
load(paste(ep@resultsDirectory,"light.changes.shuffle",sep="/"))
load(paste(ep@resultsDirectory,"maps2d.l1",sep="/"))
load(paste(ep@resultsDirectory,"maps2d.l2",sep="/"))
load(paste(ep@resultsDirectory,"bdscore.hd1",sep="/"))
load(paste(ep@resultsDirectory,"bdscore.hd2",sep="/"))
##################################################################################
# calculate theta index from power spectra of instantaneous firing rates
freq=spectrum.freq[1,]
theta.i=c()
for (i in 1:dim(spectrum)[2]){
  wf=spectrum[,i]
  th=mean(wf[freq>6 & freq<10])
  b=mean(c(wf[freq>3 & freq<5],wf[freq>11 & freq<13]))
  ti=(th-b)/(th+b)
  theta.i=c(theta.i,ti)
}

##################################################################################
# 1) plot changes in preferred direction of non-rhythmic and theta-rhythmic HD cells 

par(mfrow=c(2,4))
x=theta.i[t$hd==1]
y=light.changes$mean.dir[t$hd==1]
plot(x[x<.07],y[x<.07],col="red",pch=19,frame=F,ylim=c(0,200),xlim=c(-.05,.4),ylab = "Pref. HD difference (deg)",
     xlab = "Theta index",las=1,main="Pref. HD difference (deg)")
points(x[x>.07],y[x>.07],col="gray34",pch=19) 

ymin=0;ymax=200
boxplot(y[x<.07],y[x>.07],col=c("red","gray"),axes=F,ylim=c(ymin,ymax), #pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),whisklty = 1,medlwd=1.5,outcex=.3
        ylab="Pref. HD difference (deg)")
axis(1, labels=c("NR","TR"), at=c(1,2))
axis(2, labels=c(ymin,50,100,150,ymax), at=c(ymin,50,100,150,ymax),las=1)

##################################################################################
# 2) changes in HD score of non-rhythmic and theta-rhythmic HD cells 

y=abs(light.changes$vectorLength.diff [t$hd==1])
plot(x[x<.07],y[x<.07],col="red",pch=19,frame=F,ylim=c(0,.5),xlim=c(-.05,.4),ylab = "HD score difference",
     xlab = "Theta index",las=1,main="HD score difference")
points(x[x>.07],y[x>.07],col="gray34",pch=19) 

ymin=0;ymax=.5
boxplot(y[x<.07],y[x>.07],col=c("red","gray"),axes=F,ylim=c(ymin,ymax), #pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),whisklty = 1,medlwd=1.5,outcex=.3
        ylab="HD score difference")
axis(1, labels=c("NR","TR"), at=c(1,2))
axis(2, labels=c(ymin,.25,ymax), at=c(ymin,.25,ymax),las=1)

##################################################################################
# 3) changes in firing rate of non-rhythmic and theta-rhythmic HD cells 

y=abs(light.changes$mean.rate.diff [t$hd==1])
plot(x[x<.07],y[x<.07],col="red",pch=19,frame=F,ylim=c(0,1),xlim=c(-.05,.4),
     ylab = "Rel. rate change",xlab = "Theta index",las=1,main="Rel. rate change")
points(x[x>.07],y[x>.07],col="gray34",pch=19) 

ymin=0;ymax=1
boxplot(y[x<.07],y[x>.07],col=c("red","gray"),axes=F,ylim=c(ymin,ymax), 
        ylab="Rel. rate change")
axis(1, labels=c("NR","TR"), at=c(1,2))
axis(2, labels=c(ymin,.5,ymax), at=c(ymin,.5,ymax),las=1)


##################################################################################
# 4) changes in bidirectionality of non-rhythmic and theta-rhythmic HD cells 

bdh1=bdscore.hd1;bdh2=bdscore.hd2
y = abs(bdh1$bdscore[t$hd==1]-bdh2$bdscore[t$hd==1])
plot(x[x<.07],y[x<.07],col="red",pch=19,frame=F,ylim=c(0,.7),xlim=c(-.05,.4),ylab = "BD score difference",
     xlab = "Theta index",las=1,main="BD score difference")
points(x[x>.07],y[x>.07],col="gray34",pch=19) 
ymin=0;ymax=.7
boxplot(y[x<.07],y[x>.07],col=c("red","gray"),axes=F,ylim=c(ymin,ymax), 
        ylab="BD score difference")
axis(1, labels=c("NR","TR"), at=c(1,2))
axis(2, labels=c(ymin,.35,ymax), at=c(ymin,.35,ymax),las=1)
