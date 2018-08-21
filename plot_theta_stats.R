# plot
# distribution of theta indices
# HD scores and mean firing rates of non-rhythmic and theta-rhythmic HD cells
# average spike waveforms of non-rhythmic and theta-rhythmic HD cells

library(mclust)
t<-read.table(paste(ep@resultsDirectory,"cells.class.table",sep="/"),sep="\t")
load(paste(ep@resultsDirectory,"cells1",sep="/"))
load(paste(ep@resultsDirectory,"spectrum",sep="/"))
load(paste(ep@resultsDirectory,"spectrum.id",sep="/"))
load(paste(ep@resultsDirectory,"spectrum.freq",sep="/"))
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

x=theta.i[t$hd==1]
par(mfrow=c(2,3))
hist(x,main="Theta index distribution",ylab = "Number of cells", xlab="Theta index",15,xlim = c(-.05,.4),las=1)

x.gmm = Mclust(x)
x.s=summary(x.gmm)
print("Fit Gaussian finite mixture model")
print(paste("Number of components of best fit: ",x.s$G,sep=""))
print(paste("Log-likelhood: ",round(x.s$loglik,2),sep=""))
print(paste("BIC: ",round(x.s$bic,2),sep=""))
print("Theta index threshold = 0.07")
lines(c(0.07,0.07),c(0,14),col="red",lwd=2)
print(paste("Number of non-rhythmic (NR) HD cells (theta index threshold < 0.07): N = ",sum(x<.07),sep=""))
print(paste("Number of theta-rhythmic (TR) HD cells (theta index threshold > 0.07): N = ",sum(x>.07),sep=""))

##################################################################################
# HD score and average firing rate of non-rhythmic and theta-rhythmic HD cells

hdScoreNR=cells1$hd.vl[theta.i<.07 & t$hd==1]
hdScoreTR=cells1$hd.vl[theta.i>.07 & t$hd==1]
meanFrNR=cells1$mean.rate[theta.i<.07 & t$hd==1]
meanFrTR=cells1$mean.rate[theta.i>.07 & t$hd==1]

boxplot(hdScoreNR,hdScoreTR,col = c("red","gray"),las=1,
        ylim=c(0,1),names = c("NR","TR"),frame=F,ylab="HD score", main = "Directional tuning")
boxplot(meanFrNR,meanFrTR,col = c("red","gray"),las=1,
        ylim=c(0,15),names = c("NR","TR"),frame=F,ylab="Mean rate (Hz)",main="Firing rate")


print("HD selectivity of NR vs. TR HD cells during vp1 trials")
print(wilcox.test(hdScoreNR,hdScoreTR))
print("Mean firing rate of NR vs. TR HD cells during vp1 trials")
print(wilcox.test(meanFrNR,meanFrTR))

##################################################################################
# Spike waveform, duration and asymmetry of non-rhythmic and theta-rhythmic HD cells 

load(file=paste(ep@resultsDirectory,"/Waveform.stats.RData",sep=""))
ymin=-100;ymax=50
# calculate average waveform, correct voltage by factor 1/7
m1=apply(wf1s, 2, mean)/7
m2=apply(wf2s, 2, mean)/7


plot(time,m1,type="l",col="red",ylab = "Voltage (muV)",xlab="Time (msec)",
     ylim=c(ymin,ymax),axes = F,xlim=c(-.75,.75), main = "Spike waveform NR vs. TR HD cells")
lines(time,m2)
lines(c(time[which.max(m1[time>0])+100],time[which.max(m1[time>0])+100]),
      c(min(c(m1,m2)),max(m1[time>0])),col="red")
lines(c(time[which.max(m2[time>0])+100],time[which.max(m2[time>0])+100]),
      c(min(c(m1,m2)),max(m2[time>0])))
axis(1, labels=c(-.75,-.5,-.25,0,.25,.5,.75), at=c(-.75,-.5,-.25,0,.25,.5,.75))
axis(2, labels=c(ymin,-50,0,ymax), at=c(ymin,-50,0,ymax),las=1)

# Trough-to-peak duration
ymin=0.10;ymax=0.90
boxplot(duration1,duration2,col=c("red","gray"),axes=F,ylim=c(ymin,ymax), #pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),whisklty = 1,medlwd=1.5,outcex=.3
        ylab="Trough-to-peak duration (msec)")
axis(1, labels=c("NR","TR"), at=c(1,2))
axis(2, labels=c(ymin,.3,.5,0.7,ymax), at=c(ymin,.3,.5,.7,ymax),las=1)

# Peak amplitude asymmetry
boxplot(sasym1,sasym2,col=c("red","gray"),ylim=c(-1.05,1.05),axes=F,#pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),whisklty = 1,medlwd=1.5,outcex=.3,
        ylab="Peak amplitude asymmetry")
axis(1, labels=c("NR","TR"), at=c(1,2))
axis(2, labels=c(-1,0,1), at=c(-1,0,1),las=1)

durationNR=duration1
durationTR=duration2
peakAsymmetryNR=sasym1
peakAsymmetryTR=sasym2
print("Trough-to-peak duration non-rhythmic vs. theta-rhythmic cells")
print(wilcox.test(durationNR,durationTR))

print("Peak amplitude asymmetry non-rhythmic vs. theta-rhythmic cells")
print(wilcox.test(peakAsymmetryNR,peakAsymmetryTR))
