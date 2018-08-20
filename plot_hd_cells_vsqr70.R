t<-read.table(paste(ep@resultsDirectory,"cells.class.table",sep="/"),sep="\t")
t2<-read.table(paste(ep@resultsDirectory,"cells.sig.changes",sep="/"),sep="\t")

load(paste(ep@resultsDirectory,"hd.pairs.id",sep="/"))
load(paste(ep@resultsDirectory,"hd.pairs.hd",sep="/"))
load(paste(ep@resultsDirectory,"hd.pairs.vl",sep="/"))
load(paste(ep@resultsDirectory,"hd.pairs.ifra",sep="/"))

#########################################################################################
# plot 1) spike time autocorrelations, l1/l2 HD tuning curves, spatial
# firing rate maps of HD cells shown in Figure 4 and Figure 5 to pdf

cellid=c("jp865-04052016-0107_3","jp2159-16062016-0107_4", "jp2159-16062016-0107_6",
         "jp757-18082016-0107_8","jp759-28082016-0107_8","jp3302-20122016-0107_2",
         "jp759-18082016-0107_2", "jp759-28082016-0107_9","jp865-02062016-0107_6",
         "jp3302-21122016-0107_4")
source(paste(indir,'/plot_example_cells_hd.R',sep=""))

#########################################################################################
# plot distribution of theta indices, HD scores, firing rates and  
# spike waveforms of non-rhythmic and theta-rhythmic 
source(paste(indir,'/plot_theta_stats.R',sep=""))

#########################################################################################
# plot changes in HD score, preferred direction, firing rates and 
# BD scores between the two light conditions of non-rhythmic and 
# theta-rhythmic HD cells
source(paste(indir,'/plot_2light_changes_theta.R',sep=""))

# # percentage how many HD cells are changing between the two light conditions
# 
# x = (round(sum(t$hd==1 & (t2$hd.score==1 | t2$fr==1 | t2$hd.dir==1))/sum(t$hd==1)*100))
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in HD score, rate or direction: ",x,"%",sep=""))
# x = round(sum(t$hd==1 & t2$hd.score==1 & t2$fr!=1 & t2$hd.dir!=1)/sum(t$hd==1)*100)
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in HD score: ",x,"%",sep=""))
# x = round(sum(t$hd==1 & t2$hd.score==1 & t2$fr==1 & t2$hd.dir!=1)/sum(t$hd==1)*100)
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in HD score and rate: ",x,"%",sep=""))
# x = round(sum(t$hd==1 & t2$hd.score==1 & t2$fr!=1 & t2$hd.dir==1)/sum(t$hd==1)*100)
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in HD score and direction: ",x,"%",sep=""))
# x = round(sum(t$hd==1 & t2$hd.score==1 & t2$fr==1 & t2$hd.dir==1)/sum(t$hd==1)*100)
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in HD score, rate and direction: ",x,"%",sep=""))
# x = round(sum(t$hd==1 & t2$hd.score!=1 & t2$fr==1 & t2$hd.dir!=1)/sum(t$hd==1)*100)
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in rate: ",x,"%",sep=""))
# x = round(sum(t$hd==1 & t2$hd.score!=1 & t2$fr==1 & t2$hd.dir==1)/sum(t$hd==1)*100)
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in rate and direction: ",x,"%",sep=""))
# x = round(sum(t$hd==1 & t2$hd.score!=1 & t2$fr!=1 & t2$hd.dir==1)/sum(t$hd==1)*100)
# print(paste("Percentage of N = ",sum(t$hd==1), " hd cells with significant changes in direction: ",x,"%",sep=""))

#########################################################################################
# plot pairwise changes and reorganization scores in preferred HD, HD score and IFR
# association of simultaneously recorded HD cells

# non-rhythmic cell pairs
id1 = as.character(t$cell.id[t$hd==1 & theta.i<.07]) 
id1 = hd.pairs.id[1,] %in% id1 & hd.pairs.id[2,] %in% id1

# theta-rhythmic cell pairs
id2 = as.character(t$cell.id[t$hd==1 & theta.i>.07])  
id2 = hd.pairs.id[1,] %in% id2 & hd.pairs.id[2,] %in% id2  

# mixed non-rhythmic & theta-rhythmic cell pairs
id3 = id1==F & id2==F 
#########################################################################################
# CELL PAIRS EXAMPLES: 
# pair 1: jp2159-16062016−0107_2, jp2159-16062016-0107_4
# pair 2: jp759-25082016-0107_2, jp759-25082016-0107_4
# pair 3: jp759-28082016-0107_8, jp759-28082016-0107_9

#########################################################################################
cexn=0.8
par(mfrow=c(3,6))
for (i in 1:3){
  if(i==3){x=hd.pairs.ifra;ylab=c(0,.2,.4)}
  if(i==1){x=hd.pairs.hd;ylab=c(0,40,80)}
  if(i==2){x=hd.pairs.vl;ylab=c(0,.3,.6)}
  if(i==3){
    
    # CELL PAIR
    p1 = which(t$cell.id=="jp759-28082016-0107_8")
    p2 = which(t$cell.id=="jp759-28082016-0107_9")
    
    # SPIKE TIME AUTOCORRELATION
    xx1=st.auto[,p1]; xx2=st.auto[,p2]
    plot(xx1/sum(xx1),col ="black",axes = F,type = "l",cex.main=.4,xlab = "Time (ms)",ylab = "Spikes (norm.)")
    lines(xx2/sum(xx2),col="blue")
    axis(side = 1,  at=c(0,150,300), tck=-0.05,labels = c(-300,0,300))
    axis(side = 2, las=1,tck=-0.05)
    
    # VP1 HD TUNING CURVES OF CELL PAIR
    df1<-as.numeric(tuning.hd1[p1,])
    df2<-as.numeric(tuning.hd1[p2,])
    radlim=max(rbind(df1,df2))
    plotrix::polar.plot(rbind(df1,df2),#df1$rate,
                        polar.pos=seq(5,355,by=10),
                        labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                        clockwise=T,rp.type="p",
                        rad.col="white",
                        show.grid=T,show.radial.grid=T,
                        radial.lim=c(0,radlim),
                        # radial.lim=c(0,1.2),
                        show.grid.labels=0,
                        xlab="",ylab="",line.col=c("black","blue"),mar=c(1,0.75,1,0.75),cex.lab=0.3)
    mtext(paste(round(max(cells1$hd.peak.rate[p1])),"Hz ",round(cells1$hd.peak.rate[p2]),"Hz "),
          side=3,at=1,line=0.125,cex=.8)
    
    # VP2 HD TUNING CURVES OF CELL PAIR
    df1<-as.numeric(tuning.hd2[p1,])
    df2<-as.numeric(tuning.hd2[p2,])
    radlim=max(rbind(df1,df2))
    plotrix::polar.plot(rbind(df1,df2),#df1$rate,
                        polar.pos=seq(5,355,by=10),
                        labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                        clockwise=T,rp.type="p",
                        rad.col="white",
                        show.grid=T,show.radial.grid=T,
                        radial.lim=c(0,radlim),
                        # radial.lim=c(0,1.2),
                        show.grid.labels=0,
                        xlab="",ylab="",line.col=c("black","blue"),mar=c(1,0.75,1,0.75),cex.lab=0.3)
    mtext(paste(round(max(cells2$hd.peak.rate[p1])),"Hz ",round(cells2$hd.peak.rate[p2]),"Hz "),
          side=3,at=1,line=0.125,cex=.8)
    
    plot(x[1,],x[3,],cex=cexn,pch=20,axes=F,xlim=c(-.6,.6),ylim=c(-.6,.6),
    main = paste("Between conditions, r=",round(cor(x[1,],x[3,]),2),sep = ""),
    xlab = "IFR association vp1",ylab = "IFR association vp2")
    points(x[1,id1],x[3,id1],pch=20,cex=cexn,col="red")
    points(x[1,id3],x[3,id3],pch=20,cex=cexn,col="orange")
    axis(side=1,tck=-.05,pos=-.6, at=c(-.6,0,.6))
    axis(side=2,tck=-.05,pos=-.6, at=c(-.6,0,.6),las=1)
  }
  
  if(i==1){
    # CELL PAIR
    p1 = which(t$cell.id=="jp2159-16062016-0107_2")
    p2 = which(t$cell.id=="jp2159-16062016-0107_4")
    
    # SPIKE TIME AUTOCORRELATION
    xx1=st.auto[,p1]; xx2=st.auto[,p2]
    plot(xx1/sum(xx1),col ="black",axes = F,type = "l",cex.main=.4,xlab = "Time (ms)",ylab = "Spikes (norm.)")
    lines(xx2/sum(xx2),col="blue")
    axis(side = 1,  at=c(0,150,300), tck=-0.05,labels = c(-300,0,300))
    axis(side = 2, las=1,tck=-0.05)
    
    # VP1 HD TUNING CURVES OF CELL PAIR
    df1<-as.numeric(tuning.hd1[p1,])
    df2<-as.numeric(tuning.hd1[p2,])
    radlim=max(rbind(df1,df2))
    plotrix::polar.plot(rbind(df1/max(df1),df2/max(df2)),#df1$rate,
                        polar.pos=seq(5,355,by=10),
                        labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                        clockwise=T,rp.type="p",
                        rad.col="white",
                        show.grid=T,show.radial.grid=T,
                        # radial.lim=c(0,radlim),
                        radial.lim=c(0,1.2),
                        show.grid.labels=0,
                        xlab="",ylab="",line.col=c("black","blue"),mar=c(1,0.75,1,0.75),cex.lab=0.3)
    mtext(paste(round(max(cells1$hd.peak.rate[p1])),"Hz ",round(cells1$hd.peak.rate[p2]),"Hz "),
          side=3,at=1,line=0.125,cex=.8)
    
    # VP2 HD TUNING CURVES OF CELL PAIR
    df1<-as.numeric(tuning.hd2[p1,])
    df2<-as.numeric(tuning.hd2[p2,])
    radlim=max(rbind(df1,df2))
    plotrix::polar.plot(rbind(df1/max(df1),df2/max(df2)),#df1$rate,
                        polar.pos=seq(5,355,by=10),
                        labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                        clockwise=T,rp.type="p",
                        rad.col="white",
                        show.grid=T,show.radial.grid=T,
                        # radial.lim=c(0,radlim),
                        radial.lim=c(0,1.2),
                        show.grid.labels=0,
                        xlab="",ylab="",line.col=c("black","blue"),mar=c(1,0.75,1,0.75),cex.lab=0.3)
    mtext(paste(round(max(cells2$hd.peak.rate[p1])),"Hz ",round(cells2$hd.peak.rate[p2]),"Hz "),
          side=3,at=1,line=0.125,cex=.8)
    
  
    plot(x[1,],x[3,],cex=cexn,pch=20,axes=F,xlim=c(0,180),ylim=c(0,180),
    main = paste("Between conditions, r=",round(cor(x[1,],x[3,]),2),sep = ""),
    xlab = "Difference Pref. HD vp1 (deg) vp1",ylab = "Difference Pref. HD vp2")
    points(x[1,id1],x[3,id1],pch=20,cex=cexn,col="red")
    points(x[1,id3],x[3,id3],pch=20,cex=cexn,col="orange")
    axis(side=1,tck=-.05,pos=-1, at=c(0,90,180))
    axis(side=2,tck=-.05,pos=-1, at=c(0,90,180),las=1)
  }
  if(i==2){
    # CELL PAIR
    p1 = which(t$cell.id=="jp759-25082016-0107_2")
    p2 = which(t$cell.id=="jp759-25082016-0107_4")
    
    # SPIKE TIME AUTOCORRELATION
    xx1=st.auto[,p1]; xx2=st.auto[,p2]
    plot(xx1/sum(xx1),col ="black",axes = F,type = "l",cex.main=.4,xlab = "Time (ms)",ylab = "Spikes (norm.)")
    lines(xx2/sum(xx2),col="blue")
    axis(side = 1,  at=c(0,150,300), tck=-0.05,labels = c(-300,0,300))
    axis(side = 2, las=1,tck=-0.05)
    mtext("Cell1", at = c(260),col = c("black"),side=1,line = -1,cex=.8)
    mtext("Cell2", at = c(260),col = c("blue"),side=1,line = -2,cex=.8)

    
    # VP1 HD TUNING CURVES OF CELL PAIR
    df1<-as.numeric(tuning.hd1[p1,])
    df2<-as.numeric(tuning.hd1[p2,])
    radlim=max(rbind(df1,df2))
    plotrix::polar.plot(rbind(df1/max(df1),df2/max(df2)),#df1$rate,
                        polar.pos=seq(5,355,by=10),
                        labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                        clockwise=T,rp.type="p",
                        rad.col="white",
                        show.grid=T,show.radial.grid=T,
                        # radial.lim=c(0,radlim),
                        radial.lim=c(0,1.2),
                        show.grid.labels=0,
                        xlab="",ylab="",line.col=c("black","blue"),mar=c(1,0.75,1,0.75),cex.lab=0.3)
    mtext(paste(round(max(cells1$hd.peak.rate[p1])),"Hz ",round(cells1$hd.peak.rate[p2]),"Hz "),
          side=3,at=1,line=0.125,cex=.8)
    
    # VP1 HD TUNING CURVES OF CELL PAIR
    df1<-as.numeric(tuning.hd2[p1,])
    df2<-as.numeric(tuning.hd2[p2,])
    radlim=max(rbind(df1,df2))
    plotrix::polar.plot(rbind(df1/max(df1),df2/max(df2)),#df1$rate,
                        polar.pos=seq(5,355,by=10),
                        labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                        clockwise=T,rp.type="p",
                        rad.col="white",
                        show.grid=T,show.radial.grid=T,
                        # radial.lim=c(0,radlim),
                        radial.lim=c(0,1.2),
                        show.grid.labels=0,
                        xlab="",ylab="",line.col=c("black","blue"),mar=c(1,0.75,1,0.75),cex.lab=0.3)
    mtext(paste(round(max(cells2$hd.peak.rate[p1])),"Hz ",round(cells2$hd.peak.rate[p2]),"Hz "),
          side=3,at=1,line=0.125,cex=.8)
    
    plot(x[1,],x[3,],cex=cexn,pch=20,axes=F,xlim=c(-0.8,.8),ylim=c(-.8,.8),
         main = paste("Between conditions, r=",round(cor(x[1,],x[3,]),2),sep = ""),
         xlab = "HD score difference vp1",ylab = "HD score difference vp2")
    points(x[1,id1],x[3,id1],pch=20,cex=cexn,col="red")
    points(x[1,id3],x[3,id3],pch=20,cex=cexn,col="orange")
    axis(side=1,tck=-.05,pos=-.8, at=c(-0.8,0,.8))
    axis(side=2,tck=-.05,pos=-.8, at=c(-0.8,0,.8),las=1)
  }
  if(i==3){
    plot(x[1,],x[2,],cex=cexn,pch=20,axes=F,xlim=c(-.6,.6),ylim=c(-.6,.6),
    main = paste("Within condition, r=",round(cor(x[1,],x[2,]),2),sep = ""),
    xlab = "IFR association vp1",ylab = "IFR association vp1")
    points(x[1,id1],x[2,id1],pch=20,cex=cexn,col="red")
    points(x[1,id3],x[2,id3],pch=20,cex=cexn,col="orange")
    axis(side=1,tck=-.05,pos=-.6, at=c(-.6,0,.6))
    axis(side=2,tck=-.05,pos=-.6, at=c(-.6,0,.6),las=1)
    
    # NON-RHYTHMIC PAIRS confidence intervals for reorganization scores
    x1=c(abs(x[1,id1]-x[2,id1]));  x2=c(abs(x[1,id1]-x[3,id1]))
    n1=length(x1);ci1=2*sd(x1)/sqrt(n1);ci2=2*sd(x2)/sqrt(n1)
    m1=mean(x1);m2=mean(x2)
    # THETA-RHYTHMIC PAIRS
    x3=c(abs(x[1,id2]-x[2,id2]));  x4=c(abs(x[1,id2]-x[3,id2]))
    n1=length(x3);ci3=2*sd(x3)/sqrt(n1);ci4=2*sd(x4)/sqrt(n1)
    m3=mean(x3);m4=mean(x4)
    
    b=barplot(c(m2,m1,m4,m3),ylim = c(0,.4),col = c("darkolivegreen2","grey47","darkolivegreen2","grey47"),
            axes=F,ylab = "IFR association",main = "Reorganization IFR association")
    
    lines(c(b[2]-.125,b[2]+.125),c(m1+ci1,m1+ci1),lwd=1.0)
    lines(c(b[2]-.125,b[2]+.125),c(m1-ci1,m1-ci1),lwd=1.0)
    lines(c(b[1]-.125,b[1]+.125),c(m2+ci2,m2+ci2),lwd=1.0)
    lines(c(b[1]-.125,b[1]+.125),c(m2-ci2,m2-ci2),lwd=1.0)
    lines(c(b[2],b[2]),c(m1-ci1,m1+ci1),lwd=1.0)
    lines(c(b[1],b[1]),c(m2-ci2,m2+ci2),lwd=1.0)
    
    lines(c(b[4]-.125,b[4]+.125),c(m3+ci3,m3+ci3),lwd=1.0)
    lines(c(b[4]-.125,b[4]+.125),c(m3-ci3,m3-ci3),lwd=1.0)
    lines(c(b[3]-.125,b[3]+.125),c(m4+ci4,m4+ci4),lwd=1.0)
    lines(c(b[3]-.125,b[3]+.125),c(m4-ci4,m4-ci4),lwd=1.0)
    lines(c(b[4],b[4]),c(m3-ci3,m3+ci3),lwd=1.0)
    lines(c(b[3],b[3]),c(m4-ci4,m4+ci4),lwd=1.0)
    
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,0.1,0),at=b,labels = c("Between","Within","Between","Within"))
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,2,0),at=c((b[2]+b[1])/2,(b[4]+b[3])/2),labels = c("NR","TR"))
    axis(side = 2, las=2, pos=0,tck=-0.05,at=ylab,mgp=c(0.15,0.7,0))
    
  }
  
  if(i==1){
    plot(x[1,],x[2,],cex=cexn,pch=20,axes=F,xlim=c(0,180),ylim=c(0,180),
    main = paste("Within condition, r=",round(cor(x[1,],x[2,]),2),sep = ""),
    xlab = "Difference Pref. HD vp1",ylab = "Difference Pref. HD vp1")
    points(x[1,id1],x[2,id1],pch=20,cex=cexn,col="red")
    points(x[1,id3],x[2,id3],pch=20,cex=cexn,col="orange")
    axis(side=1,tck=-.05,pos=-1, at=c(0,90,180))
    axis(side=2,tck=-.05,pos=-1, at=c(0,90,180),las=1)
    # NON-RHYTHMIC PAIRS confidence intervals for reorganization scores
    x1=c(abs(x[1,id1]-x[2,id1]));  x2=c(abs(x[1,id1]-x[3,id1]))
    n1=length(x1);ci1=2*sd(x1)/sqrt(n1);ci2=2*sd(x2)/sqrt(n1)
    m1=mean(x1);m2=mean(x2)
    # THETA-RHYTHMIC PAIRS
    x3=c(abs(x[1,id2]-x[2,id2]));  x4=c(abs(x[1,id2]-x[3,id2]))
    n1=length(x3);ci3=2*sd(x3)/sqrt(n1);ci4=2*sd(x4)/sqrt(n1)
    m3=mean(x3);m4=mean(x4)

    b=barplot(c(m2,m1,m4,m3),ylim = c(0,80),col = c("darkolivegreen2","grey47","darkolivegreen2","grey47"),
              axes=F,ylab = "Pref. HD (deg)",main = "Reorganization pref. HD ")

    lines(c(b[2]-.125,b[2]+.125),c(m1+ci1,m1+ci1),lwd=1.0)
    lines(c(b[2]-.125,b[2]+.125),c(m1-ci1,m1-ci1),lwd=1.0)
    lines(c(b[1]-.125,b[1]+.125),c(m2+ci2,m2+ci2),lwd=1.0)
    lines(c(b[1]-.125,b[1]+.125),c(m2-ci2,m2-ci2),lwd=1.0)
    lines(c(b[2],b[2]),c(m1-ci1,m1+ci1),lwd=1.0)
    lines(c(b[1],b[1]),c(m2-ci2,m2+ci2),lwd=1.0)
    
    lines(c(b[4]-.125,b[4]+.125),c(m3+ci3,m3+ci3),lwd=1.0)
    lines(c(b[4]-.125,b[4]+.125),c(m3-ci3,m3-ci3),lwd=1.0)
    lines(c(b[3]-.125,b[3]+.125),c(m4+ci4,m4+ci4),lwd=1.0)
    lines(c(b[3]-.125,b[3]+.125),c(m4-ci4,m4-ci4),lwd=1.0)
    lines(c(b[4],b[4]),c(m3-ci3,m3+ci3),lwd=1.0)
    lines(c(b[3],b[3]),c(m4-ci4,m4+ci4),lwd=1.0)
    
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,0.1,0),at=b,labels = c("Between","Within","Between","Within"))
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,2,0),at=c((b[2]+b[1])/2,(b[4]+b[3])/2),labels = c("NR","TR"))
    axis(side = 2, las=2, pos=0,tck=-0.05,at=ylab,mgp=c(0.15,0.7,0))
  }
  if(i==2){
    plot(x[1,],x[2,],cex=cexn,pch=20,axes=F,xlim=c(-.8,.8),ylim=c(-.8,.8),
         main = paste("Within condition, r=",round(cor(x[1,],x[2,]),2),sep = ""),
         xlab = "HD score difference vp1",ylab = "HD score difference vp1")
    points(x[1,id1],x[2,id1],pch=20,cex=cexn,col="red")
    points(x[1,id3],x[2,id3],pch=20,cex=cexn,col="orange")
    axis(side=1,tck=-.05,pos=-.8, at=c(-.8,0,.8))
    axis(side=2,tck=-.05,pos=-.8, at=c(-.8,0,.8),las=1)
    mtext("NR", at = c(-.5),col = c("red"),side=3,line = -1,cex=.8)
    mtext("mixed", at = c(-.5),col = c("orange"),side=3,line = -2,cex=.8)
    mtext("TR", at = c(-.5),col = c("black"),side=3,line = -3,cex=.8)
    
    x1=c(abs(x[1,id1]-x[2,id1]));  x2=c(abs(x[1,id1]-x[3,id1]))
    n1=length(x1);ci1=2*sd(x1)/sqrt(n1);ci2=2*sd(x2)/sqrt(n1)
    m1=mean(x1);m2=mean(x2)
    
    x3=c(abs(x[1,id2]-x[2,id2]));  x4=c(abs(x[1,id2]-x[3,id2]))
    n1=length(x3);ci3=2*sd(x3)/sqrt(n1);ci4=2*sd(x4)/sqrt(n1)
    m3=mean(x3);m4=mean(x4)
    
    b=barplot(c(m2,m1,m4,m3),ylim = c(0,.6),col = c("darkolivegreen2","grey47","darkolivegreen2","grey47"),
              axes=F,ylab = "IFR association",main = "Reorganization IFR association")
    
    lines(c(b[2]-.125,b[2]+.125),c(m1+ci1,m1+ci1),lwd=1.0)
    lines(c(b[2]-.125,b[2]+.125),c(m1-ci1,m1-ci1),lwd=1.0)
    lines(c(b[1]-.125,b[1]+.125),c(m2+ci2,m2+ci2),lwd=1.0)
    lines(c(b[1]-.125,b[1]+.125),c(m2-ci2,m2-ci2),lwd=1.0)
    lines(c(b[2],b[2]),c(m1-ci1,m1+ci1),lwd=1.0)
    lines(c(b[1],b[1]),c(m2-ci2,m2+ci2),lwd=1.0)
    
    lines(c(b[4]-.125,b[4]+.125),c(m3+ci3,m3+ci3),lwd=1.0)
    lines(c(b[4]-.125,b[4]+.125),c(m3-ci3,m3-ci3),lwd=1.0)
    lines(c(b[3]-.125,b[3]+.125),c(m4+ci4,m4+ci4),lwd=1.0)
    lines(c(b[3]-.125,b[3]+.125),c(m4-ci4,m4-ci4),lwd=1.0)
    lines(c(b[4],b[4]),c(m3-ci3,m3+ci3),lwd=1.0)
    lines(c(b[3],b[3]),c(m4-ci4,m4+ci4),lwd=1.0)
    
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,0.1,0),at=b,labels = c("Between","Within","Between","Within"))
    axis(side = 1, pos=0,tck=-0.05,mgp=c(0.5,2,0),at=c((b[2]+b[1])/2,(b[4]+b[3])/2),labels = c("NR","TR"))
    axis(side = 2, las=2, pos=0,tck=-0.05,at=ylab,mgp=c(0.15,0.7,0))
  }
}
