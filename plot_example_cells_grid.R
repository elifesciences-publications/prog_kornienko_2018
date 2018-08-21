plot.hist.one.distribution.scale<-function(x, axis.y.las=2,
                                           outma=c(2,1,1.5,0),margin=c(2,2,1,1),
                                           mgp.x=c(0.5,0.1,0),mgp.y=c(0.8,0.2,0),
                                           probability=1, ylab="",xlab="",vline="",...)
{
  if(vline==""){y=x}else{y=c(x,vline)}
  
  m<-min(y);M<-max(y);r<-M-m;r20<-r*0.2;m<-m-r20
  M<-M+r20;int=(M-m)/40;plotxlim=c(m,M)
  h<-hist(x,breaks=seq(m,M,int),plot=F);plotylim=c(0,max(h$counts))

  par(oma=outma,mar=margin,cex.lab=0.6,cex.axis=0.6)
  hist(x,xlim=plotxlim,breaks=seq(m,M,int),axes=F,labels=F,main="",xlab="",ylab="",...)
  if(vline!="")
  {
    lines(c(vline,vline),c(plotylim[1],plotylim[2]),col="red",lty=1,lwd=1)
  }
  axis(side = 1, pos=plotylim[1],tck=-0.05,cex.axis=0.60,mgp=c(0.5,0.1,0))
  axis(side = 2, las=axis.y.las, pos=plotxlim[1],tck=-0.05,cex.axis=0.60,mgp=c(0.8,0.2,0))
  title(xlab=xlab,mgp=mgp.x)
  title(ylab=ylab,mgp=mgp.y)
}



########################################################################################################

plot.page.full.maps<-function(fn="map.changes.stats.pdf",
                              cell.list=cellid, 
                              cells1=cells1,
                              cells2=cells2,
                              tuning.hd1=tuning.hd1,
                              tuning.hd2=tuning.hd2,
                              maps2d.l1=maps2d.l1,
                              maps2d.l2=maps2d.l2,
                              st.auto=st.auto,
                              lc=light.changes,
                              lcs=light.changes.shuffle){
  
  jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  
  num.cols<- 7
  num.rows<-8
  plot.per.page<-num.cols*num.rows
  m<-matrix(c(rep(seq(0,1-(1/num.cols),1/num.cols),num.rows),
              rep(seq(1/num.cols,1,1/num.cols),num.rows),
              rep(seq(1-(1/num.rows),0,0-1/num.rows),each=num.cols),
              rep(seq(1,1/num.rows,0-1/num.rows),each=num.cols)),ncol=4)
  
  pdf(file=fn,onefile=TRUE,paper="a4",height=10, width = 7.5)
  index=1; 
  for ( i in 1:length(cell.list)) ## any cell list
  {
    if(index==1){split.screen(m)}
    #print(i)
    
    id1=which(light.changes$cell.id==cellid[i])
    id2=which(light.changes.shuffle$cell.id==cellid[i])

    mgp.x=c(0.5,0.0,0);mgp.y=c(1,0.2,0)
    
    ################################ spike-time autocorrelations
    screen(index)
    # par(oma = c(2,1,1.5,0), mar = c(2,2,1,1),cex.lab=0.6,cex.axis=0.6)
    # xx=st.auto[,id1]
    # plot(xx,col ="black",axes = F,type = "l",main = cellid[i],cex.main=.4)
    # axis(side = 1,  at=c(0,150,300), tck=-0.05,cex.axis=0.60,labels = c(-300,0,300),mgp=mgp.x)
    # axis(side = 2, las=1,tck=-0.05,cex.axis=0.60,mgp=mgp.y)
    # title(xlab="Time (ms)",mgp=mgp.x)
    # title(ylab="Spikes",mgp=mgp.y)
    index=index+1

    ################################ HD polar plots l1, l2
    df1<-as.numeric(tuning.hd1[id1,])
    df2<-as.numeric(tuning.hd2[id1,])
    radlim=max(rbind(df1,df2))
    
    screen(index)
   
    par(cex.lab=0.5,cex.axis=0.5)
    plotrix::polar.plot(rbind(df1/max(df1),df2/max(df2)),#df1$rate,
                        polar.pos=seq(5,355,by=10),
                        labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                        clockwise=T,rp.type="p",
                        rad.col="white",
                        show.grid=T,show.radial.grid=T,
                        #radial.lim=c(0,radlim),
                        radial.lim=c(0,1.2),
                        show.grid.labels=0,
                        xlab="",ylab="",line.col=c("black","cyan4"),mar=c(1,0.75,1,0.75),cex.lab=0.3)
    mtext(paste(round(max(cells1$hd.peak.rate[id1])),"Hz ",round(cells2$hd.peak.rate[id1]),"Hz "),
          side=3,at=0,line=0.125,cex=0.4)
    index=index+1
    
    ################################ spatial firing rate maps l1, l2
    x=maps2d.l1[seq(dim(maps2d.l1)[1],1,-1),,id1]
    y=maps2d.l2[seq(dim(maps2d.l1)[1],1,-1),,id1]
    
    screen(index)  
    firingRateMapPlot(x)
    index=index+1
    
    screen(index)
    firingRateMapPlot(y)
    index=index+1
  
   
    ################################ differences in preferred HD, HD score and firing rate
    
    xh<-lc[id1,]
    xhs<-try(lcs[id2,],TRUE)
    
     screen(index)
     plot.hist.one.distribution.scale(x=xhs$map.cor,
         xlab=expression(paste("Map sim.")),ylab="Shuffles",vline=xh$map.cor)
     index=index+1
     
     screen(index)
     plot.hist.one.distribution.scale(x=xhs$mean.dir,
                                      xlab=expression(paste(Delta," Direction")),ylab="Shuffles",vline=xh$mean.dir)
     index=index+1
     
     screen(index)
     plot.hist.one.distribution.scale(x=xhs$mean.rate.diff,
                                      xlab=expression(paste(Delta," Rate (Hz)")),ylab="Shuffles",vline=xh$mean.rate.diff)
     
     if(index==plot.per.page){close.screen(all = TRUE);index=0}
     index=index+1
  }
  close.screen(all = TRUE)
  dev.off()
  

}

######################################################################################


load(paste(ep@resultsDirectory,"light.changes",sep="/"))
load(paste(ep@resultsDirectory,"light.changes.shuffle",sep="/"))
load(paste(ep@resultsDirectory,"maps2d.l1",sep="/"))
load(paste(ep@resultsDirectory,"maps2d.l2",sep="/"))
load(paste(ep@resultsDirectory,"tuning.hd1",sep="/"))
load(paste(ep@resultsDirectory,"tuning.hd2",sep="/"))
load(paste(ep@resultsDirectory,"cells1",sep="/"))
load(paste(ep@resultsDirectory,"cells2",sep="/"))
load(paste(ep@resultsDirectory,"st.auto",sep="/"))


print(paste("Plot polar plots, spatial firing rate maps and changes in pref. HD, HD score and rate in  ",ep@resultsDirectory,"/plot_example_cells_grid.pdf",sep=""))
plot.page.full.maps(fn=paste(ep@resultsDirectory,"plot_example_cells_grid.pdf",sep="/"),
                    cell.list=cellid, 
                    cells1=cells1,
                    cells2=cells2,
                    tuning.hd1=tuning.hd1,
                    tuning.hd2=tuning.hd2,
                    maps2d.l1=maps2d.l1,
                    maps2d.l2=maps2d.l2,
                    st.auto=st.auto,
                    lc=light.changes,
                    lcs=light.changes.shuffle)



