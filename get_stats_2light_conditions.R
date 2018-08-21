### function doing all the work on individual sessions
sqr70_spatial_properties_2_conditions<-function(rs){
  print(rs@session)
 
  ## load the data in R
  myList<-getRecSessionObjects(rs)
  st<-myList$st
  pt<-myList$pt
  cg<-myList$cg
  sp<-myList$sp
  hd<-myList$hd
  
  ## get the position data for sqr70, but only with l1 on
  file=paste(rs@path,"/",rs@session,".light_intervals",sep="")
  if(!file.exists(file))
    stop(paste("sqr70.cell.properties:",rs@session, "needs",file))
  lightInt<-read.table(file,header=F)
  #################################################################################################
  ## remove trials
  file=paste(rs@path,"/",rs@session,".light_intervals.trial_id",sep="")
  if(!file.exists(file))
    stop(paste("sqr70.cell.properties:",rs@session, "needs",file))
  removeTrl<-read.csv(file,header=T,sep="\t")
  lightInt<-lightInt[which(!removeTrl[,2]),]
  #################################################################################################
  ## number of bins in firing rate maps
  xbins=40;  ybins=40
  maps2d<-array(dim=c(xbins,ybins,0))
  autos2d<-array(dim=c(xbins*2+1,ybins*2+1,0))
  stats2d<-data.frame()
  tuning.hd<-data.frame()
  
  for(condition in c("l1","l2")){
    m<-matrix(data=as.numeric(c(lightInt[which(lightInt$V1==condition),2], 
                                lightInt[which(lightInt$V1==condition),3])),ncol=2)
    ptsqr70<-speedFilter(pt,minSpeed=3,maxSpeed = 100)
    st<-setIntervals(st,s=m)
    ## get the spatial properties during the sqr70 in order to identify the neurons
    sp<-getMapStats(sp,st,ptsqr70) ## get info score, sparsity, border, grid score
    df<-statsAsDataFrame(sp,shuffle=FALSE)    
    ## get spatial firing rate maps for each trial type (l1,l2) and their autocorrelations  
    sp<-firingRateMap2d(sp,st,ptsqr70,nRowMap = xbins,nColMap = ybins)
    sp<-mapSpatialAutocorrelation(sp)
    maps2d<-abind::abind(maps2d,sp@maps)
    autos2d<-abind::abind(autos2d,sp@autos)
    ## get firing rate
    st<-meanFiringRate(st)
    ## directional tuning
    hd<-headDirectionHisto(hd,st,ptsqr70)
    hdSt<-headDirectionStats(hd,st,ptsqr70)
    tuning.hd<-rbind(tuning.hd,data.frame(firingRate=t(hd@histo)))
    
    ## paste together the information about each cell
    ## stats2d data.frame can be used as index lookup for maps2d
    stats2d<-rbind(stats2d,data.frame(clu.id=df$clu.id,
                                      tet.id=cg@tetrode,
                                      condition=condition,
                                      peakRate=df$peakRate,
                                      infoScore=df$infoScore,
                                      sparsity=df$sparsity,
                                      gridScore=df$gridScore,
                                      meanRate=st@meanFiringRate,
                                      hd.peakRates=hdSt@peakRates,
                                      hd.meanDirection=hdSt@meanDirection,
                                      hd.vectorLength=hdSt@vectorLength))
  }
  return(list(maps2d=maps2d,stats2d=stats2d,autos2d=autos2d))
}
