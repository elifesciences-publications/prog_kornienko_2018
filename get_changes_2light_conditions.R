### function doing all the work on individual sessions
get_changes_2light_conditions<-function(rs){
  print(paste(rs@session,rs@path))
  
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
  colnames(lightInt)<-c("condition","s","e")
  #################################################################################################
  ## remove trials
  file=paste(rs@path,"/",rs@session,".light_intervals.trial_id",sep="")
  if(!file.exists(file))
    stop(paste("sqr70.cell.properties:",rs@session, "needs",file))
  removeTrl<-read.csv(file,header=T,sep="\t")
  lightInt<-lightInt[which(!removeTrl[,2]),]
  #################################################################################################
  
  xbins=40; ybins=40
  
  ## create empty data.frame to store the data
  light.changes.shuffle<-data.frame()
  
  ## select data when the animal ran between 3 and 100
  ptsqr70<-speedFilter(pt,minSpeed=3,maxSpeed = 100)
  
  
  ######## real conditions ##########
  ## get intervals for l1
  m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l1"),2], 
                              lightInt[which(lightInt$condition=="l1"),3])),ncol=2)
  st<-setIntervals(st,s=m)
  
  sp1<-getMapStats(sp,st,ptsqr70)
  sp1<-firingRateMap2d(sp1,st,ptsqr70,nRowMap = xbins,nColMap = ybins)
  sp1<-mapSpatialAutocorrelation(sp1)
  maps2d1<-sp1@maps
  autos2d1<-sp1@autos
  r1<-meanFiringRate(st)@meanFiringRate

  
  hd0<-headDirectionHisto(hd,st,ptsqr70)
  hdSt1<-headDirectionStats(hd0,st,ptsqr70)
  tuning.hd1<-data.frame(firingRate=t(hd0@histo))
  dir1<-hdSt1@meanDirection
  v1<-hdSt1@vectorLength

  
  ## get intervals for l2
  m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l2"),2], 
                              lightInt[which(lightInt$condition=="l2"),3])),ncol=2)
  st<-setIntervals(st,s=m)
  
  sp2<-getMapStats(sp,st,ptsqr70)
  sp2<-firingRateMap2d(sp2,st,ptsqr70,nRowMap = xbins,nColMap = ybins)
  sp2<-mapSpatialAutocorrelation(sp2)
  maps2d2<-sp2@maps
  autos2d2<-sp2@autos
  r2<-meanFiringRate(st)@meanFiringRate
  
  hd0<-headDirectionHisto(hd,st,ptsqr70)
  hdSt2<-headDirectionStats(hd0,st,ptsqr70)
  tuning.hd2<-data.frame(firingRate=t(hd0@histo))
  dir2<-hdSt2@meanDirection
  v2<-hdSt2@vectorLength

  
  light.changes<-data.frame(cell.id=cg@id, 
                          map.cor=firingRateMapCorrelation(sp1,sp2),
                          mean.rate.diff=(r1-r2)/(r1+r2),
                          mean.dir=180-abs(abs(dir1 - dir2) - 180),
                          vectorLength.diff=v1-v2)


  
  
  ##### shuffling procedure #####
  n=500
  for(i in 1:n)
  {
    ## reassign l1 and l2 randomly
    lightInt$condition<-sample(lightInt$condition) 
    ## get intervals for l1
    m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l1"),2], 
                                lightInt[which(lightInt$condition=="l1"),3])),ncol=2)
    st<-setIntervals(st,s=m)
    # obtain all scores for l1
    r1<-meanFiringRate(st)@meanFiringRate
    sp1<-getMapStats(sp,st,ptsqr70)
    hd0<-headDirectionHisto(hd,st,ptsqr70)
    hdSt1<-headDirectionStats(hd0,st,ptsqr70)
    dir1<-hdSt1@meanDirection
    v1<-hdSt1@vectorLength

    
    ## get intervals for l2
    m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l2"),2], 
                                lightInt[which(lightInt$condition=="l2"),3])),ncol=2)
    st<-setIntervals(st,s=m)
    ## obtain all scores for l2
    r2<-meanFiringRate(st)@meanFiringRate
    sp2<-getMapStats(sp,st,ptsqr70)
    hd0<-headDirectionHisto(hd,st,ptsqr70)
    hdSt2<-headDirectionStats(hd0,st,ptsqr70)
    dir2<-hdSt2@meanDirection
    v2<-hdSt2@vectorLength

    
    ###########################################################################
    # calculate differences between scores from light 1 and light 2
    light.changes.shuffle<-rbind(light.changes.shuffle,
                               data.frame(cell.id=cg@id, 
                                          shuffle=i,
                                          map.cor=firingRateMapCorrelation(sp1,sp2),
                                          mean.rate.diff=(r1-r2)/(r1+r2),
                                          mean.dir=180-abs(abs(dir1 - dir2) - 180),
                                          vectorLength.diff=v1-v2))
  }
  
  return(list(light.changes=light.changes, light.changes.shuffle=light.changes.shuffle,
              maps2d.l1=maps2d1,maps2d.l2=maps2d2,autos2d.l1=autos2d1,autos2d.l2=autos2d2,
              tuning.hd1=tuning.hd1,tuning.hd2=tuning.hd2))
}
