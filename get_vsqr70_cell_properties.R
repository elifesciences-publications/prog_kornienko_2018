### function doing all the work on individual sessions
get_vsqr70_cell_properties<-function(rs){
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
  for (ll in 1:2){
    light=ll
    print(ll)
    id<-which(lightInt$V1==paste("l",light,sep=""))
    m<-matrix(data=as.numeric(c(lightInt[id,2], lightInt[id,3])),ncol=2)
    
    st<-setIntervals(st,s=m);
    ptsqr70<-setInvalidOutsideInterval(pt,s=m)
    ptsqr70<-speedFilter(ptsqr70,minSpeed=3,maxSpeed = 100)
    
    ## get the spatial properties during the sqr70 in order to identify the neurons
    sp<-getMapStats(sp,st,ptsqr70) ## get info score, sparsity, border, grid score
    # browser()
    # sp=autocorrelationDoughnut(sp)
    sp<-getMapStatsShuffle(sp,st,ptsqr70) ## get the shuffling values
    
    
    ## need the head direction selectivity for open field, histo, vector length and shuffled ##
    hd<-headDirectionStats(hd,st,ptsqr70)
    hd<-headDirectionStatsShuffle(hd,st,ptsqr70)
    
    ## get speed scores
    st<-ifr(st)
    sp<-speedScore(sp,st,ptsqr70,minSpeed=3,maxSpeed=100,runLm=F)
    sp<-speedScoreShuffle(sp,st,ptsqr70,minSpeed=3,maxSpeed=100)
    
    ## get mean firing rate of the neurons
    st<-meanFiringRate(st)
    
    ## create a data frame containing the data for each cell ##
    cells<-data.frame(mouse=rs@animalName,session=rs@session,cell.id=cg@id,tetrode.id=cg@tetrodeId,region=cg@brainRegion,
                      clu.to.tet=cg@cluToTetrode,
                      mean.rate=st@meanFiringRate, 
                      info.score=sp@infoScore, 
                      sparsity=sp@sparsity,
                      border.score=sp@borderScore, 
                      grid.score=sp@gridScore,
                      speed.score=sp@speedScore,
                      hd.vl=hd@vectorLength,hd.peak.rate=hd@peakRates)
    ## create a data frame with the shuffled data
    shuf<-data.frame(mouse=rs@animalName,session=rs@session,
                     cell.id=rep(cg@id,sp@nShufflings),
                     info.score=sp@infoScoreShuffle,
                     sparsity=sp@sparsityShuffle,
                     border.score=sp@borderScoreShuffle,
                     grid.score=sp@gridScoreShuffle, 
                     speed.score=sp@speedScoreShuffle,
                     hd.vl=hd@vectorLengthShuffle)
    
    if(ll==1){cells1=cells;shuf1=shuf}
    if(ll==2){cells2=cells;shuf2=shuf}
  }
  # return(list(cells1=cells1,shuf1=shuf1,cells2=cells2,shuf2=shuf2))
  return(list(cells1.test=cells1,shuf1.test=shuf1,cells2.test=cells2,shuf2.test=shuf2))
}
