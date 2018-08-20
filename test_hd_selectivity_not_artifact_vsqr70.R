test_hd_selectivity_not_artifact_vsqr70<-function(rs){
  rs<-loadRecSession(rs)
  print(rs@session)
  myList<-getRecSessionObjects(rs)
  pt<-myList$pt
  st<-myList$st
  sp<-myList$sp
  hd<-myList$hd
  cg<-myList$cg
  ptsqr70<-speedFilter(pt,minSpeed=3,maxSpeed = 100)
  #################################################################################################
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
   DR1=c();DR2=c();DR.id=c();
    for (condition in 1:2) {
      l=paste("l",condition,sep="")
      m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition==l),2], 
                                lightInt[which(lightInt$condition==l),3])),ncol=2)
      hd<-headDirectionHisto(hd,st,ptsqr70);hd<-headDirectionStats(hd,st,ptsqr70)
      st<-setIntervals(st,s=m)
      # set parameters for directional distributive ratio
      sp@cmPerBin=5
      sp@smoothOccupancySd=0
      sp@smoothRateMapSd=0
      hd@smoothOccupancySd=0
      hd@smoothRateHistoSd=0
      hd@degPerBin = 10
      
      hd<-headDirectionHisto(hd,st,ptsqr70) # observed histo
      sp<-firingRateMap2d(sp,st,pt,nRowMap=45,nColMap=45)
      dr<-directionalDistributiveRatioFromHdHisto(sp,st,ptsqr70,hd,nRowMap=45,nColMap=45) 
      
      if(condition==1){DR1=c(DR1,dr);DR.id=c(DR.id,cg@id)}
      if(condition==2){DR2=c(DR2,dr)}
    }
  ##############################################################################################
   return(list(DR1=t(DR1),DR2=t(DR2),DR.id=t(DR.id)))
}
