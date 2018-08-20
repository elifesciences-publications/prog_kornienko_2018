get_bdscores_vsqr70<-function(rs){
  ################################################
  # calculate local maxima
  localMaxima <- function(x) {
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
      y <- y[-1]
    }
    y
  }
  ################################################
  # calculate birectional score of HD tuning curves
  get_bdscore <- function(thd){
    n = length(thd)
    #smooth head direction tuning curve and identify local maxima
    idp=localMaxima(smooth.spline(thd,spar = 0.3)$y)
    if(length(idp)==0){idp=which.max(smooth.spline(thd,spar = 0.3)$y)}
    if(length(idp)> 2){idp=localMaxima(smooth.spline(thd,spar = 0.4)$y)}
    lmax=idp
    
    # if local maxima at both edges of tuning curve take the larger one
    if(1%in%lmax & n%in%lmax){lmax=lmax[lmax!=n]}
    if(1%in%lmax & (n-1)%in%lmax){lmax=lmax[lmax!=(n-1)]}
    if(2%in%lmax & n%in%lmax){lmax=lmax[lmax!=n]}
    if(n%in%lmax & thd[1]>thd[n]){lmax=lmax[lmax!=n]}
    
    # sort local maxima in decreasing order
    idx=order(thd[lmax],decreasing=T)
    # firing rates of local maxima
    idxmax=thd[lmax[idx]]
    # number of local maxima, can be one or two
    numpeaks=((length(idxmax)>1)+0)*2 + ((length(idxmax)==1)+0)
    # if more than one take the first two highest values
    if(numpeaks>1){
      lmax=lmax[idx[1:2]]
      bdscore=thd[lmax[2]]/thd[lmax[1]]
      idxdeg=seq(5,355,by=10)[lmax]
      degdiff=180-abs(abs(idxdeg[1] - idxdeg[2]) - 180)
      #get minimum between maxima
      lmax1=min(lmax);lmax2=max(lmax);
      if((lmax2-lmax1)<=(n/2)){
      lmin=  localMaxima(-thd[lmax1:lmax2])
      lmin = min(thd[lmax1:lmax2][lmin])
      }else{
        lmin =localMaxima(-c(thd[lmax2:n],thd[1:lmax1]))
        lmin = min(c(thd[lmax2:n],thd[1:lmax1])[lmin])
      }
      pks=c(thd[lmax[1]],thd[lmax[2]],lmin)
    }else{bdscore=0;degdiff=0;pks=c(0,0,0)}
    return(c(bdscore,degdiff,pks))
  }
  #################################################################################################
  print(rs@path)
  
  ## load the data in R
  myList<-getRecSessionObjects(rs)
  pt<-myList$pt
  cg<-myList$cg
  sp<-myList$sp
  hd<-myList$hd
  #################################################################################################
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
  ## create some empty data.frames to store the data
  
  ## select data when the animal ran between 3 and 100
  ptsqr70<-speedFilter(pt,minSpeed=3,maxSpeed = 100)
  xbins=40;  ybins=40
  #################################################################################################

  st<-myList$st
  lightInt=lightInt[apply(lightInt[,2:3],1,diff)>=118*20000,]
  bdscore.hd<-data.frame()
  
  for(condition in c("l1","l2")){
    m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition==condition),2], 
                                lightInt[which(lightInt$condition==condition),3])),ncol=2)
    
    st0<-setIntervals(st,s=m);
    
    ## get the hd properties during the sqr70 in order to identify the neurons
    hd@smoothRateHistoSd<-20
    hd0<-headDirectionHisto(hd,st0,ptsqr70)
    #hdSt<-headDirectionStats(hd0,st0,ptsqr70)

    thd=hd0@histo
    bds=    apply(thd, 2, get_bdscore)[1,]
    degdiff=apply(thd, 2, get_bdscore)[2,]
    pks1=   apply(thd, 2, get_bdscore)[3,]
    pks2=   apply(thd, 2, get_bdscore)[4,]
    pksmin= apply(thd, 2, get_bdscore)[5,]
    

    
    if(condition=="l1"){
    bdscore.hd1<-data.frame(clu.id=cg@id,
                            bdscore=bds,degdiff=degdiff,
                            peak.r1=pks1,peak.r2=pks2,
                            peak.min=pksmin)
    }
    if(condition=="l2"){
      bdscore.hd2<-data.frame(clu.id=cg@id,
                              bdscore=bds,degdiff=degdiff,
                              peak.r1=pks1,peak.r2=pks2,
                              peak.min=pksmin)
    }

  }
  ##### shuffling procedure ######################################################################
  # bdscore.hd.shuffle=data.frame()
  # n=500
  # for(i in 1:n)
  # {
  #   lightInt$condition<-sample(lightInt$condition) ## shuffle the conditions
  #   ## get intervals for l1
  #   m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l1"),2],
  #                               lightInt[which(lightInt$condition=="l1"),3])),ncol=2)
  #   st0<-setIntervals(st,s=m)
  #   ## get the hd properties during the sqr70 in order to identify the neurons
  #   hd0<-headDirectionHisto(hd,st0,ptsqr70)
  #   #hdSt<-headDirectionStats(hd0,st0,ptsqr70)
  # 
  #   thd=hd0@histo#[,1]
  #   #thd[thd==-1]=0
  #   if(length(cid0)==1){
  #     bds1=get_bdscore(thd)[1]
  #     degdiff1=get_bdscore(thd)[2]
  #   }else{
  #     bds1=apply(thd, 2, get_bdscore)[1,]
  #     degdiff1=apply(thd, 2, get_bdscore)[2,]
  #   }
  # 
  #   #################################################################################################
  #   m<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l2"),2],
  #                               lightInt[which(lightInt$condition=="l2"),3])),ncol=2)
  #   st0<-setIntervals(st,s=m)
  #   ## get the hd properties during the sqr70 in order to identify the neurons
  #   hd@smoothRateHistoSd<-20
  #   hd0<-headDirectionHisto(hd,st0,ptsqr70)
  #   #hdSt<-headDirectionStats(hd0,st0,ptsqr70)
  # 
  #   thd=hd0@histo#[,1]
  #   #thd[thd==-1]=0
  #   if(length(cid0)==1){
  #     bds2=get_bdscore(thd)[1]
  #     degdiff2=get_bdscore(thd)[2]
  #   }else{
  #     bds2=apply(thd, 2, get_bdscore)[1,]
  #     degdiff2=apply(thd, 2, get_bdscore)[2,]
  #   }
  #   ############################################
  # 
  #   bdscore.hd.shuffle<-rbind(bdscore.hd.shuffle,data.frame(clu.id=cid0,
  #                                                        bdscore=bds1-bds2,
  #                                                        degdiff=180-abs(abs(degdiff1 - degdiff2) - 180)))
  # }
   #return(list(bdscore.hd.shuffle=bdscore.hd.shuffle,bdscore.hd=bdscore.hd))
  return(list(bdscore.hd1=bdscore.hd1,bdscore.hd2=bdscore.hd2))
}
