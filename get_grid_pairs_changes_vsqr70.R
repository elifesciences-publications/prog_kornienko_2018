### function doing all the work on individual sessions
grid.pairs.changes.vsqr70<-function(rs,cid){
  print(rs@path)
  ## load the data in R
  myList<-getRecSessionObjects(rs)
  pt<-myList$pt
  cg<-myList$cg
  sp<-myList$sp
  hd<-myList$hd
  st<-myList$st
  sid<-unlist(sapply(cid,strsplit,"_",simplify = T,USE.NAMES = F))[seq(1,2*length(cid),2)]
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
  ## select data when the animal ran between 3 and 100
  ptsqr70<-speedFilter(pt,minSpeed=3,maxSpeed = 100)
  xbins=40;  ybins=40
  #################################################################################################
  ## SET CELL LIST
  cid0<- cid[which(rs@session==sid)]
  id<-which(cg@id%in%cid0)+1
  st<-setCellList(st,id)
  lightInt=lightInt[apply(lightInt[,2:3],1,diff)>=118*20000,]
  ########## get intervals ########################################################################
  m1<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l1"),2], 
                               lightInt[which(lightInt$condition=="l1"),3])),ncol=2)
  m3<-matrix(data=as.numeric(c(lightInt[which(lightInt$condition=="l2"),2], 
                               lightInt[which(lightInt$condition=="l2"),3])),ncol=2)
  n1=round(dim(m1)[1]/2)
  s=sample(dim(m1)[1]);
  s1=sort(s[1:n1]);
  s2=sort(s[(n1+1):dim(m1)[1]]);
  s3=sort(sample(dim(m3)[1])[1:n1])
  
  m2<-m1[s2,];m1<-m1[s1,];m3<-m3[s3,]
  
  st1<-setIntervals(st,s=m1)
  st2<-setIntervals(st,s=m2)
  st3<-setIntervals(st,s=m3)
  ######## MAPs  ########################################################################
  sp1<-getMapStats(sp,st1,ptsqr70)
  sp1<-firingRateMap2d(sp1,st1,ptsqr70,nRowMap = xbins,nColMap = ybins)
  sp2<-getMapStats(sp,st2,ptsqr70)
  sp2<-firingRateMap2d(sp2,st2,ptsqr70,nRowMap = xbins,nColMap = ybins)
  sp3<-getMapStats(sp,st3,ptsqr70)
  sp3<-firingRateMap2d(sp3,st3,ptsqr70,nRowMap = xbins,nColMap = ybins)
  maps1=sp1@maps;maps2=sp2@maps; maps3=sp3@maps;
  
  maps1[maps1==-1]=NA;maps2[maps2==-1]=NA;maps3[maps3==-1]=NA;
  ######## ifr ########################################################################
  st1<-ifr(st1,windowSizeMs = 100, spikeBinMs = 1,kernelSdMs = 200)#10000 
  st2<-ifr(st2,windowSizeMs = 100, spikeBinMs = 1,kernelSdMs = 200)#10000 
  st3<-ifr(st3,windowSizeMs = 100, spikeBinMs = 1,kernelSdMs = 200)#10000 
  st1<-meanFiringRate(st1);st2<-meanFiringRate(st2);st3<-meanFiringRate(st3)
  
  ifr1=st1@ifr;ifr2=st2@ifr;ifr3=st3@ifr
  #######################################################################################
  cor1=c();cor2=c();cor3=c();
  ifrcor1=c();ifrcor2=c();ifrcor3=c();

  pairs.id=c()
  if(length(id)==2){
    pairs.id=cg@id[id-1]
    #maps
    cor1=cor(as.vector(maps1[,,1]),as.vector(maps1[,,2]),use="pairwise.complete.obs")
    cor2=cor(as.vector(maps2[,,1]),as.vector(maps2[,,2]),use="pairwise.complete.obs")
    cor3=cor(as.vector(maps3[,,1]),as.vector(maps3[,,2]),use="pairwise.complete.obs")
    #ifr
    ifrcor1=cor(ifr1[1,],ifr1[2,],use="pairwise.complete.obs")
    ifrcor2=cor(ifr2[1,],ifr2[2,],use="pairwise.complete.obs")
    ifrcor3=cor(ifr3[1,],ifr3[2,],use="pairwise.complete.obs")
  }else{
    for(i in 1:(length(id)-1)){
      for(j in (i+1):length(id)){
        pairs.id=cbind(pairs.id,cg@id[c(id[i]-1,id[j]-1)])
        #maps
        cor1=c(cor1,cor(as.vector(maps1[,,i]),as.vector(maps1[,,j]),use="pairwise.complete.obs"))
        cor2=c(cor2,cor(as.vector(maps2[,,i]),as.vector(maps2[,,j]),use="pairwise.complete.obs"))
        cor3=c(cor3,cor(as.vector(maps3[,,i]),as.vector(maps3[,,j]),use="pairwise.complete.obs"))
        #ifr
        ifrcor1=c(ifrcor1,cor(ifr1[i,],ifr1[j,],use="pairwise.complete.obs"))
        ifrcor2=c(ifrcor2,cor(ifr2[i,],ifr2[j,],use="pairwise.complete.obs"))
        ifrcor3=c(ifrcor3,cor(ifr3[i,],ifr3[j,],use="pairwise.complete.obs"))
     
      }
    }
  }
  pairs.ifra<- rbind(t(ifrcor1),t(ifrcor2),t(ifrcor3)) #rows: l1.1,l1.2,l2; cols:cell pairs
  pairs.map<-rbind(t(cor1),t(cor2),t(cor3)) 

  return(list(grid.pairs.ifra=pairs.ifra,
              grid.pairs.map=pairs.map,
              grid.pairs.id=as.matrix(pairs.id)))
}