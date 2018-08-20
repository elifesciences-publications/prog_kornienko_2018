## identify cells with significant changes between the two light conditions
## in mean firing rate, map similarity, preferred direction or HD score

load(paste(ep@resultsDirectory,"light.changes",sep="/"))
load(paste(ep@resultsDirectory,"light.changes.shuffle",sep="/"))

lc <-light.changes
lcs<-light.changes.shuffle
prob=c(0.01,0.99)
sig.cells<-data.frame()

for (cellid in lc$cell.id) 
{
  x<-lc[which(lc$cell.id==cellid),]
  xs<-lcs[which(lcs$cell.id==cellid),]
  
  # mean firing rate
  xxs<-xs$mean.rate.diff;q.xxs<-quantile(xxs,prob,na.rm = T) 
  xx<-x$mean.rate.diff
  sig.fr<-as.numeric(xx<q.xxs[1] |  xx>q.xxs[2])+0
  
  # map similarity
  xxs<-xs$map.cor;q.xxs<-quantile(xxs,prob,na.rm = T) 
  xx<-x$map.cor
  sig.mapcor<-as.numeric(xx<q.xxs[1])+0
  
  # preferred direction
  xxs<-xs$mean.dir;q.xxs<-quantile(xxs,prob,na.rm = T) 
  xx<-x$mean.dir
  sig.hd.dir<-as.numeric(xx>q.xxs[2])+0
  
  # HD score
  xxs<-xs$vectorLength.diff;q.xxs<-quantile(xxs,prob,na.rm = T) 
  xx<-x$vectorLength.diff
  sig.vectorLength<-as.numeric(xx<q.xxs[1] |  xx>q.xxs[2])+0

  sig.cells<-rbind(sig.cells,c(sig.fr,sig.mapcor,sig.hd.dir,sig.vectorLength))
}
sig.cells<-cbind(lc$cell.id,sig.cells)
colnames(sig.cells)<-c("cell.id","fr","mapcor","hd.dir","hd.score")

write.table(sig.cells,file=paste(ep@resultsDirectory,"cells.sig.changes",sep="/"),sep="\t")