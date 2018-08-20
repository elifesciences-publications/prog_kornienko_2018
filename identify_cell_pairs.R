t<-read.table(paste(ep@resultsDirectory,"cells.class.table",sep="/"),sep="\t")
load(paste(ep@resultsDirectory,"cells1",sep="/"))

#####################################################################################################
# IDENTIFY HD CELL PAIRS
a<-sapply(rss,function(x){x@session})
A<-cbind(as.character(cells1$session),as.character(cells1$cell.id))
A<-A[t$hd==1,]
f<-(table(A[,1]));f<-names(f[f>1])

cidl<-list()
for (i in 1:length(f)){
  cells<- A[f[i]==A[,1],2]
  if(length(cells)<15){
    cells<-c(cells, rep(" ",15-length(cells)))}
  cidl[[i]]<-cells
}
id.hd.pairs=c()
for(i in 1:length(cidl)) {
  cc=cidl[[i]];
  cc=cc[cc!=" "]
  if(length(cc)==2){ 
    id.hd.pairs=rbind(id.hd.pairs,cc)
  }else{
    for (j in 1:(length(cc)-1)) {
      for (jj in (j+1):(length(cc))) {
        id.hd.pairs=rbind(id.hd.pairs,c(cc[j],cc[jj]))
      }
    }
  }
}
rownames(id.hd.pairs)=NULL
cid<-unlist(cidl);cid<-cid[cid!=" "]
sid<-sapply(cid,function(x){strsplit(x,"_")[[1]][1]})
rss.hd<-rss[which((a %in% unique(sid)))]

######################################################################################################
# IDENTIFY GRID CELL PAIRS

A<- cbind(as.character(cells1$session),as.character(cells1$cell.id))
A<-A[t$grid==1 ,]
f<-(table(A[,1]));f<-names(f[f>1])

cidl<-list()
for (i in 1:length(f)){
  cells<- A[f[i]==A[,1],2]
  if(length(cells)<15){
    cells<-c(cells, rep(" ",15-length(cells)))}
  cidl[[i]]<-cells
}

id.grid.pairs=c()
for(i in 1:length(cidl)) {
  cc=cidl[[i]];
  cc=cc[cc!=" "]
  if(length(cc)==2){ 
    id.grid.pairs=rbind(id.grid.pairs,cc)
  }else{
    for (j in 1:(length(cc)-1)) {
      for (jj in (j+1):(length(cc))) {
        id.grid.pairs=rbind(id.grid.pairs,c(cc[j],cc[jj]))
      }
    }
  }
}
rownames(id.grid.pairs)=NULL
cid<-unlist(cidl);cid<-cid[cid!=" "]
sid<-sapply(cid,function(x){strsplit(x,"_")[[1]][1]})
rss.grid<-rss[which((a %in% unique(sid)))]

rm(cid,sid,cc,cells,A,a,f)
