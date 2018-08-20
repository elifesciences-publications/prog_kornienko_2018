
# load(paste(ep@resultsDirectory,"cells1.test",sep="/"))
# load(paste(ep@resultsDirectory,"cells2.test",sep="/"))
# load(paste(ep@resultsDirectory,"shuf1.test",sep="/"))
# load(paste(ep@resultsDirectory,"shuf2.test",sep="/"))
# shuf1=shuf2.test
# shuf2=shuf2.test
# cells1=cells1.test
# cells2=cells2.test

load(paste(ep@resultsDirectory,"cells1",sep="/"))
load(paste(ep@resultsDirectory,"cells2",sep="/"))
load(paste(ep@resultsDirectory,"shuf1",sep="/"))
load(paste(ep@resultsDirectory,"shuf2",sep="/"))

load(paste(ep@resultsDirectory,"DR.id",sep="/"))
load(paste(ep@resultsDirectory,"DR1",sep="/"))
load(paste(ep@resultsDirectory,"DR2",sep="/"))


cells.class<-data.frame(matrix(NA, nrow = dim(cells1)[1], ncol = 3))
colnames(cells.class)<-c("cell.id","grid","hd")
cells.class$cell.id<-cells1$cell.id 


## significance levels
prob=0.99
##################################################################################
t.grid.score<-quantile(c(shuf1$grid.score,shuf2$grid.score),prob,na.rm=T)
t.info.score<-max(quantile(c(shuf1$info.score,shuf2$info.score),prob,na.rm=T)) 


print(paste("threshold from shuffling with a probability of",prob))
print(paste("GridScore threshold:",round(t.grid.score,3)))
print(paste("InfoScore threshold:",round(t.info.score,3)))



### grid cells
#Neurons with a significant grid score & significant info score in l1 or l2
cells.class$grid<-(cells1$grid.score>t.grid.score | cells2$grid.score>t.grid.score) &
 (cells1$info.score>t.info.score & cells2$info.score>t.info.score)
cells.class$grid<-cells.class$grid+0
print(paste("Number of grid cells (above GridScore & InfoScore threshold), N =",sum(cells.class$grid,na.rm=T)))

### hd cells
# Head direction cells had a significant vector 
# length and a peak firing rate above 5Hz 
cells.class$hd<-((cells1$hd.vl>0.4 & cells1$hd.peak.rate>5) | (cells2$hd.vl>0.4 & cells2$hd.peak.rate>5))+0
tmp<-(cells.class$hd==1 & cells.class$grid!=1)
cells.class$hd[DR1<.2 & DR2<.2 & tmp]<- 0

print(paste("Number of head-direction cells, HDscore > 0.4 and PeakRate > 5 Hz, N =",sum(cells.class$hd==1 & cells.class$grid!=1,na.rm = T)))
print(paste("Number of conjunctive grid x head-direction cells (excluded), N =",sum(cells.class$hd==1 & cells.class$grid==1,na.rm = T)))
print(paste("Number of head-direction cells with directional distributive ratio < 0.2 (excluded), N =",sum(DR1[tmp]<.2 & DR2[tmp]<.2)))

cells.class$hd[cells.class$hd==1 & cells.class$grid==1]<- 0

##################################################################################


rm(shuf1,shuf2,prob,DR1,DR2,tmp,t.grid.score,DR.id)

write.table(cells.class,file=paste(ep@resultsDirectory,"cells.class.table",sep="/"),sep="\t")
      








