### function doing all the work on individual sessions
get_stime_autocorrelation<-function(rs,windowSize,binSize){
  print(paste(rs@session,rs@path))
  ## load the data in R
  myList<-getRecSessionObjects(rs)
  st<-myList$st
  cg<-myList$cg
  st<-spikeTimeAutocorrelation(st,binSizeMs = binSize,windowSizeMs = windowSize,probability = F)
  st<-spikeTimeCrosscorrelation(st,binSizeMs = binSize,windowSizeMs = windowSize,probability = F)
  st.auto<-(st@auto);st.cross<-st@cross;st.pair.list=data.frame(rs@session,st@cellPairList);
  return(list(st.auto=st.auto,st.cross=st.cross,st.pair.list=st.pair.list))
}
  
