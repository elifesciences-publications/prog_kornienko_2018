get_frequency_spectrum<-function(rs){
  print(rs@session)
  myList<-getRecSessionObjects(rs)
  st<-myList$st
  pt<-myList$pt
  cg<-myList$cg

  wf=c();wfid=c()
  m<-getIntervalsAtSpeed(pt,5,100)
  for (cc in 1:length(cg@id)) {
    st<-myList$st
    st<-setCellList(st,cc+1)
  ##########################################################
    st1<-setIntervals(st,s=m)
    st1<-ifr(st1,kernelSdMs = 5,windowSizeMs = 2)
    Fs=1000/2
    x=st1@ifr[1,]
    xts <- ts(x, frequency=Fs)
    w <- oce::pwelch(xts,nfft=512*2, plot=FALSE,log="no")
    wf0=w$spec*w$freq
    wf=rbind(wf,wf0)
    wfid=cbind(wfid,cg@id[cc])
  }
  return(list(spectrum=t(wf),spectrum.id=wfid,spectrum.freq=w$freq))
  }

