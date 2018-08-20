## Spike wave form, peak amplitude asymmetry, trough-to-peak duration
#################################################################
## Average spike wave forms from 3-4 recording channels of 
## N = 93 HD cells extracted from raw data

#load("~/repo/Waveforms.Rdata")
load(paste(ep@resultsDirectory,"Waveform.RData",sep=""))

#################################################################
#################################################################
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
#################################################################
## TEST DIFFERENCES FOR THETA RHYTHMIC & NON-RHYTHMIC CELLS
#################################################################

ii1=which(theta.i[t$hd==1]<0.07)
ii2=which(theta.i[t$hd==1]>0.07)
# par(mfrow=c(3,4))
# for (ii in ii1) {
#   plot(wfs[[ii]][,1],type="l",ylim=c(min(as.vector(wfs[[ii]])),max(as.vector(wfs[[ii]]))),
#        main = round((theta.i[t$hd==1])[ii],2),xlab=ii) #,ylab=cid[ii]
#   lines(wfs[[ii]][,2],col="red")
#   lines(wfs[[ii]][,3],col="blue")
#   try(lines(wfs[[ii]][,4],col="green"))
# }

# identify channel with largest spike waveform deflection 
wf1=sapply(wfs[ii1], function(x){x[,which.min(apply(x, 2, min))]})
wf2=sapply(wfs[ii2], function(x){x[,which.min(apply(x, 2, min))]})

# exclude inverted waveforms 23, 25, 27, 33, 35, 39, 40 
wf1=wf1[,!(1:dim(wf1)[2]%in%c(3, 5, 7, 30))]
wf2=wf2[,!(1:dim(wf2)[2]%in%c(23, 25, 27, 33, 35, 39,40))]

##############################################################################
# trough-to-peak duration AND peak amplitude asymmetry of NON-RHYTHMIC HD CELLS

par(mfrow=c(3,4))
duration1=c()
sasym1=c()
wf1s=c()
tt=seq(-1.475,1.475,.05)
time=seq(-1,1,0.01)


for (ii in 1:dim(wf1)[2]) {
  # approximate and smooth waveform on finer time scale
  y0=approx(x=tt,y=wf1[,ii],xout = time)
  voltage=smooth.spline(y0$y,spar=.3)$y
  # detect spike time (around 0 ms), and second peak in spike waveform
  spike_min=which.min(voltage)
  spike_max2=find_peaks(voltage[spike_min:length(time)])[1]+spike_min-1
  # detect first peak in spike waveform
  n1=((length(time)-1)/4)
  spike_max1=tail(find_peaks(voltage[n1:spike_min]),n=1)+(n1-1)
  # if no clear first peak take saddle point
  if(length(spike_max1)==0){
    spike_max1 = tail(which(diff(diff(voltage[n1:spike_min]))>0),n=1)+(n1-1)
    }
  
  # calculate peak amplitude asymmetry
  a = voltage[spike_max1]
  b = voltage[spike_max2]  
  sasym1=c(sasym1,(b-a)/(abs(b)+abs(a)))
  # calculate trough-to-peak duration
  duration1 =c(duration1,time[spike_max2]-time[spike_min]) 
  wf1s=rbind(wf1s,voltage)
  
  # plot(time,voltage,type="l", main=paste("t=",round(time[spike_max2],2)," msec, asym=",round((b-a)/(abs(b)+abs(a)),2),sep=""),ylab=ii)
  # lines(c(time[spike_max1],time[spike_max1]),c(min(voltage),max(voltage)),col="red")
  # lines(c(time[spike_max2],time[spike_max2]),c(min(voltage),max(voltage)),col="green")
  
}


############################################################
# B) RHYTHMIC HD CELLS
duration2=c()
sasym2=c()
wf2s=c()

for (ii in 1:dim(wf2)[2]) {
  y0=approx(x=tt,y=wf2[,ii],xout = time)
  voltage=smooth.spline(y0$y,spar=.3)$y
  # detect spike time (around 0 ms), and second peak in spike waveform
  spike_min=which.min(voltage)
  spike_max2=find_peaks(voltage[spike_min:length(time)])[1]+spike_min-1
  # detect first peak in spike waveform
  n1=((length(time)-1)/4)
  spike_max1=tail(find_peaks(voltage[n1:spike_min]),n=1)+(n1-1)
  # if no clear first peak take saddle point
  if(length(spike_max1)==0){
    spike_max1= tail(which(diff(diff(voltage[n1:spike_min]))>0),n=1)+(n1-1)
    }
  
  # calculate peak amplitude asymmetry
  a = voltage[spike_max1]
  b = voltage[spike_max2]  
  sasym2=c(sasym2,(b-a)/(abs(b)+abs(a)))
  # calculate trough-to-peak duration
  duration2 =c(duration2,time[spike_max2]-time[spike_min]) 
  wf2s=rbind(wf2s,voltage)
  
  # plot
  # plot(time,voltage,type="l",main=paste("t=",round(time[spike_max2],2),
  #                                      " msec, asym=",round((b-a)/(abs(b)+abs(a)),2),sep=""),ylab=ii)
  # lines(c(time[spike_max1],time[spike_max1]),c(min(voltage),max(voltage)),col="red")
  # lines(c(time[spike_max2],time[spike_max2]),c(min(voltage),max(voltage)),col="green")
}

###################################################################################
# PLOT SPIKE WAVEFORM

save(file=paste(ep@resultsDirectory,"/Waveform.stats.RData",sep=""),wf1s,wf2s,time,duration1,duration2,sasym1,sasym2)

