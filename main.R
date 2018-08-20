################################################################
# Code to run the analysis for the manuscript "Non-rhythmic   ##
# head-direction cells in the parahippocampal region are not  ##  
# constrained by attractor network dynamics"                  ##
# submitted at eLife, July 2018                               ##
################################################################

# load database .res and .clu files 
rm(list=ls())
library(relectro)
library(mclust)

indir<-"~/repo/hd_mec_2017_github/"
ep<-new("ElectroProject",directory="~/HD_MEC_2017_DATABASE/")
ep<-setSessionList(ep)
rss<-getSessionList(ep,clustered=T,env="sqr70")
rss<-sortRecSessionListChronologically(rss)
print(paste("Number of recording sessions",length(rss)))

################################################################
# use parallel processing otherwise set parallel=F
parallel = T
workers<-c(rep("localhost",6))
cl<-snow::makeCluster(workers, type = "SOCK",outfile="")
snow::clusterEvalQ(cl,library(relectro))


#################################################################
# #  obtain HD scores, grid scores, fix parameters
# source(paste(indir,"get_vsqr70_cell_properties.R",sep=""))
# runOnSessionList(ep,sessionList=rss,fnct=get_vsqr70_cell_properties,
#                  save=T,overwrite=T,parallel=parallel,cluster=cl)
# rm(get_vsqr70_cell_properties)

##################################################################
# calculate directional distributive ratio (Cacucci et al., 2004) 
# to ensure directional tuning is not a by-product of spatial 
# selectivity coupled with unequal HD sampling

source(paste(indir,"test_hd_selectivity_not_artifact_vsqr70.R",sep="/"))
runOnSessionList(ep,sessionList=rss,fnct=test_hd_selectivity_not_artifact_vsqr70,
                 save=T,overwrite=T,parallel=parallel,cluster=cl)
rm(test_hd_selectivity_not_artifact_vsqr70)

###################################################################
# obtain spike-time autocorrelations
windowSize<-300; binSize<-2
source(paste(indir,"get_stime_autocorrelation.R",sep="/"))
runOnSessionList(ep,sessionList=rss,fnct=get_stime_autocorrelation,
                 save=T,overwrite=T,parallel=T,cluster=cl,windowSize,binSize)
rm(get_stime_autocorrelation)  

###################################################################
# obtain power spectra from instantaneous firing rates
source(paste(indir,"get_frequency_spectrum.R",sep="/"))
runOnSessionList(ep,sessionList=rss,fnct=get_frequency_spectrum,
                 save=T,overwrite=T,parallel=parallel,cluster=cl)
rm(get_frequency_spectrum)

###################################################################
# detect significant changes between the two light conditions
source(paste(indir,"get_changes_2light_conditions.R",sep="/"))
runOnSessionList(ep,sessionList=rss,fnct=get_changes_2light_conditions,
                 save=T,overwrite=T,parallel=parallel,cluster=cl)
rm(get_changes_2light_conditions)

###################################################################
# calculate bidirectionality of HD tuning curves
source(paste(indir,"get_bdscores_vsqr70.R",sep="/"))
runOnSessionList(ep,sessionList=rss,fnct=get_bdscores_vsqr70,
                 save=T,overwrite=T,parallel=parallel,cluster=cl)
rm(get_bdscores_vsqr70)

###################################################################
# identify functional cell types: HD and grid cells
source(paste(indir,'/identify_cell_types.R',sep=""))

###################################################################
# identify simultanously recorded HD or grid cell pairs
source(paste(indir,'/identify_cell_pairs.R',sep=""))

###################################################################
# calculate pairwise changes in pref. HD, HD score and instantaneous
# firing rate association of simultanously recorded HD cell pairs

source(paste(indir,"get_hd_pairs_changes_vsqr70.R",sep="/"))
runOnSessionList(ep,sessionList=rss.hd,fnct=hd.pairs.changes.vsqr70,
save=T,overwrite=T,parallel=parallel,cluster=cl,unique(as.vector(id.hd.pairs)))
rm(hd.pairs.changes.vsqr70)

###################################################################
# calculate pairwise changes in map similarity and instantaneous
# firing rate association of simultanously recorded grid cell pairs

source(paste(indir,"get_grid_pairs_changes_vsqr70.R",sep="/"))
runOnSessionList(ep,sessionList=rss.grid,fnct=grid.pairs.changes.vsqr70,
save=T,overwrite=T,parallel=parallel,cluster=cl,unique(as.vector(id.grid.pairs)))
rm(grid.pairs.changes.vsqr70)

###################################################################
## identify cells with significant changes between the two light conditions
## in mean firing rate, map similarity, preferred direction or HD score
source(paste(indir,'/identify_significant_changes_2light_conditions.R',sep=""))

###################################################################
########################   PLOT RESULTS    ########################
###################################################################
## HD cells

## Plot 1 (Related to Figure 3)
# Distribution of theta indices, HD scores, firing rates and  
# spike waveforms of non-rhythmic and theta-rhythmic

## Plot 2 (Related to Figure 4 and Figure 5)
# pdf file with HD cell examples including spike time 
# autocorrelations, vp1/vp2 HD tuning curves, spatial firing 
# rate maps and stats on significant changes

## Plot 3 (Related to Figure 4, Figure 5 and Figure 6)
# Changes in HD score, preferred direction, firing rates and 
# BD scores between the two light conditions of non-rhythmic and 
# theta-rhythmic HD cells 

## Plot 4 (Related to Figure 7)
# Pairwise changes and reorganization scores in preferred HD,
# HD score and IFR association of simultaneously recorded HD cells

source(paste(indir,'/plot_hd_cells_vsqr70.R',sep=""))

###################################################################
## Grid cells

## Plot 5 (Related to Figure 8)
# Pdf file with grid cell examples including polar plots, 
# spatial firing rate maps and changes map  
# simililarity, HD score and rate of single grid cells

## Plot 6 (Related to Figure 8)
# change pairwise map similarity, IFR association and 
# pairwise reorganization scores  

source(paste(indir,'/plot_grid_cells_vsqr70.R',sep=""))









