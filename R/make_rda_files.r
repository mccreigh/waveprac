## basic steps for creating the .rda files from text files. 

## lees ferry annual flows, units of acre*ft
lf=as.data.frame(matrix(scan("LF.txt"), ncol=2, byrow=T))
names(lf) <- c("year","annual.acreft")
lf$annual.acreft= lf$annual.acreft - mean(lf$annual.acreft)
lees.ferry <- lf
save(lees.ferry, file='~/R/jlm_lib/waveprac/data/lees.ferry.rda')

#####################################################################
## nino3 example
nino3=as.data.frame(matrix(scan("sst_nino3.dat"), ncol=1, byrow=T))
names(nino3) <- c("nino3")
nino3$year <- 1871.0 + (0:(length(nino3$nino3)-1))*.25
save(nino3, file='~/R/jlm_lib/waveprac/data/nino3.rda')
