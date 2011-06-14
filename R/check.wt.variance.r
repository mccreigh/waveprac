## name: check.wt.variance
## author: james mccreight = mccreigh ^ gmail * com
## info:
##   part of "A practical guide to wavelet analysis" for R,
##   adapted from Torrence & Compo:
##   http://paos.colorado.edu/research/wavelets/

## General Purpose:
##  Perform some tests on wavelet transforms to verify preservation
##  of variance. Parseval's theorem and reconstruction are used.

## Return Value:
##   List of 6
##    variance         : num : original time series variance
##    parseval.variance: num 
##    pe.parceval      : num  percent error in parceval variance
##    recon.variance   : num 
##    pe.recon         : num  percent error in reconstructed timeseries variance
##    rmse.recon       : num  root mean square error in reconstructed timeseries.

## Inputs
##   wt - a wavelet transform object
##   slient=FALSE Prints all the returned information with units by default.

check.wt.variance <- function(wt, silent=FALSE) {

  ntime=length(a.POSIXct(a.input(wt)))
  nscale=length(a.scale(wt))

  variance = var(a.data(a.input(wt)))
  if (!silent) {
    print("")
    print(paste(a.data.name(a.input(wt)), ' variance: ', variance, a.data.units(a.input(wt)),'^2',sep=''))
  }
  
  ## Parseval's theorem [Eqn(14)] - check variance
  scale.avg = matrix(a.scale(wt), ncol=ntime, nrow=nscale )  ## expand scale-->(J+1)x(N) array
  power.norm = abs(a.transform(wt))^2/scale.avg
  parseval.variance = a.dj(wt)*a.dt(wt)/(a.Cdelta(wt)*ntime)*sum(power.norm)  ## [Eqn(14)]

  if (!silent) {
    print(paste("Parseval variance: ", parseval.variance, a.data.units(a.input(wt)),'^2',sep='')) 
    print(paste("Percent error in Parseval variance: ", (1-(variance/parseval.variance))*100,sep=''))
  }
  
  ## Reconstruction - check variance and RMS errors
  ## if there's no reconstruction in the wt, make ie
  if (typeof(a.recon(wt))=='logical') {
    if (a.Cdelta(wt) == -1) {
      recon = -1
      warning('Cdelta undefined, cannot reconstruct with this wavelet.')
    } else {
      recon=a.dj(wt)*sqrt(a.dt(wt))/(a.Cdelta(wt)*a.psi0(wt))*(t(Re(a.transform(wt))) %*% (1./sqrt(a.scale(wt)))) 
      recon = recon[1:ntime]
    }
  } else recon=a.recon(wt) 

  recon.variance = var(recon)
  rmse = sqrt(sum((a.data(a.input(wt)) - recon)^2)/ntime) ## RMS of Reconstruction [Eqn(11)]
  
  if (!silent) {
    print(paste("Reconstruction variance: ", recon.variance, a.data.units(a.input(wt)),'^2',sep='')) 
    print(paste("Percent error in reconstruction variance: ", (1-(variance/recon.variance))*100,sep=''))
    print(paste("RMS Error in reconstruction timeseries: ", rmse,' ', a.data.units(a.input(wt)), sep=''))
    if (a.mother(wt) == 'dog') print('Note: for better reconstruction with the DOG, you need to use a very small s0.')
    print("")
  }

  invisible( list( variance=variance, parseval.variance=variance, pe.parceval=(1-(variance/parseval.variance))*100,
                   recon.variance=recon.variance, pe.recon=(1-(variance/recon.variance))*100, rmse.recon=rmse )
            )
  
}
