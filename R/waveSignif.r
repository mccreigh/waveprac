## name: waveSignif method (on wavelet transform object)
## author: james mccreight = mccreigh ^ gmail * com
## info:
##   adapted from Torrence and Compo: A Practical Guide to Wavelet
##   analysis. see license and credit below.

## General Purpose:
## Compute the significance levels for a wavelet transform,
## including wavelet transform, time-average (global wavelet)
## spectrum, and scale-average spectrum significances.

## Return Value:
## A list with the following entries:
##   List of 13
##    signif    : num : vector of significance as a function of scale. A single value for scale average significance.
##    confidence: num : matrix of confidence interval (min,max) corresponding to signif. 
##    ratio     : NULL or num: for wavlet transform, a ratio of the wavelet transform to the significance at each scale.
##    sigtest   : num : which test was asked for.
##    siglvl    : num : the significance level asked for.
##    dof       : num : the degrees of freedom.
##    mother    : chr : the mother wavlet.
##    fft.theor : num : the theoretical fourier spectrum
##    gws       : logi : if the theoretical fourier specturm was set to the global wavelet spectrum
##    lag1      : num : the alpha parameter in eqn 16.
##    Cdelta    : num : Cdelta reconstruction param for the wavelet.
##    Savg      : num :parameter used in scale-averaging
##    Smid      : num :parameter used in scale-averaging

## Required inputs:
##    wt - a wavelet transform object.
##    sigtest = 0, 1, or 2. Default is 0.
##           0, then just do a regular chi-square test,
##       		 i.e. Eqn (18) from Torrence & Compo.
##           1, then do a "time-average" test, i.e. Eqn (23).
##       		 In this case, DOF should be set to NA, the number
##       		 of local wavelet spectra that were averaged together.
##       		 For the Global Wavelet Spectrum, this would be NA=N,
##       		 where N is the number of points in your time series.
##           2, then do a "scale-average" test, i.e. Eqns (25)-(28).
##       		 In this case, DOF should be set to a
##       		 two-element vector [S1,S2], which gives the scale
##       		 range that was averaged together.
##       		 e.g. if one scale-averaged scales between 2 and 8,
##                       then DOF=[2,8].

## Optional inputs:
##    lag1=0.0  - parameter in eqn 16, used for calculating red noise. 0 gives white noise.
##    siglvl=.95  - the desired significance.
##    dof=dof     - the degrees of freedom.
##    gws=FALSE   - use global wavelet spectrum as theoretical fourier spectrum.
##    fft.theor = (1-lag1^2)/(1-2*lag1*cos(freq*2*pi)+lag1^2) ## eqn 16, the default theoretical spectrum. you could supply your own.

##----------------------------------------------------------------------------
## Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo,
## University of Colorado, Program in Atmospheric and Oceanic Sciences.
## This software may be used, copied, or redistributed as long as it is not
## sold and this copyright notice is reproduced on each copy made.  This
## routine is provided as is without any express or implied warranties whatsoever.
##
## Notice: Please acknowledge the use of the above software in any publications:
##    ``Wavelet software was provided by C. Torrence and G. Compo,
##      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
##
##----------------------------------------------------------------------------

waveSignif <- function(wt, sigtest=sigtest,    ##*** required inputs
                       lag1=0.0, siglvl=.95,
                       dof=dof,   
                       gws=FALSE,
                       fft.theor = (1-lag1^2)/(1-2*lag1*cos(freq*2*pi)+lag1^2) ## eqn 16
                       )
{
    
  ## time series variance
  variance <- sd(a.data(a.input(wt)))^2
    
  ## default NULL values for these output items
  ratio <- NULL; Savg <- NULL; Smid <- NULL
   
  ## I could just calculate this.... does ken have red and white options?..
  lag1 = lag1[1]
  
  ##############################
  ## Wavelet specific setup
  J = length(a.scale(wt))-1
  J1 = J+1  ## this is silly, a hold over from IDL code.
  nscale=J1; ntime=length(a.POSIXct(a.input(wt)))
  s0 = a.s0(wt) ## min(a.scale(wt))
  dj = a.dj(wt) ## log(a.scale(wt)[2]/a.scale(wt)[1])/log(2)
  dt= a.dt(wt) ## not defined anywhere but used... getting from parent level! yikes.
  
  if (a.mother(wt)=='morlet') {
    if (a.param(wt) == -1) k0=6 else k0=a.param(wt)
    fourier_factor = (4*pi)/(k0 + sqrt(2+k0^2)) ## [Sec.3h]
    empir = c(2.,-1,-1,-1)
    if (k0 == 6) empir[2:4]=c(0.776,2.32,0.60)
  }

  if (a.mother(wt)=='paul') {
    if (a.param(wt) == -1) m=4 else m=a.param(wt)
    fourier_factor = 4*pi/(2*m+1)
    empir = c(2.,-1,-1,-1)
    if (m == 4) empir[2:4]=c(1.132,1.17,1.5)
  }

  if (a.mother(wt)=='dog') {
    if (a.param(wt) == -1) m=2 else m=a.param(wt)
    fourier_factor = 2*pi*sqrt(2./(2*m+1))
    empir = c(1.,-1,-1,-1)
    if (m == 2) empir[2:4] = c(3.541,1.43,1.4)
    if (m == 6) empir[2:4] = c(1.966,1.37,0.97)
  }

  ## first three should already exist in the wt object... somewhat redundant.
  period = a.scale(wt)*fourier_factor 
  dofmin = empir[1] ## Degrees of freedom with no smoothing
  Cdelta = empir[2] ## reconstruction factor 
  gamma = empir[3]  ## time-decorrelation factor
  dj0 = empir[4]    ## scale-decorrelation factor

  ##############################
  ## set up significance params
  freq = dt/period  ## normalized frequency
  ## Assign a theoretical background fourier spectrum. Mostly done in defaults above.
  ## default: fft.theor = (1-lag1^2)/(1-2*lag1*cos(freq*2*pi)+lag1^2)  ## [Eqn(16)], red-noise spectrum base on AR1 process.
  fft.theor <- variance*fft.theor
  ## If use global wavelet specturm as theoretical background, then calculate it.
  if (gws) fft.theor <- rowMeans( abs(a.transform(wt))^2 )

  ## should check the fft.theor (now since the promise can be kept)
  if (length(fft.theor) != length(a.scale(wt)))
    warning("The theoretical fourier spectrum is not the appropriate dimension.", immediate.=TRUE)

  ##############################
  ## wavelet transform w/o averaging
  if (sigtest==0) {   ## no smoothing, DOF=dofmin
    dof = dofmin
    ## signif
    signif = fft.theor*qchisq(1. - siglvl, dof, lower.tail=FALSE)/dof 
    ## confidence
    sig = (1. - siglvl)/2.
    chisqr = dof/ c( qchisq(sig, dof, lower.tail=FALSE),  qchisq(1-sig, dof, lower.tail=FALSE) )
    conf = as.matrix(fft.theor) %*% chisqr
    ## power/signif ratio, >1 means significant power.
    ratio <- abs(a.transform(wt))^2 / matrix( signif, ncol=ntime, nrow=nscale )   
  }

  ##############################
  ## time-averaged, DOFs depend upon scale [Sec.5a]
  if (sigtest==1) {   
    if (missing(dof)) dof = dofmin
    if (gamma == -1) warning(paste('Gamma (decorrelation factor) not defined for ',a.mother(wt),' with param=',a.param(wt), sep=''))
    if (length(dof) == 1) dof = (0:J) + dof
    dof = pmax(dof , 1)
    dof = dofmin*sqrt( 1 + (dof*dt/gamma/a.scale(wt))^2 ) ## [Eqn(23)]
    dof = pmax(dof , dofmin)   ## minimum DOF is dofmin
    signif = fft.theor
    ## signif
    for (a1 in 1:J1) {
      chisqr = qchisq( 1. - siglvl, dof[a1], lower.tail=FALSE)/dof[a1]  
      signif[a1] = fft.theor[a1]*chisqr
    }
    ## conf
    conf = matrix(NA, ncol=J1, nrow=2)
    sig = (1. - siglvl)/2.
    for ( a1 in 1:J1 ) {
      chisqr = dof[a1]/ c( qchisq(sig, dof[a1], lower.tail=FALSE),  qchisq(1-sig, dof[a1], lower.tail=FALSE) )
      conf[,a1] = fft.theor[a1]*chisqr
    }
  }

  ##############################
  ## scale-averaged, DOFs depend upon scale range [Sec.5b]
  if (sigtest==2) {  
    if (length(dof) != 2) warning('DOF must be set to [S1,S2], the range of scale-averages')
    if (Cdelta == -1) warning(paste('Cdelta & dj0 not defined for ',a.mother(wt),' with param=',a.param(wt),sep=''))
    s1 = dof[1]
    s2 = dof[2]
    avg = which((a.scale(wt) >= s1) & (a.scale(wt) <= s2))
    navg=length(avg)
    if (navg < 1) warning(paste('No valid scales between ',s1,' and ',s2,'.',sep=''))
    s1 = min(a.scale(wt)[avg])
    s2 = max(a.scale(wt)[avg])
    Savg = 1./sum(1./a.scale(wt)[avg])       ## [Eqn(25)]
    Smid = exp((log(s1)+log(s2))/2.)   ## power-of-two midpoint
    dof = (dofmin*navg*Savg/Smid)*sqrt(1 + (navg*dj/dj0)^2)  ## [Eqn(28)]
    fft.theor <- Savg*sum(fft.theor[avg]/a.scale(wt)[avg])   ## [Eqn(27)]
    ## signif
    chisqr = qchisq(1-siglvl, dof, lower.tail=FALSE)/dof 
    signif = (dj*dt/Cdelta/Savg)*fft.theor*chisqr  ## [Eqn(26)]
    ## conf
    sig = (1. - siglvl)/2.
    chisqr = dof/ c( qchisq(sig, dof, lower.tail=FALSE),  qchisq(1-sig, dof, lower.tail=FALSE) )
    conf = (dj*dt/Cdelta/Savg)*fft.theor*chisqr  ## [Eqn(26)]
  }  
   
  ##############################
  ## return list
  return( list(signif=signif,
               confidence=conf,
               ratio=ratio,
               sigtest=sigtest,
               siglvl=siglvl,
               dof=dof, 
               mother=a.mother(wt),
               fft.theor=fft.theor,
               gws=gws,
               lag1=lag1,               
               Cdelta=Cdelta,
               Savg=Savg,
               Smid=Smid)
         )
  
}
  
## todos
## ken red and white options?? 
