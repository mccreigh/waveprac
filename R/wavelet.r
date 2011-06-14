#******************************************************************* 
#+
# NAME:   WAVELET
#
# PURPOSE:   Compute the WAVELET transform of a 1D time series.
#       
#
# CALLING SEQUENCE:
#
#      wave = WAVELET(Y,DT)
#
#
# INPUTS:
#
#    Y = the time series of length N.
#
#    DT = amount of time between each Y value, i.e. the sampling time.
#
#
# OUTPUTS:
#
#    WAVE is the WAVELET transform of Y. This is a complex array
#    of dimensions (N,J+1). FLOAT(WAVE) gives the WAVELET amplitude,
#    ATAN(IMAGINARY(WAVE),FLOAT(WAVE)) gives the WAVELET phase.
#    The WAVELET power spectrum is ABS(WAVE)^2.
#
#
# OPTIONAL KEYWORD INPUTS:
#
#    S0 = the smallest scale of the wavelet.  Default is 2*DT.
#
#    DJ = the spacing between discrete scales. Default is 0.125.
#         A smaller # will give better scale resolution, but be slower to plot.
#
#    J = the # of scales minus one. Scales range from S0 up to S0*2^(J*DJ),
#        to give a total of (J+1) scales. Default is J = (LOG2(N DT/S0))/DJ.
#
#    MOTHER = A string giving the mother wavelet to use.
#            Currently, 'Morlet','Paul','DOG' (derivative of Gaussian)
#            are available. Default is 'Morlet'.
#
#    PARAM = optional mother wavelet parameter.
#            For 'Morlet' this is k0 (wavenumber), default is 6.
#            For 'Paul' this is m (order), default is 4.
#            For 'DOG' this is m (m-th derivative), default is 2.
#
#    PAD = if 1, then pad the time series with enough zeroes to get
#         N up to the next higher power of 2. This prevents wraparound
#         from the end of the time series to the beginning, and also
#         speeds up the FFT's used to do the wavelet transform.
#         This will not eliminate all edge effects (see COI below).
#         if 2 then the original timeseries is mirriored (excluding the final point)
#         to replace the 0s above.
#
#    LAG1 = LAG 1 Autocorrelation, used for SIGNIF levels. Default is 0.0
#
#    SIGLVL = significance level to use. Default is 0.95
#
#    VERBOSE = if set, then print out info for each analyzed scale.
#
#    RECON = if set, then reconstruct the time series, and store in Y.
#            Note that this will destroy the original time series,
#            so be sure to input a dummy copy of Y.
#
#    FFT_THEOR = theoretical background spectrum as a function of
#                Fourier frequency. This will be smoothed by the
#                wavelet function and returned as a function of PERIOD.
#
#
# OPTIONAL KEYWORD OUTPUTS:
#
#    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
#           to the SCALEs.
#
#    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J
#            where J+1 is the total # of scales.
#
#    COI = if specified, then return the Cone-of-Influence, which is a vector
#        of N points that contains the maximum period of useful information
#        at that particular time.
#        Periods greater than this are subject to edge effects.
#        This can be used to plot COI lines on a contour plot by doing:
#            IDL>  CONTOUR,wavelet,time,period
#            IDL>  PLOTS,time,coi,NOCLIP=0
#
#    YPAD = returns the padded time series that was actually used in the
#         wavelet transform.
#
#    DAUGHTER = if initially set to 1, then return the daughter wavelets.
#         This is a complex array of the same size as WAVELET. At each scale
#         the daughter wavelet is located in the center of the array.
#
#    SIGNIF = output significance levels as a function of PERIOD
#
#    FFT_THEOR = output theoretical background spectrum (smoothed by the
#                wavelet function), as a function of PERIOD.
#
#----------------------------------------------------------------------------
# Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
#
# This software may be used, copied, or redistributed as long as it is not
# sold and this copyright notice is reproduced on each copy made.
# This routine is provided as is without any express or implied warranties
# whatsoever.
#
# Notice: Please acknowledge the use of the above software in any publications:
#    ``Wavelet software was provided by C. Torrence and G. Compo,
#      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
#
# Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
#            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
#
# Please send a copy of such publications to either C. Torrence or G. Compo:
#  Dr. Christopher Torrence               Dr. Gilbert P. Compo
#  Research Systems, Inc.                 Climate Diagnostics Center
#  4990 Pearl East Circle                 325 Broadway R/CDC1
#  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
#  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
#----------------------------------------------------------------------------

morlet <- function(k0,scale,k) 
{ 
  if (k0 == -1) k0 <- 6
  n = length(k)
  expnt = (-1*(scale*k - k0)^2)/2*(k > 0)
  dt = 2*pi/(n*k[2])
  norm = sqrt(2*pi*scale/dt)*(pi^(-0.25))           ## total energy=N   [Eqn(7)]
  morlet = norm*exp( pmax(expnt , -100) )
  morlet = morlet*(expnt > -100)                   ## avoid underflow errors
  morlet = morlet*(k > 0)                           ## Heaviside step function (Morlet is complex)
  fourier_factor = (4*pi)/(k0 + sqrt(2+k0^2)) ## Scale-->Fourier [Sec.3h]
  period = scale*fourier_factor
  coi = fourier_factor/sqrt(2)   ## Cone-of-influence [Sec.3g]
  dofmin = 2   ## Degrees of freedom with no smoothing
  Cdelta = -1
  if (k0 == 6) Cdelta <- 0.776 ## reconstruction factor
  psi0 = pi^(-0.25)
  ##	PRINT,scale,n,sqrt(TOTAL(ABS(morlet)^2,/DOUBLE))
  return(list(mother=morlet, param=k0, period=period, coi=coi, dofmin=dofmin, Cdelta=Cdelta, psi0=psi0))
}

paul <- function(m,scale,k)
{
  if (m == -1) m <- 4
  n = length(k)
  expnt = -1*(scale*k)*(k > 0)
  dt = 2*pi/(n*k[2])
  norm = sqrt(2*pi*scale/dt)*(2^m/sqrt(m*factorial(2*m-1)))
  paul = norm*((scale*k)^m)*exp(pmax(expnt ,-100))*(expnt > -100)
  paul = paul*(k > 0)
  fourier_factor = 4*pi/(2*m+1)
  period = scale*fourier_factor
  coi = fourier_factor*sqrt(2)
  dofmin = 2   ## Degrees of freedom with no smoothing
  Cdelta = -1
  if (m == 4) Cdelta = 1.132 ## reconstruction factor
  psi0 = 2.^m*factorial(m)/sqrt(pi*factorial(2*m))
  ##	PRINT,scale,n,norm,sqrt(TOTAL(paul^2,/DOUBLE))*sqrt(n)
  return(list(mother=paul, param=m, period=period, coi=coi, dofmin=dofmin, Cdelta=Cdelta, psi0=psi0))  
}

dog <- function(m,scale,k)
{
  if (m == -1) m = 2
  n = length(k)
  expnt = -(scale*k)^2/2
  dt = 2*pi/(n*k(1))
  norm = sqrt(2*pi*scale/dt)*sqrt(1/gamma(m+0.5))
  I = 0+1i
  gauss = -1*norm*(I^m)*(scale*k)^m*exp(pmax(expnt > -100))*(expnt > -100)
  fourier_factor = 2*pi*sqrt(2/(2*m+1))
  period = scale*fourier_factor
  coi = fourier_factor/sqrt(2)
  dofmin = 1   ## Degrees of freedom with no smoothing
  Cdelta = -1
  psi0 = -1
  if (m == 2) {
    Cdelta = 3.541 ## reconstruction factor
    psi0 = 0.867325
  }
  if (m == 6) {
    Cdelta = 1.966 ## reconstruction factor
    psi0 = 0.88406
  }
  ##	PRINT,scale,n,norm,sqrt(TOTAL(ABS(gauss)^2,/DOUBLE))*sqrt(n)
  return(list(mother=paul, param=m, period=period, coi=coi, dofmin=dofmin, Cdelta=Cdelta, psi0=psi0))
}
  
## ####################################################################
## class: waveletTransform
if (!isClassUnion("vectorOrLogical")) setClassUnion("vectorOrLogical", c("vector","logical"))
if (!isClassUnion("arrayOrLogical")) setClassUnion("arrayOrLogical", c("array","logical"))

setClass("waveletTransform",
         representation(transform="array",
                        s0="numeric",
                        dj="numeric",
                        j="numeric",
                        dt="numeric",
                        mother="character",
                        param="numeric",
                        Cdelta="numeric",
                        psi0="numeric",
                        pad="vector",
                        lag1="numeric",
                        siglvl="numeric",
                        recon="vectorOrLogical",
                        period="vector",
                        scale="vector",
                        coi="vector",
                        signif="vector",
                        fft.theor="vector",
                        daughter="arrayOrLogical",
                        dofmin="numeric",
                        input= "timeSeries")
         )

## spaceTime accessor functions
set.accessors( 'waveletTransform' )

######################################################################
## method: waveletTransform

## waveletTransform method, generic function
if (!isGeneric("waveletTransform")) {  ## creates a generic function
  if (is.function("waveletTransform"))
    fun <- subset else fun <- function(Y, ...) standardGeneric("waveletTransform")
  setGeneric("waveletTransform", fun)
}

setMethod("waveletTransform", c("timeSeries"), 
          function(Y, s0=2*dt, dj=.125, j=floor((log(ntime*dt/s0)/log(2))/dj),
                   mother='morlet', param=-1, pad=0,
                   lag1=0.0, siglvl=.95, recon=FALSE, fft_theor=NULL,
                   daughter=FALSE, verbose=TRUE, ...) {

            do.daughter <- daughter
            ntime <- length(a.POSIXct(Y))  ## ntime is the length of time, not padded
            ## this puts dt to the nearest fraction of a year.
            dt <- 1/round(1/mean( (as.numeric(a.POSIXct(Y)[2:ntime])-as.numeric(a.POSIXct(Y)[1:(ntime-1)]))/365/24/60/60 ))
            if (dt==0) warning('dt is zero, something wrong with input', immediate.=TRUE)
            ## check for uniformity of timeseries
            ## this is a bit tricky b/c of leap year. probably need to go with a fraction tolerance. 
            ## check that times steps are all uniform
            ##if (!all(as.numeric(a.POSIXct(Y)[2:ntime])-as.numeric(a.POSIXct(Y)[1:(ntime-1)]))/24/60/60==dt*))
            ##warning("waveletTransform: Timesteps are not uniform in input timeseries.", immediate.=TRUE)
          
            ##....construct time series to analyze, pad if necessary
            ypad <- a.data(Y) - mean(a.data(Y))

            if (pad >0){
              base2 <- trunc(log(ntime)/log(2) + 0.4999)   # power of 2 nearest to N
              ypad <- c(ypad, rep(0, 2^(base2 + 1) - ntime))
              if (pad >= 2) {
                pad.length <- min( ntime-2, length(ypad)-ntime-2 )
                ypad[(ntime+1):(ntime+1+pad.length)] <-
                  if (pad==2) rev(ypad[1:ntime])[2:(2+pad.length)] else
                              predict( ar( ypad[1:ntime] ), n.ahead=pad.length+1 )$pred
                
              }
              
            }
            n <- length(ypad)

            ##....construct SCALE array & empty PERIOD & WAVE arrays
            j1 <- j+1
            scale <- s0*2^((0:j)*dj)
            period <- scale*NA;
            wave <- (matrix(data=as.complex(0), nrow=j1, ncol=n))  # define the wavelet array
            if (do.daughter) daughter <- wave

            ##....construct wavenumber array used in transform [Eqn(5)]
            k <- (1:trunc(n/2)) * ((2*pi)/(n*dt))
            k <- c(0, k, -rev(k[1:floor((n-1)/2)]))
            
            ## ....compute FFT of the (padded) time series
            yfft = fft(ypad)    # [Eqn(3)]
            
            if (verbose) {
              print(mother)
              print(paste('#points=',n,'   s0=',s0,'   dj=',dj,'   J=',floor(j),sep=''))
              if (n != ntime) print(paste('(padded with ',n-ntime,' zeroes)',sep=''))
              print(paste('   j',' scale','   period  variance',sep='  '))  ## more to follow later...
            }

            if (length(fft_theor) == n) {
              fft_theor_k = fft_theor
            } else {
              fft_theor_k = (1-lag1^2)/(1-2*lag1*cos(k*dt)+lag1^2)  # [Eqn(16)]
              fft_theor = (1:j1) * 0 
            }

            ## loop through all scales and compute transform
            for(a1 in 1:j1){

              psi.fft <- do.call( mother, list(param, scale[a1], k))
              ## returns: $mother, $period. $coi, $dofmin, $Cdelta, $psi0
              wave[a1,] <- fft((yfft*psi.fft$mother), inverse = TRUE) / n  ## wavelet transform[Eqn(4)]
              period[a1] <- psi.fft$period
              fft_theor[a1] <- sum( (abs( psi.fft$mother)^2 ) * fft_theor_k ) / n
              if (do.daughter) daughter[a1,] <- fft( psi.fft$mother, inverse=TRUE)
              if (verbose) print(paste( format(a1,digits=1,nsmall=2,width=4),
                                        format(scale[a1],digits=1,nsmall=2, width=6),
                                        format(period[a1],digits=1,nsmall=2, width=8),
                                        sum(abs(wave[a1,])^2),
                                        sep='  '))
                
            }
            
            coi <- psi.fft$coi*c( (1:(floor(ntime + 1)/2))-1, rev( (1:floor(ntime/2))-1 )) * dt  ## COI [Sec.3g]

            if (do.daughter) daughter <- cbind( daughter[,(n-(ntime/2)+1):n], daughter[,1:(ntime/2)])
                        
            ## ....significance levels [Sec.4]
            fft_theor = var(a.data(Y))*fft_theor  # include time-series variance (it's var not sd as per the original var names)
            dof = psi.fft$dofmin
            signif = fft_theor*qchisq(1. - siglvl, dof, lower.tail=FALSE)/dof   ## [Eqn(18)]

            ## Reconstruction [Eqn(11)]
            if (recon){  
              if (psi.fft$Cdelta == -1) {
                recon = -1
                warning('Cdelta undefined, cannot reconstruct with this wavelet.')
              } else {
                recon=dj*sqrt(dt)/(psi.fft$Cdelta*psi.fft$psi0)*(t(Re(wave)) %*% (1./sqrt(scale)))
                recon = recon[1:ntime]
              }
            }

            return( new( 'waveletTransform', transform=wave[,1:ntime],
                        s0=s0, dj=dj, j=j,
                        dt=dt,
                        mother=mother, param=psi.fft$param,
                        Cdelta=psi.fft$Cdelta, psi0=psi.fft$psi0,
                        pad=ypad,
                        lag1=lag1, siglvl=siglvl,
                        recon=recon, fft.theor=fft_theor,
                        period=period, scale=scale, coi=coi,
                        signif=signif,  daughter=daughter, dofmin=dof,
                        input=Y
                        )
                 )
                 
               }
          )


