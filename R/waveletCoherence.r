##
##   WAVE_COHERENCE
##
## PURPOSE:   Compute the wavelet coherency between two time series.
##
##
## CALLING SEQUENCE:
##
##     WAVE_COHERENCY, $
##         wave1,time1,scale1, $
##         wave2,time2,scale2, $
##         WAVE_COHER=wave_coher,WAVE_PHASE=wave_phase, $
##         TIME_OUT=time_out,SCALE_OUT=scale_out
##
##
## INPUTS:
##
##    WAVE1 = wavelet power spectrum for time series #1
##    TIME1 = a vector of times for time series #1
##    SCALE1 = a vector of scales for time series #1
##    WAVE2 = wavelet power spectrum for time series #2
##    TIME2 = a vector of times for time series #2
##    SCALE2 = a vector of scales for time series #2
##
##
## OPTIONAL KEYWORD INPUTS:
##
##    DT = amount of time between each Y value, i.e. the sampling time.
##         If not input, then calculated from TIME1(1)-TIME1(0)
##
##    DJ = the spacing between discrete scales.
##         If not input, then calculated from SCALE1
##
##   VERBOSE = if set, then print out the scales and system time
##
##   NOSMOOTH = if set, then just compute the GLOBAL_COHER, GLOBAL_PHASE,
##              and the unsmoothed CROSS_WAVELET and return
##
##
## OPTIONAL KEYWORD OUTPUTS:
##
##   WAVE_COHER = the wavelet coherency, as a function of
##       TIME_OUT and SCALE_OUT
##
##   TIME_OUT = the time vector, given by the overlap of TIME1 and TIME2
##
##   SCALE_OUT = the scale vector of scale indices, given by the overlap
##               of SCALE1 and SCALE2
##
##   COI_OUT = the vector of the cone-of-influence
##
##	GLOBAL_COHER = the global (or mean) coherence averaged over all times.
##
##   GLOBAL_PHASE = the global (or mean) phase averaged over all times
##
##	CROSS_WAVELET = the cross wavelet between the time series
##
##   POWER1 = the wavelet power spectrum; should be the same as WAVE1
##            if TIME1 and TIME2 are identical, otherwise it is only the
##            overlapping portion. If NOSMOOTH is set,
##            then this is unsmoothed, otherwise it is smoothed.
##
##   POWER2 = same as POWER1 but for time series #2
##
##
##----------------------------------------------------------------------------
## Copyright (C) 1998-2005, Christopher Torrence
## This software may be used, copied, or redistributed as long as it is not
## sold and this copyright notice is reproduced on each copy made.  This
## routine is provided as is without any express or
## implied warranties whatsoever.
##
## Reference: Torrence, C. and P. J. Webster, 1999: Interdecadal changes in the
##            ENSO-monsoon system. <I>J. Climate</I>, 12, 2679-2690.
##
## Please send a copy of any publications to C. Torrence:
##  Dr. Christopher Torrence
##  Research Systems, Inc.
##  4990 Pearl East Circle
##  Boulder, CO 80301, USA
##  E-mail: chris[AT]rsinc[DOT]com
##----------------------------------------------------------------------------
##-

## R version 4/5/11 - james mccreight - mccreigh 2 gmail > com

setClass("waveletCoherence",
         representation(coherence="array",
                        phase="array",
                        POSIXct="POSIXct",
                        scale="vector",
                        dt='numeric', dj='numeric', 
                        coi="vector",
                        global_coherence="vector",
                        global_phase="vector",
                        cross.wavelet="array",
                        wt1.name="character",
                        wt1.units="character",
                        wt2.name="character",
                        wt2.units="character",
                        nosmooth="logical",
                        monteCarlo="character",
                        confidence="numericNULL",
                        nsimulations="numericNULL",
                        qmc.coherence="arrayNULL",
                        qmc.global_coherence="vectorNULL"
                        )
         )
set.accessors( 'waveletCoherence' )

##****************************************************************** WAVELET
## waveletTransform method, generic function
if (!isGeneric("waveletCoherence")) {  ## creates a generic function
  if (is.function("waveletCoherence"))
    fun <- subset else fun <- function(wt1, wt2, ...) standardGeneric("waveletCoherence")
  setGeneric("waveletCoherence", fun)
}

setMethod("waveletCoherence", c("waveletTransform", "waveletTransform"), 
          function( wt1, wt2, verbose=FALSE, nosmooth=FALSE,
                   monteCarlo='white', confidence=.95, nsimulations=500, nprocessors=1 )
{

  ## Define a coherence function which can be called for both the timeseries coherence and the
  ## montecarlo confidence simulations.           
  

  coherence <- function( wt1, wt2, verbose=FALSE, nosmooth=FALSE ) {
    ##  wave1,time1,scale1,wave2,time2,scale2, $   ##*** required inputs
    ##  COI1=coi1, $
    ##  DT=dt,DJ=dj, $
    ##  WAVE_COHER=wave_coher,WAVE_PHASE=wave_phase, $
    ##  TIME_OUT=time_out,SCALE_OUT=scale_out,COI_OUT=coi_out, $
    ##  GLOBAL_COHER=global_coher,GLOBAL_PHASE=global_phase, $
    ##  CROSS_WAVELET=cross_wavelet,POWER1=power1,POWER2=power2, $
    ##  NOSMOOTH=nosmooth, $
    ##  VERBOSE=verbose
    time1=a.POSIXct(a.input(wt1))
    time2=a.POSIXct(a.input(wt2))
    scale1=a.scale(wt1)
    scale2=a.scale(wt2)
    
    ## find overlapping times
    timeStart  = max(min(time1),min(time2))
    timeEnd    = min(max(time1), max(time2))
    time1Start = min(which((time1 >= timeStart)))
    time1End   = max(which((time1 <= timeEnd)))
    time2Start = min(which((time2 >= timeStart)))
    time2End   = max(which((time2 <= timeEnd)))
    
    ## find overlapping scales
    scaleStart  = max(min(scale1),min(scale2))
    scaleEnd    = min(max(scale1),max(scale2))
    scale1Start = min(which((scale1 >= scaleStart)))
    scale1End   = max(which((scale1 <= scaleEnd)))
    scale2Start = min(which((scale2 >= scaleStart)))
    scale2End   = max(which((scale2 <= scaleEnd)))
    
    ##*** cross wavelet & individual wavelet power
    cross_wavelet <-
      a.transform(wt1)[scale1Start:scale1End,time1Start:time1End] *
        Conj(a.transform(wt2)[scale2Start:scale2End,time2Start:time2End])
    power1 = Mod(a.transform(wt1)[scale1Start:scale1End,time1Start:time1End])^2
    power2 = Mod(a.transform(wt2)[scale2Start:scale2End,time2Start:time2End])^2
    
    dt=a.dt(wt1)
    ntime = time1End - time1Start + 1
    nj = scale1End - scale1Start + 1
    dj = a.dj(wt1) 
    scale = scale1[scale1Start:scale1End]
    if (verbose) print(paste(dt,ntime,dj,nj,sep=' '))
    time_out = time1[time1Start:time1End]
    scale_out = scale1[scale1Start:scale1End]
    coi_out = a.coi(wt1)[time1Start:time1End]
    
    ## calculate global coherency before doing local smoothing
    global1 = rowSums(power1)
    global2 = rowSums(power2)
    global_cross = rowSums(cross_wavelet)
    global_coher = Mod(global_cross)^2/(global1*global2)
    global_phase = 180/pi*atan2(Im(global_cross),Re(global_cross))
    
    ## time-smoothing
    for ( j in 1:nj ) {
      ## st1 = SYSTIME(1)
      nt = floor(4*scale[j]/dt)/2*4 + 1
      
      time_wavelet = ( (0:(nt-1)) - floor(nt/2) ) * dt/scale[j] 
      wave_function = exp(-time_wavelet^2/2)   ##*** Morlet
      wave_function = Re(wave_function/sum(wave_function)) ## normalize
      nz = floor(nt/2)
      
      zeros = array(as.complex(0), nz )
      cross_wave_slice = c(zeros,cross_wavelet[j,],zeros)
      cross_wave_slice = convolve(cross_wave_slice,wave_function, type='filter')
      cross_wavelet[j,] = cross_wave_slice ##[(nz+1):(ntime+nz)]
      
      zeros = Re(zeros)
      power_slice = c(zeros,power1[j,],zeros)
      power_slice = convolve(power_slice,wave_function, type='filter')
      power1[j,] = power_slice #[(nz+1):(ntime + nz)]
      power_slice = c(zeros,power2[j,],zeros)
      power_slice = convolve(power_slice,wave_function, type='filter')
      power2[j,] = power_slice #[(nz+1):(ntime + nz)]
      if (verbose) print(paste(j,scale[j],sep=' '))  #,SYSTIME(1)-st1##,FORMAT='(I4,$)'
    }
    
    ## normalize by scale
    scales = matrix(scale, ncol=ntime, nrow=nj)
    cross_wavelet = cross_wavelet/scales
    power1 = power1/scales
    power2 = power2/scales
    
    nweights = floor(0.6/dj/2 + 0.5)*2 - 1   ## closest (smaller) odd integer
    weights = rep(1.,nweights)
    weights = weights/sum(weights) # normalize
    
    ##*** scale-smoothing
    if (nweights>1) {  ## point less for nweights=1
      pad=rep(0,floor(nweights/2))
      for (i in 1:ntime) {
        cross_wavelet[,i] = c(pad, convolve(cross_wavelet[,i],weights, type='filter'), pad)
        power1[,i] = c(pad, convolve(power1[,i],weights, type='filter'), pad)
        power2[,i] = c(pad, convolve(power2[,i],weights, type='filter'), pad)
      }
    }
    
    wave_phase = 180./pi*atan2(Im(cross_wavelet),Re(cross_wavelet))
    wave_coher = (abs(cross_wavelet)^2)/pmax(power1*power2, 1E-9)
    
    
    return( new("waveletCoherence",
                coherence=wave_coher,
                phase=wave_phase,
                POSIXct=time_out,
                scale=scale_out,
                dt=dt, dj=dj, 
                coi=coi_out,
                global_coherence=global_coher,
                global_phase=global_phase,
                cross.wavelet=cross_wavelet,
                wt1.name=a.data.name(a.input(wt1)),
                wt1.units=a.data.units(a.input(wt1)),
                wt2.name=a.data.name(a.input(wt2)),
                wt2.units=a.data.units(a.input(wt2)),
                nosmooth=nosmooth,
                monteCarlo=monteCarlo,  ## these get set before the return of waveletCoherence
                confidence=NULL,
                nsimulations=NULL,
                qmc.coherence=NULL,
                qmc.global_coherence=NULL
                )
           )
    
  }


  ## main routine
  if (!(monteCarlo=='white' | monteCarlo=='red' | monteCarlo=='spectrum')) monteCarlo <- 'none'
  
  ## coherence of input timeseries
  waveletCoherence.out <- coherence( wt1, wt2, verbose=verbose, nosmooth=nosmooth )

  ## confidence simulations
  if (monteCarlo != 'none')   {
    
    if (nprocessors >1) {
      have.multicore <- require(multicore)      
      if (!have.multicore) {
        print(paste("multicore package is NOT present, will run",nsimulations,"simulations on 1 processor."))
        print("Do you wish to continue? (y/n)")
        if (tolower(substr(readLines(n=1),1,1))=='n') return(NULL) else nprocessors <- 1
      }
    }    
    
    ## pre-compute the statistics for generating noise timeseries
    ## this could probably be sped up by precomputing the noise timeseries... could be indexed by dum in the call... maybe.
    pre.noise.1 <- pre.noise( a.data(a.input(wt1)), type=monteCarlo )
    pre.noise.2 <- pre.noise( a.data(a.input(wt2)), type=monteCarlo )
    
    sim.loop <- function( dum ) {
      ## could use dum to give verbose progress diagnostics?? maybe?
      
      ## generate 2 pseudo-random timeseries replacing it in the waveletTransforms passed to the calling level
      a.data(a.input(wt1)) <- noise( pre.noise.1 ) 
      a.data(a.input(wt2)) <- noise( pre.noise.2 )
      
      ## wavelet transform the timeseries' (wt1 and wt2 are in the calling level)
      wt1.mc <- waveletTransform(a.input(wt1), s0=a.s0(wt1), dj=a.dj(wt1), j=a.j(wt1),
                                 mother=a.mother(wt1), param=a.param(wt1), pad=length(a.pad(wt1))!=length(a.data(a.input(wt1))),
                                 lag1=a.lag1(wt1), siglvl=a.siglvl(wt1), verbose=verbose)
      
      wt2.mc <- waveletTransform(a.input(wt2), s0=a.s0(wt2), dj=a.dj(wt2), j=a.j(wt2),
                                 mother=a.mother(wt2), param=a.param(wt2), pad=length(a.pad(wt2))!=length(a.data(a.input(wt2))),
                                 lag1=a.lag1(wt2), siglvl=a.siglvl(wt2), verbose=verbose)
      
      ## calculate wavelet coherence, setting monteCarlo='none' avoids infinte recursion.
      wc.mc <- waveletCoherence( wt1.mc, wt2.mc, monteCarlo='none' )
      
      ## local and global coherence
      mc.coher <- a.coherence(wc.mc)
      mc.global_coher <- a.global_coherence(wc.mc)
      
      list( coherence=mc.coher, global_coherence=mc.global_coher )      
    }
    
    
    sim.list <- if (nprocessors == 1) lapply( 1:nsimulations, sim.loop) else {
      print("Killing this process may result in zombie-like children")
      ## mc.preschedule is faster for shorter calls. mc.set.seed=TRUE is **essential**, otherwise the same seed is used on all cores.
      mclapply( 1:nsimulations, sim.loop, mc.set.seed=TRUE, mc.cores=nprocessors, mc.preschedule=TRUE )
    }

    ## calculate the quantiles in a distributed fashion if multicore is in use, significant speed up.
    coherence.quantiles <- if (nprocessors == 1) {
      lapply( as.list( data.frame( t( sapply( sim.list, "[[", 'coherence') ) ) ), quantile, probs=confidence) 
    } else {
      mclapply( as.list( data.frame( t( sapply( sim.list, "[[", 'coherence') ) ) ), quantile, probs=confidence,
                mc.cores=nprocessors, mc.preschedule=TRUE)
    }
      
    global_coherence.quantiles <- if (nprocessors == 1) {
      lapply( as.list( data.frame( t( sapply( sim.list, "[[", 'global_coherence') ) ) ), quantile, probs=confidence )
    } else {
      mclapply( as.list( data.frame( t( sapply( sim.list, "[[", 'global_coherence') ) ) ), quantile, probs=confidence,
               mc.cores=nprocessors, mc.preschedule=TRUE)
    }

    ## reshape to the wavelet transform dims
    coherence.quantiles <- array(unlist(coherence.quantiles), dim=dim(a.coherence(waveletCoherence.out)) )
    global_coherence.quantiles <- array(unlist(global_coherence.quantiles), dim=length(a.global_coherence(waveletCoherence.out)) )

    ## set the output
    a.confidence(waveletCoherence.out) <- confidence
    a.nsimulations(waveletCoherence.out) <- nsimulations
    a.qmc.coherence(waveletCoherence.out) <- coherence.quantiles
    a.qmc.global_coherence(waveletCoherence.out) <- global_coherence.quantiles

  } 

  invisible(waveletCoherence.out)
}  
)


## noise functions. need amended for red noise and perhaps phase shuffling to preserve the global spectrum?
pre.noise <- function( ts, type='white', ... )
{
  length=length(ts)
  mean=mean(ts)
  sd=sd(ts)  
  if (type=='white') return( list(type=type, length=length, mean=mean, sd=sd) )
}

noise <- function( pre.list ){
  if (pre.list$type=='white') noise <- ( (rnorm(pre.list$length))*pre.list$sd + pre.list$mean )
  return(noise)
}

