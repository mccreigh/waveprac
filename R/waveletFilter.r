## name: waveletFilter
## author: james mccreight = mccreigh ^ gmail * com
## info:
##   part of "A practical guide to wavelet analysis" for R,
##   adapted from Torrence & Compo .
##   http://paos.colorado.edu/research/wavelets/

## General Purpose:
## Get scale-averaged wavelet power from a wavelet transform.

## Return Value
##   filter - A numeric vector time series of scale averaged wavelet power.
##   windows.in - the windows passedin
##   windows.labels.in - convenient labels for the above
##   windows.out - the actual windows (inclusive) depending on
##                the discrete wt scale/period
##   window.labels.out - labels
##   unit - was "scale" or "period" used?

## Inputs
##   wt - a wavelet transform object
##   windows - either a vector with the bounds/breaks on scales of
##            filtering or a matrix of such vectors which allows
##            gaps.
##            Note that the lower bound is not included in each window
##            while the upper is. Though, you almost have to know
##            the discrete periods exactly for this to matter much.
##            The return slot windows.out gives the inclusive bounds
##            actually used.
##            EXS: windows=c(a,b,c) gives filtered timeseries on the
##                                 intervals (a,b], (b,c].
##                 windows= a b
##                         c d   gives filtering on (a,b] and (c,d].
##                 
##  unit="period" - or "scale", compares windows against these in the wt.
##  plot=FALSE - if true calls ggFilter
##  ... options to ggFilter plotting routine.

waveletFilter <- function(wt, windows=NULL, unit='period', plot=FALSE, ...) {

  if (unit!='period' & unit!='scale')
    warning('scale.avg.wt: unit must be either of "period" or "scale".',immediate.=TRUE)

  ntime <- length(a.POSIXct(a.input(wt)));   nscale=length(a.scale(wt))
  
  ## construct the normalized real part of the WT at all scales
  real.norm <- Re(a.transform(wt)) / sqrt(matrix(a.scale(wt), ncol=ntime, nrow=nscale))

  if (!is.matrix(windows) & !is.vector(windows)) warning('windows must either be a matrix or a vector.', immediate.=TRUE)

  ## construct a list from breaks
  if (!is.matrix(windows)) {
    windows.in <- matrix(NA, ncol=2, nrow=length(windows)-1)
    for (ww in 1:(length(windows)-1)) { windows.in[ww,] <- c(windows[ww],windows[ww+1]) }
  } else windows.in <- windows
  windows.out <- windows.in
  
  nwindows <- dim(windows.in)[1]
  filtered <- matrix(NA, nrow=nwindows, ncol=ntime)
  scale.dum <- if (unit=="period") a.period(wt) else a.scale(wt)
  for (ww in 1:nwindows) {
    wh.filter <- if (ww==1)  ##include lowest if first window.
                            which( scale.dum >= windows.in[ww,1] & scale.dum <= windows.in[ww,2] ) else
                            which( scale.dum > windows.in[ww,1] & scale.dum <= windows.in[ww,2] ) 
    if (length(wh.filter)==0) { filtered=filtered[1:(ww-1),]; break }
    windows.out[ww,] <- c( min(scale.dum[wh.filter]), max(scale.dum[wh.filter]) )
    ## make sure it's a matrix or colSums will complain.
    filtered[ww,] <- ( (a.dj(wt)*a.dt(wt)^.5)/(a.Cdelta(wt)*a.psi0(wt)) ) *
      colSums( matrix(real.norm[wh.filter,],nrow=length(wh.filter)) )  ## [(Eqn 29)]
  }


  fmt.cut <- function(x) gsub(" ","",paste('[', format(x[1],nsmall=2,digits=1),',',format(x[2],nsmall=2,digits=1),']',sep=''))
  windows.labels.in <-  apply( windows.in, 1, fmt.cut)
  windows.labels.out <- apply( windows.out,1, fmt.cut)

  out <- list( filter=filtered,
               windows.in= windows.in,
               windows.labels.in=windows.labels.in,
               windows.out= windows.out,
               windows.labels.out=windows.labels.out,
               unit=unit
              )

  if (plot) ggFilter( wt, out, ... )

  invisible(out)
  
}


