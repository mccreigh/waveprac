## name
## author: james mccreight = mccreigh ^ gmail * com
## info:
##   part of "A practical guide to wavelet analysis" for R,
##   adapted from Torrence & Compo .
##   http://paos.colorado.edu/research/wavelets/

## General Purpose:
## Get scale-averaged wavelet power from a wavelet transform.

## Return Value
## A numeric vector time series of scale averaged wavelet power.

## Inputs
##   wt - a wavelet transform object
##   low, high - the bounds on scales of averaging. Lower is
##               included but upper is not, e.g. [low, high).
##   unit - default scale is specified in period, can also
##          be specified in "scale"

scale.avg.wt <- function(wt, low, high, unit='period') {

  if (unit!='period' & unit!='scale')
    warning('scale.avg.wt: unit must be either of "period" or "scale".',immediate.=TRUE)

  ntime <- length(a.POSIXct(a.input(wt)));   nscale=length(a.scale(wt))
  power.norm <- abs(a.transform(wt))^2/matrix(a.scale(wt), ncol=ntime, nrow=nscale )
    
  wh.avg <- if (unit=="period") which((a.period(wt) > low) & (a.period(wt) <= high)) else 
                               which((a.scale(wt) > low) & (a.scale(wt) <= high))

  if (length(wh.avg)==0) warning('No scales in the specified range.', immediate.=TRUE)
  scale.avg <-
    if (length(wh.avg)>1) a.dj(wt)*a.dt(wt)/a.Cdelta(wt)*colSums(power.norm[wh.avg,])  else # [Eqn(24)]
                          a.dj(wt)*a.dt(wt)/a.Cdelta(wt)*power.norm[wh.avg,]
  return(scale.avg)
}
