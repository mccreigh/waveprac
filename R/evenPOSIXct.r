## name: evenPOSIXct
## purpose: to make an "evenly" spaced POSIXct timeseries with in each YEAR.
## Because years are not equal lengths of seconds, this is not trivial.

## ts - either:
##              a year.fraction timeseries *with potential gaps*
##              (e.g. c(1980.0, 1980.25, 1980.75, 1981.0) )
##      OR when used with origin and dY:
##              the length of an origin+dY generated timseries.
## origin, dY - must be supplied together
##      origin is a single year.fraction
##      dY - is the fraction of a year of the timestep.

evenPOSIXct <- function( nts , origin=NULL, dY=NULL, format='%d/%m/%Y', tz='UTC' ) {

  if (length(nts) > 1 ) warning("ts must be supplied as a scalar length of the timeseries.", .immediate=TRUE)
  if (missing(origin) | missing(dY)) warning('Both origin and dY required', immediate.=TRUE)

  ## The prupose here is to NOT build up round off error over more than one year. And to
  ## remove round off errors associated with the offset.
  origin.offset <- (origin-floor(origin))*round(1/dY)  ## how many dYs the first time is from the beginning of a year
  origin.offset.fract <- 1/round(1/(origin.offset - floor(origin.offset))) ## fractional dY, rounding here is key
  origin.y <- floor(origin) ## the base year 

  ## calculating from the beginning of each year prevents round off error.
  ## the base year + number of years to add + the number of dYs in a given year, including fractional offset.
  ts <- origin.y + ((origin.offset+(0:(nts-1))) %/% round(1/dY)) +
    (floor((origin.offset+(0:(nts-1))) %% round(1/dY)) + origin.offset.fract)/round(1/dY)

  ## now separate the whole years from the fractions of years
  y <- floor(ts) ; f <- ts-y  ## y.f

  ## now need to get the number of seconds in each year
  nsec=as.numeric(as.POSIXct(strptime(paste('01/01',y+1,sep='/'),format='%d/%m/%Y'),tz=tz)) -
        as.numeric(as.POSIXct(strptime(paste('01/01',y,sep='/'),format='%d/%m/%Y'),tz=tz))

  invisible(as.POSIXct(strptime(paste('01/01',y,sep='/'),format=format),tz=tz) + (nsec * f))
  
}

is.POSIXct <- function(x) any("POSIXct" %in% class(x))

## this calculates the temporal step size as a fraction of a year
## of a POSIXct variale or extracts one from a timeSeries object to do so on it.
get.dY <- function( ts ) {
  if (!is.POSIXct(ts)) ts <- a.POSIXct(ts) ## if it's not POSIXct then assume a timeseries was passed.
  ntime <-  length(ts)
  1/round(1/mean((as.numeric(ts[-1])-as.numeric(ts[-ntime]))/365/24/60/60 ))
}
