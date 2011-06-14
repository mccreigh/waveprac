## name: filterField (though this should just be waveletFilter method on a spaceTime object)
## general prupose:
## take a spaceTime object.
## return a SpaceTimeBands object.

## inputs:
## st  - a spaceTime object.
## wt - a wavelet transform object.
## filter - a filter on the above wavelet transform.
## (this setup assumes you'll have a predictand timeseries in wt who's
## global wavelet power spectrumm you have looked at and
## defined a filter upon - thus you just transfer the info from that
## analysis to the spaceTime object.)
## array - lets you get back an array of dims nband,nspace,ntime
##         rather than a list.

######################################################################

## filterField method, generic function
if (!isGeneric("filterField")) {  ## creates a generic function
  if (is.function("filterField"))
    fun <- subset else fun <- function(st, wt=wt, ...) standardGeneric("filterField")
  setGeneric("filterField", fun)
}

#####################################################################

setMethod("filterField", c("spaceTime","waveletTransform"),
function( st, wt=wt, filter=filter, array=FALSE, nprocessors=1, ... )  {
  
  ## ingest the wt, take the necessary params, then rm(wt)
  s0=a.s0(wt);  dj=a.dj(wt);  j=a.j(wt)
  dt=a.dt(wt) ## check dt later
  rm(wt)  ## discard all the excess we dont need
  
  ## ingest the filter, take window.out, rm(filter)
  window.out=filter$window.out
  window.labels.out=filter$window.labels.out
  rm(filter) ## discard the excess

  nspace=length(a.lon(st)); ntime=length(a.POSIXct(st)); nbands=length(window.labels.out)
  
  ## "loop" over spatial points
  filter.ts <- function( ts ) {
    ##   check for completness of each timeseries
    if (any(is.na(ts))) return( ts*NA )

    ## have to make it a timeSeries object
    ts <- new( 'timeSeries', data = ts, data.name = "", data.units = "",
               POSIXct = a.POSIXct(st) )

    ##   wavelet transform
    ts.wt <- waveletTransform( ts, s0=s0, dj=dj, j=j, verbose=FALSE )
    if (a.dt(ts.wt) != dt) warning('Timesteps dont match')  ## maybe do this some where else?
    ##   filter timeseries
    d <- waveletFilter( ts.wt, window=window.out )$filter
  }

  ## the syntax gets a little funky, but it's worth the potential speed up.
  ff <- choose.lapply( as.list( as.data.frame(t(a.data(st)))), filter.ts, nprocessors=nprocessors,... )
  ffl <- choose.lapply( as.list(1:nbands),
                        function(b) matrix(t(as.matrix(as.data.frame( lapply( ff, function(x) x[b,])))),
                                           nrow=nspace),
                        nprocessors=nprocessors, ...)
  names(ffl) <- window.labels.out
  
  new("spaceTimeBands",
      data = ffl,
      data.name = a.data.name(st),
      data.units = a.data.units(st) ,
      lon = a.lon(st), 
      lat = a.lat(st), 
      POSIXct = a.POSIXct(st),
      bands= window.out,
      s0=s0,
      j=j,
      dj=dj
      )
  
}
          
)
