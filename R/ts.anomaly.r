## a function to coerce a timeSeries object to a data.frame
ts.as.df <- function( ts ) {
  data <- a.data(ts)
  attributes(data) <- list( data.name=a.data.name(ts), data.units=a.data.units(ts) )
  invisible( data.frame( POSIXct=a.POSIXct(ts), data=data ) )
}

## could probably be more succinctly written in plyr though would probably loose multicore lapply....

## accepts both timeSeries objs (univariate) or unnivariate and multivariate data frames
## with a single (common in the case of multivariate) POSIXct variable.
ts.anomaly <- function(ts.in, format=NULL, period=NULL, ...)  {

  ## if class is timeSeries, convert to a data.frame. covert back at the end
  ts <- if ('timeSeries' %in% class(ts.in)) ts.as.df(ts.in) else ts.in

  ## identify the POSIXct variable and others.
  ts.names <- names(ts)
  wh.POSIXct <-  which( ts.names=='POSIXct' )
  wh.data <- which( ts.names != 'POSIXct' )
  if (length(wh.POSIXct)<1)
    warning("There does not appear to be a POSIXct variable in the input timeSeries.")
  if (length(wh.data)<1) warning("There does not appear to be data in the input timeSeries.")

  ## can calculate this outside the following function since all data on the same POSIXct.
  levs=as.numeric(format(ts$POSIXct, format=format))
  ## function to actually calculate the anomalies.
  calc.anoms <- function( ts.list.sub, ... ) {
    lev.means <- aggregate( ts.list.sub, list( levs ), mean, ... )
    for (ll in  unique(lev.means$Group.1) ) {
      ts.list.sub[which(levs==ll)] <-
        ts.list.sub[which(levs==ll)] - lev.means$x[which(lev.means$Group.1==ll)]
    }
    ts.list.sub
  }

  ts.anom.out <- as.data.frame( choose.lapply( as.list( ts[-wh.POSIXct] ), calc.anoms, ...) )

  ts.anom.out$POSIXct <-
    if ('timeSeries' %in% class(ts.in)) ts.in@POSIXct else ts.in$POSIXct ## put it back on there.
  
  ## if it was a timeSeries, give output as such.
  if ('timeSeries' %in% class(ts.in))
    { a.data(ts.in) <- ts.anom.out[,which(colnames(ts.anom.out)!='POSIXct')]; ts.anom.out <- ts.in }

  invisible(ts.anom.out)
}



