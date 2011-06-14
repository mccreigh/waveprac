ts.subset <- function(ts, start, end) {

  ## temporal range subsetting
  wh.start <- which(format(a.POSIXct(ts),'%m/%d/%Y') == start)
  wh.end <- which(format(a.POSIXct(ts),'%m/%d/%Y') == end)
  if (length(wh.start)==0 | length(wh.end)==0)
    warning( 'Either start or end time not found in the timeSeries, check their values.', immediate.=TRUE)
  
  a.data(ts) <- a.data(ts)[wh.start:wh.end]
  a.POSIXct(ts) <- a.POSIXct(ts)[wh.start:wh.end]
  ts
  
}  
