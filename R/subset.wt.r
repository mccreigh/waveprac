if (!isGeneric("subset")) {  ## creates a generic function
  if (is.function("subset"))
    fun <- subset else fun <- function(x, ...) standardGeneric("subset")
  setGeneric("subset", fun)
}

setMethod("subset", "waveletTransform", 
          function(x, start, end)
 {
   wt <- x; rm(x)
   inds <- if (missing(end)) start else
              which( a.POSIXct(a.input(wt))==start ) : which( a.POSIXct(a.input(wt))==end )

   a.transform(wt) <- a.transform(wt)[,inds]
   a.pad(wt) <- a.pad(wt)[inds]
   a.coi(wt) <- a.coi(wt)[inds]
   a.data(a.input(wt)) <- a.data(a.input(wt))[inds]
   a.POSIXct(a.input(wt)) <- a.POSIXct(a.input(wt))[inds]
   invisible(wt)
 })
          
