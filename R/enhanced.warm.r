## ... is for options to ar() and predict.ar()
## could also accept a function argument
## must differntiate between simulation and prediction... by n.ahead==0???
enhanced.warm <- function( wt, windows, n.ahead=0, nprocessors=1, nsims=1, ... )
{

  wt.filter <- waveletFilter(wt, windows)

  nwindows= dim( as.matrix(wt.filter$filter) )[1]
  band.list <- as.list(1:nwindows)
  names(band.list) <- wt.filter$windows.labels.out
  
  warm.band <- function( band ) {
    sawp <- scale.avg.wt( wt, wt.filter$windows.in[band,1], wt.filter$windows.in[band,2] ) 
    band.norm <- wt.filter$filter[band,]/sawp
    fit.ar <- ar( band.norm )

    ## simulate
    if (n.ahead == 0) {
      do.sim <-  function( dum ) {
        samp.start <- sample( length(band.norm)-fit.ar$order ,1 )
        sim <- band.norm[samp.start:(samp.start+fit.ar$order)]
        while (length(sim) < length(band.norm))
          sim <-c( sim, predict(fit.ar, tail(sim,fit.ar$order), n.ahead=1 )$pred )
        sim * sawp
      }
      out <- choose.lapply( as.list(1:nsims), do.sim, nprocessors=nprocessors, ... )
    }
    
    ## predict
    if (n.ahead > 0) {
      pred.sim <- predict( fit.ar, n.ahead=n.ahead, ... )
      
      ##  model the FUTURE SAWP with past SAWP, though not totally wild
      ##band.pred <- as.numeric( pred.sim$pred * rev(sawp)[2:(n.ahead+1)] )
      ntime <- length(sawp)
      future.sawp <- predict( ar( sawp[(ntime-ar(sawp)$order):ntime] ), n.ahead=n.ahead )
      band.pred <- as.numeric( pred.sim$pred ) * as.numeric( future.sawp$pred )
      
      ## additive error model recursive function
      ##band.se <- as.numeric( pred.sim$se * rev(sawp)[2:(n.ahead+1)] )
      band.se <- as.numeric( pred.sim$se ) * as.numeric( future.sawp$pred )
      out <- list(pred=band.pred, se=band.se)
    }

    out
  }

  ## apply over all windows
  warm.bands <- llply( band.list, warm.band )    

  if (n.ahead == 0) {
    matricize <- function( list ) matrix( unlist(list), byrow=T, ncol=length( list[[1]] ) )
    sims <- laply( warm.bands, matricize )
    sims <- aaply( dum, c(2,3), sum)
  }
  
  if (n.ahead > 0) {
    pred <- colSums(laply( warm.bands , '[[', 'pred' )) 
    se.band <- laply( warm.bands, '[[', 'se')
    additive.error <-  ## could probably replace this via Reduce()
      function( vec ) if (length(vec)>1) sqrt(vec[1]^2 + additive.error(vec[2:length(vec)])^2) else vec[1]
    se <- aaply( se.band, 2, additive.error)    
    pred <- list(pred=pred, se=se)
  }

  if (n.ahead == 0) sim else pred
}
  
