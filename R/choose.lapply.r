
## This chooses between single and multi core versions of lapply (lapply and
## mclapply in the multicore package).

## Arguments can be passed to both via ... . The nprocessors replaces the mc.cores
## keyword to mclapply when nprocessors >1.  The defaults for mc.preschedule and mc.set.seed
## are set to true. See the multicore manual for further details:
## mc.preschedule is faster for shorter calls.
## mc.set.seed=TRUE is **essential**, otherwise the same seed is used on all cores and
## all simulations identical.

choose.lapply <- function( list, function.name, nprocessors=1,
                           mc.preschedule=TRUE, mc.set.seed=TRUE, ... ) {
  
  if (nprocessors >1) {
    have.multicore <- if ('package:multicore' %in% search()) TRUE else require(multicore)      
    if (!have.multicore) {
      print(paste("multicore package is NOT present, will run",
                  nsimulations,"simulations on 1 processor."))
      print("Do you wish to continue? (y/n)")
      if (tolower(substr(readLines(n=1),1,1))=='n') return(NULL) else nprocessors <- 1
    }
  }

  result <-
    if (nprocessors==1) 
      { lapply( list,  function.name, ... ) } else 
      { mclapply( list, function.name,
                  mc.cores=nprocessors, 
                  mc.set.seed=mc.set.seed, mc.preschedule=mc.preschedule,
                  ... )
      }

  invisible(result)
    
}
