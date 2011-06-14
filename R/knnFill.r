## fill in timeseries in a spaceTime object
## reguire: knnflex package
## I'm probably not handling cases where all observations are missing.

knnFill <- function( st ) {

  require(knnflex)

  ntime <- length(a.POSIXct(st)); space <- length(a.lon(st))

  wh.all.good <- which( !is.na( colSums( a.data(st) ) ) )
  
  ## same training set for all fixes
  train <- t(a.data(st)[,wh.all.good])
  
  ## fill the individual timeseries
  for (ss in 1:nspace ) {
    
    wh.miss <- which(is.na(a.data(st)[ss,]))     
    if (length(wh.miss)==0) next
    test <- t(a.data(st)[,wh.miss])
    
    response.all <- a.data(st)[ss,wh.all.good]  ## a tad redundant, since this is in the training set as well.
    
    kdist <- knn.dist(rbind(train,test))
    
    ltrain=length(train[,1]); ltot=ltrain+length(test[,1])
    pred <- knn.predict(1:ltrain, (ltrain+1):ltot,
                        response.all,
                        kdist, k=25, agg.meth='mean')

    a.data(st)[ss,wh.miss] <- pred
   
    ##dum <- data.frame(x=a.POSIXct(st), y=a.data(st)[ss,])    
    ##print(ggplot( dum , aes(x=x,y=y) ) + geom_point() +
    ##      geom_point( data=dum[whObsMiss,], colour='red' ) +opts(title=theObs) )
    
  }

  return(st)
  
}

  
