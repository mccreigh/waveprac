require(devtools)
load_all("waveprac")
load_data('waveprac')

#---------------------------------------------------------------------
## lees ferry example
## could debate about where to put these with in each year. should probably be 10/1, or 12/31 or even 9/30/year+1
lf <- new( 'timeSeries', data = lees.ferry$annual.acreft,
          data.name = "Lee's Ferry Annual Flows",
          data.units = "acre*ft",
          POSIXct = evenPOSIXct( length(lees.ferry$year), origin=lees.ferry$year[1]+1/2, dY=1) )
lf.wt <- waveletTransform(lf, pad=TRUE, dj=.25, mother="morlet")

lf.wt.check <- check.wt.variance(lf.wt)
gg.list <- ggWavelet( lf.wt, siglvl=.90, siglvl.global=.90, siglvl.scale=.90,
                      gws=TRUE, coi=TRUE, plot=TRUE )
gg.list <- ggWavelet( lf.wt, siglvl=.90, siglvl.global=.95, siglvl.scale=.90,
                     scale.avg=c(9,16), gws=TRUE, coi=TRUE, plot=TRUE, use.current=T )

##still working on these
ewarm.sim <- enhanced.warm( lf.wt, windows=matrix( c(9,16, 50,70), byrow=T, nrow=2 ), n.ahead=0, nsim=70, npr=30 )
ewarm.pred <- enhanced.warm( lf.wt, windows=matrix( c(9,16, 50,70), byrow=T, nrow=2 ), n.ahead=5, nsim=70, npr=30 )

#---------------------------------------------------------------------
## nino3 example

nino3 <- new( 'timeSeries', data = nino3$nino3,
          data.name = "Nino 3 SST index",
          data.units = "deg C",
          POSIXct = evenPOSIXct( length(nino3$year), origin=nino3$year[1]+1/8, dY=1/4) )

nino3.wt <- waveletTransform(nino, pad=TRUE, dj=.25, s0=.25, j=9/.25 , mother="morlet")
nino3.wt.check <- check.wt.variance(nino3.wt)

gg.nino3.wt <- ggWavelet( nino3.wt, siglvl=.90, siglvl.global=.90, siglvl.scale=.90,
                         gws=TRUE, coi=TRUE, plot=TRUE )
gg.nino3.wt <- ggWavelet( nino3.wt, siglvl=.90, siglvl.global=.90, siglvl.scale=.90,
                         scale.avg=c(2,8), gws=TRUE, coi=TRUE, plot=TRUE )
ggFilter( nino3.wt )
nino3.filter <- waveletFilter( nino3.wt, window=c(2,8,45) )
gg.nino3.filter <- ggFilter( nino3.wt, nino3.filter )

## pull some data from the ggplot objects
orig <- gg.nino3.filter$orig.recon.ts$data
orig <- orig[which(orig$variable=='original'),]
orig$variable <- as.character(orig$variable)

low.freq <- gg.nino3.filter$filter.ts$data
low.freq <- low.freq[which(as.numeric(low.freq$variable)==2),]
low.freq$variable <- as.character(low.freq$variable)

orig.low <- rbind(orig, low.freq)
orig.low$variable <- factor(orig.low$variable)
ggplot( orig.low, aes(x=POSIXct, y=value, color=variable) ) + geom_line()
