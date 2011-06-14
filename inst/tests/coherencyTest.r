src("pgwa.r")

## nino sst
sst <- read.table('sst_nino3.dat')[,1]
sst <- (sst - sum(sst)/length(sst))
dt <- 0.25
POSIXct <- evenPOSIXct( length(sst), origin=1871.0, dY=.25)
## define a timeseries object
nino <- new( 'timeSeries', data = sst,
            data.name = "nino3 sst index",
            data.units = "deg.C",
            POSIXct = POSIXct)
nino.wt <- waveletTransform(nino, s0=.25, pad=TRUE, dj=.25, j=9./.25, mother="morlet", recon=TRUE, daughter=TRUE, verbose=FALSE)

## monsoon india rainfall
rain <- read.table('monsoon.dat')[,1]
rain <- (rain - sum(rain)/length(rain))
dt <- 0.25
POSIXct <- evenPOSIXct( length(rain), origin=1871.0, dY=.25)
## define a timeseries object
rain <- new( 'timeSeries', data = rain,
            data.name = "monsoon rainfall",
            data.units = "mm",
            POSIXct = POSIXct)
rain.wt <- waveletTransform(rain, s0=.25, pad=TRUE, dj=.25, j=9./.25, mother="morlet", recon=TRUE, daughter=TRUE, verbose=FALSE)

## coherency
start
nino.rain.coherence <- waveletCoherence( nino.wt, rain.wt, verbose=FALSE, monteCarlo='white', nsim=120, nproc=12) 

ggCoh <- ggCoherence( nino.rain.coherence, coi=TRUE )
