## name: wavetest
## author: james mccreight = mccreigh ^ gmail * com
## info:
##   adapted from Torrence and Compo: A Practical Guide to Wavelet
##   analysis. See http://paos.colorado.edu/research/wavelets/
## general purpose: Example R program using NINO3 SST dataset for wavelet analysis.

src("pgwa.r")

sst <- read.table('sst_nino3.dat')[,1]

## normalize by standard deviation (not necessary, but makes it easier
## to compare with plot on Interactive Wavelet page, at
## "http://paos.colorado.edu/research/wavelets/plot/"
sst <- (sst - sum(sst)/length(sst))
dt <- 0.25
POSIXct <- evenPOSIXct( length(sst), origin=1871.0, dY=.25)
## alternatively
## year.frac <- (1:(length(sst))-1)*dt + 1871.0  # construct time array
## POSIXct <- evenPOSIXct( year.frac )

## define a timeseries object
nino <- new( 'timeSeries', data = sst,
            data.name = "nino3 sst index",
            data.units = "deg.C",
            POSIXct = POSIXct)

## Wavelet transform on timeseries object:
nino.wt <- waveletTransform(nino, s0=.25, pad=TRUE, dj=.25, j=9./.25, mother="morlet", recon=TRUE, daughter=TRUE)

## Parseval and reconstruction check.
nino.wt.check <- check.wt.variance(nino.wt)

## Estimate lag-1 autocorrelation, for red-noise significance tests
## Note that we actually use the global wavelet spectrum (GWS)
## for the significance tests, but if you wanted to use red noise,
## here's how you could calculate it... though this isnt used for anything.... 
sst.acf<-acf(sst,lag.max=2,plot=FALSE)$acf
lag1 <- (sst.acf[2]+sqrt(sst.acf[3]))/2

gg.list <- ggWavelet( nino.wt, siglvl=.90, siglvl.global=.95, siglvl.scale=.90,
                      gws=TRUE, coi=TRUE, plot=TRUE, scale.avg=c(2,8),
                      wavelet.breaks=c(.5,1,2,4) )
## can omit scales of averaging and still visual output.

## some of the things under the hood of ggWavelet which can be called
## by the user.
if (FALSE) {

  ## wavelet transform significance
  nino.wt.power <- abs(nino.wt@transform)^2
  signif = waveSignif(nino.wt, sigtest=0, gws=TRUE, siglvl=.90)  ## note this also returns the gws in the fft.theor

  ## GWS & significance levels (time-averaged):
  gws <- if (signif$gws) signif$fft.theor else rowMeans( nino.wt.power )
  dof = length(nino@POSIXct) - nino.wt@scale   ## the -scale corrects for padding at edges - should this be rolled into waveSignif?
  global.signif = waveSignif(nino.wt, sigtest=1, dof=dof, siglvl=.95)

  nino.wt.2_8 <- scale.avg.wt( nino.wt, 2, 8)
  scaleavg.signif <- waveSignif(nino.wt, sigtest=2, gws=TRUE, siglvl=0.90, dof=c(2,7.9))
  
}
