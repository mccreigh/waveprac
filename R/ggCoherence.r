## name: ggCoherence (method on wavelet coherence object)
## author: james mccreight = mccreigh ^ gmail * com
## info: no warranties or gurantees of any kind.

## General Purpose:

## Return Value:

## Inputs(=defaults):
## :: wc ::
## A wavelet transform object. 

## :: plot=TRUE, coi=TRUE, wavelet.breaks=NULL, font.size=10 ::
## These are all plotting, which only apply if plot=TRUE. This option produces
## graphical output which mimics the IDL output of Torrence & Compo.
## NOTE: This option aims to construct a landscape 8.5x11 inch graphic window
## unless the current device is any of 'pdf','postscript','cairo_pdf','cairo_ps',
## in which case the user is responsible for sizing the output.
## The coi option plots the cone of influence.
## Color levels are specified on the wavelet transform by wavelet.breaks.
## The overall font size in the default plotting is controlled by font.size.

y.fmt <- function(x)  format(x, digits=1, nsmall=2, width=y.width, justify='right') ## y.width lives in the calling level

## ggCoherence method, generic function
if (!isGeneric("ggCoherence")) {  ## creates a generic function
  if (is.function("ggCoherence"))
    fun <- subset else fun <- function(wc, ...) standardGeneric("ggCoherence")
  setGeneric("ggCoherence", fun)
}

setMethod("ggCoherence", c("waveletCoherence"), 
          function( wc, phase.threshold=.5, 
                    plot=TRUE, coi=TRUE, coherence.breaks=c(0.01,.5,.8,.9),
                    phase.subsamp=4,
                    global.coherence=FALSE, global.phase=FALSE,font.size=10, ...) {

  ##############################
  ## y-axis label lenghts calculated up-front
  breaks= rev(a.scale(wc))

  wh.breaks <- which( (breaks %% floor(breaks)) ==0 )
  breaks <- if  (length(wh.breaks)!=0) breaks[wh.breaks] else breaks <- 2^(floor(log(breaks[length(breaks)],base=2)):ceiling(log(breaks[1],base=2)))
  lab.breaks=format(breaks, digits=1, nsmall=2)
  y.width <<- max(nchar(lab.breaks))+1

  coherence.breaks

  ##############################
  ## local coherence
  ntime=length(a.POSIXct(wc))
  nscale=length(a.scale(wc))
  
  coh.frame <-
    data.frame( coherence=as.vector(a.coherence(wc)),
                POSIXct=rep(a.POSIXct(wc),each=nscale),
                period=  rep(a.scale(wc), ntime),
                phase= as.vector(a.phase(wc))
               )

  ## need to NA out approx 3/4 times/columns so the vector field is legible...
  ## also, only show phase where coh>phase.threshold
  ## not sure this is the best approach... a smoothing approach might be better than subsampling... hard to say
  phase.sub <- which( (rep((1:ntime)-1,each=nscale) %% phase.subsamp)==0 & coh.frame$coherence>phase.threshold)
  phase.frame <- coh.frame[phase.sub,]
  
  if (coi) {
    coh.frame$coi <- as.numeric(! (matrix(a.coi(wc), ncol=ntime, nrow=nscale, byrow=TRUE) <
                                    matrix(a.scale(wc), ncol=ntime, nrow=nscale) ) )
  }

  if (a.monteCarlo(wc)!='none') coh.frame$sig <-  as.numeric(as.vector( a.coherence(wc) >= a.qmc.coherence(wc) ))
  
  tile <- if (coi) geom_tile( aes(x=POSIXct, y=period, fill=coherence, alpha=coi))  else
                   geom_tile( aes(x=POSIXct, y=period, fill=coherence) )

  fill.scale <- scale_colour_gradientn(colour = rainbow(7), breaks=coherence.breaks) 
  ##scale_fill_continuous(low="black", high="red", breaks=coherence.breaks) 

  gg.coh.out <- ggplot( coh.frame )

  
  if (plot) TransRevLog2 <<- Trans$new("revlog2", function(x) -(log(x, 2)), function(x) 2^(-x), function(x) bquote(2^.(-x)))  
  if (plot) gg.coh <- ggplot( coh.frame ) +
                         tile +
                         geom_contour( aes(x=POSIXct, y=period, z=coherence), color='white', breaks=c(coherence.breaks) ) +
                         fill.scale +
                         geom_text( data=phase.frame, aes(x=POSIXct,y=period,angle=-1*phase), label='^', size=4, colour='white' ) +
                         scale_y_log2( name='Period ( years )', trans='RevLog2',
                                      labels=lab.breaks, breaks=breaks ) +
                         scale_x_datetime(name='Time',limits=c(min(coh.frame$POSIXct),max(coh.frame$POSIXct))) + 
                         theme_grey(base_size=font.size) +
                         opts(title=paste("Wavelet Coherence and Phase of", a.wt1.name(wc), "&", a.wt2.name(wc),
                                if (a.monteCarlo(wc)!='none')
                                paste('\n(',100*a.confidence(wc),'% confidence from ',a.nsimulations(wc),
                                      ' Monte Carlo sims inside black contour)',sep='') else ''),
                              panel.grid.minor=theme_blank())


  if (plot & coi) gg.coh <- gg.coh + scale_alpha( breaks=c(0,1), labels=c('invalid','valid'), limits=c(-2,1) )

  if (plot & (a.monteCarlo(wc)!='none') ) gg.coh <- gg.coh + geom_contour( aes(x=POSIXct, y=period, z=sig), color='black',
                                                                           breaks=c(.01) )
  
  ##############################
  ## global wavelet spectrum
  global.frame <-
    data.frame( period = a.scale(wc),
               coherence =  a.global_coherence(wc),
               phase= a.global_phase(wc)
               )
  if (a.monteCarlo(wc)!='none') global.frame$confidence <- a.qmc.global_coherence(wc) 
  
  gg.global.out <- ggplot( global.frame )
  gg.global <- gg.global.out
  if (plot & (a.monteCarlo(wc)!='none')) gg.global <- gg.global + geom_line( aes(x=period, y=confidence), linetype=2 )
  if (plot) gg.global <- gg.global +
                           geom_line( aes(x=period, y=coherence) ) +
                           geom_text( aes(x=period,y=coherence,angle=-1*phase), label='^', size=6, colour='red' ) +
                           scale_x_continuous(breaks=breaks, labels=lab.breaks,
                                              name="Period (years)",trans="RevLog2") +
##                           scale_y_continuous(name=paste('Power (',a.data.units(a.input(wc)),'^2)',sep='')) +
                           coord_flip() +
                           theme_grey(base_size=font.size) +
                           opts(title="Global Coherence and Phase")


  ##############################
  ## plot? 
  if (plot) {

    ## open a new window of controlled size for non-file output. this list could need updated.
    if ( !(.Device %in% c('pdf','postscript','cairo_pdf','cairo_ps')) ) dev.new( width=11, height=8.5 )

    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 5)))
    vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)

    ## plot into view ports
    print(gg.coh,   vp=vplayout( 1 , 1:4 ) )
    print(gg.global, vp=vplayout( 1 , 5 ) )
  }

  ##############################
  ## return list of ggplot objects.
  out <- list( wave=gg.coh.out, global=gg.global.out )

  invisible(out)                  
  
}

)
