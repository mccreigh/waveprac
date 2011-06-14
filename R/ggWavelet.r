#' name: ggWavelet (method on wavelet transform object)
#' author: james mccreight = mccreigh ^ gmail * com
#' info: no warranties or gurantees of any kind.

#' General Purpose:
#' This function takes a wavelet transform object (wt) and prepares it for
#' display using ggplot2, returning a list og ggplot 2 objects which can
#' be subjected to.
#' It gives optional, default graphical output mimicing that of Torrence
#' and Compo's IDL routines.

#' Return Value:
#' A list with fundamental (atmoic) ggplot plot structures, which can
#' then be customized by the user with ggplot syntax. The list contains
#' "input" , "wave" and "global" entries, corresponding to the dataset
#' subject to wavelet transform, its wavelet transform, and the global
#' wavelet spectra. Examples of how these can be use are provided in
#' this code by the default plot options.

#' Inputs(=defaults):
#' :: wt ::
#' A wavelet transform object. 

#' :: scale.avg=FALSE ::
#' If scale.avg is supplied as a length 2 numeric vector, a scale entry
#' is also returned in the output list which contains the information for
#' plotting the scale average wavelet power specturm.

#' :: siglvl=.95, siglvl.global=siglvl, siglvl.scale=siglvl ::
#' Significance levels are calculated for each list entry using a siglvl
#' which applies to the wavelet transform. If siglvl.global and siglvl.scale
#' are not specified individually, siglvl is used in their absence.

#' :: gws==FALSE ::
#' When computing the significance levels for the wavelet transform, use the
#' global wavelet spectrum as a theoretical background spectrum. If not true,
#' then a spectrum must/can be supplied in its place.

#' :: plot=TRUE, coi=TRUE, wavelet.breaks=NULL, font.size=10 ::
#' These are all plotting, which only apply if plot=TRUE. This option produces
#' graphical output which mimics the IDL output of Torrence & Compo.
#' NOTE: This option aims to construct a landscape 8.5x11 inch graphic window
#' unless the current device is any of 'pdf','postscript','cairo_pdf','cairo_ps',
#' in which case the user is responsible for sizing the output.
#' The coi option plots the cone of influence.
#' Color levels are specified on the wavelet transform by wavelet.breaks.
#' The overall font size in the default plotting is controlled by font.size.

y.fmt <- function(x)  format(x, digits=1, nsmall=2, width=y.width, justify='right') ## y.width lives in the calling level


ggWavelet <- function( wt, siglvl=.95, siglvl.global=siglvl, siglvl.scale=siglvl, scale.avg=FALSE, gws=FALSE,
                      plot=TRUE, coi=TRUE, wavelet.breaks=NULL, font.size=10, use.current.window=FALSE, ...) {

  ##############################
  ## y-axis label lenghts calculated up-front
  breaks= rev(a.scale(wt)[which( (0:a.j(wt)) %% (1/a.dj(wt))==0)])
  lab.breaks=format(breaks, digits=1, nsmall=2)
  y.width <<- max(nchar(lab.breaks))+1
  
  ##############################
  ## the input timeseries
  input.frame <- data.frame( time=a.POSIXct(a.input(wt)),
                             data=a.data(a.input(wt)) )

  gg.input.out <- ggplot( input.frame )
  if (plot) gg.input <- gg.input.out + 
                          geom_point( aes(x=time, y=data), color='red') +
                          geom_line( aes(x=time, y=data) ) +
                          scale_y_continuous(name=paste(a.data.name(a.input(wt)),
                                                        ' (', a.data.units(a.input(wt)),')',sep=''),
                                             formatter=y.fmt) +
                          theme_grey(base_size=font.size) +
                          opts(title=paste(a.data.name(a.input(wt))), legend.position='none')

  ##############################
  ## local wavelet spectrum
  if (!exists('waveSignif')) source("waveSignif.r")
  ntime=length(a.POSIXct(a.input(wt)))
  nscale=length(a.scale(wt))
  
  wave.frame <-
    data.frame( POSIXct=rep(a.POSIXct(a.input(wt)),each=nscale),
                period=  rep(a.period(wt), ntime),
                power=   as.vector(abs(a.transform(wt))^2),
                signif=  as.vector(waveSignif(wt, sigtest=0, gws=gws, siglvl=siglvl, ...)$ratio)
               )
  wave.frame$time <- wave.frame$POSIXct

  if (coi) {
    wave.frame$coi <- as.numeric(! (matrix(a.coi(wt), ncol=ntime, nrow=nscale, byrow=TRUE) <
                                    matrix(a.period(wt), ncol=ntime, nrow=nscale) ) )
  }

  tile <- if (coi) geom_tile( aes(x=time, y=period, fill=power, alpha=coi)) else
                   geom_tile( aes(x=time, y=period, fill=power) )
  
  fill.scale <- if (is.null(wavelet.breaks)) {
    scale_fill_continuous(low="black", high="red")
  } else {
    scale_fill_continuous(low="black", high="red", breaks=wavelet.breaks) 
  }

  gg.wave.out <- ggplot( wave.frame )

  if (plot) TransRevLog2 <<- Trans$new("revlog2", function(x) -(log(x, 2)), function(x) 2^(-x), function(x) bquote(2^.(-x)))  
  if (plot) gg.wave <- ggplot( wave.frame ) +
                         tile +
                         geom_contour( aes(x=time, y=period, z=signif), color='white', breaks=c(siglvl) ) +
                         fill.scale +
                         scale_y_log2( name='Period ( years )', trans='RevLog2',
                                      labels=lab.breaks, breaks=breaks ) +
                         scale_x_datetime(limits=c(min(wave.frame$time),max(wave.frame$time))) + 
                         theme_grey(base_size=font.size) +
                         opts(title=paste("Wavelet Power Spectrum (", a.data.units(a.input(wt)),'^2),',
                                          ' contours at ', siglvl,' significance.',sep=''),
                              legend.position = "none", panel.grid.minor=theme_blank())

  if (plot & coi) gg.wave <- gg.wave+ scale_alpha( breaks=c(0,1), labels=c('invalid','valid'), limits=c(-2,1) )   

  ##############################
  ## global wavelet spectrum
  dof = length(a.POSIXct(a.input(wt))) - a.scale(wt)   ## the -scale corrects for padding at edges
  global.frame <-
    data.frame( period= a.period(wt),
                power=  rowMeans(abs(a.transform(wt))^2),
                signif= waveSignif(wt, sigtest=1, gws=FALSE, dof=dof, siglvl=siglvl.global, ...)$signif
               )

  gg.global.out <- ggplot( global.frame )
  if (plot) gg.global <- gg.global.out +
                           geom_line( aes(x=period, y=power) ) +
                           geom_line( aes(x=period, y=signif), linetype=2 ) +
                           scale_x_continuous(breaks=breaks, labels=lab.breaks,
                                              name="Period (years)",trans="RevLog2") +
                           scale_y_continuous(name=paste('Power (',a.data.units(a.input(wt)),'^2)',sep='')) +
                           coord_flip() +
                           theme_grey(base_size=font.size) +
                           opts(title=paste("Global, ", siglvl.global,' significance.',sep=''))   

  ##############################
  ## scale-average power time series : if scale.avg is NOT logical (FALSE by default)
  do.scale <- !(typeof(scale.avg)=='logical')  ## this should probably be a length and type=numeric check

  if (do.scale) {
    if (!exists('scale.avg.wt')) source("scale.avg.wt.r")  
    scale.frame <-
      data.frame( time=a.POSIXct(a.input(wt)),
                 variance=scale.avg.wt( wt, scale.avg[1], scale.avg[2]),
                 signif=waveSignif(wt, sigtest=2, gws=gws, siglvl=siglvl.scale, dof=scale.avg, ...)$signif
                 )
    gg.scale.out <- ggplot( scale.frame )
    if (plot) gg.scale <- gg.scale.out +
                            geom_line( aes(x=time, y=variance) ) +
                            geom_line( aes(x=time, y=signif), linetype=2) +
                            scale_y_continuous(name=paste("Average variance (",a.data.units(a.input(wt)),'^2)',sep=''),
                                               formatter=y.fmt) +
                            theme_grey(base_size=font.size) +
                            opts(title=paste('Scale average (',scale.avg[1],'-',scale.avg[2],' years) power, ',
                                              siglvl.scale,' significance.',sep=''), legend.position='none')    
  }

  ##############################
  ## plot? 
  if (plot) {
    ## credits viewport, yay chris and gil!
    ct.gc <- ggplot() + annotate("text", size=font.size*.37, x = 0, y = 0,
                                 label = paste("A Practical Guide to\nWavelet Analysis",
                                               " \n\nby\nC. Torrence & G.P. Compo",
                                               " \n\nhttp://paos.colorado.edu/\nresearch/wavelets",sep='')) +
               theme_bw() +
               opts(axis.ticks.length=unit(0,'lines'), axis.text.x=theme_blank(), axis.text.y=theme_blank(),
                    axis.title.x=theme_blank(), axis.title.y=theme_blank(),
                    panel.grid.major=theme_blank(), panel.grid.minor=theme_blank()
                    )


    ## open a new window of controlled size for non-file output. this list could need updated.
    if (!use.current.window) check.new.dev()  ## open a new device?

    ## set up view ports
    ncol=4 
    nrow <- if (do.scale) 4 else 3
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow, ncol)))
    vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)

    ## plot into view ports
    print(gg.input,  vp=vplayout( 1 , 1:3   ) )   
    print(ct.gc, vp=vplayout( 1, 4))    
    print(gg.wave,   vp=vplayout( 2:3 , 1:3 ) )
    print(gg.global, vp=vplayout( 2:3 , 4 ) )
    if (do.scale) print(gg.scale,  vp=vplayout( 4 , 1:3 ) )
    
  }

  ##############################
  ## return list of ggplot objects.
  out <- list( input=gg.input.out, wave=gg.wave.out, global=gg.global.out )
  if (do.scale) out$scale <- gg.scale.out

  invisible(out)                  
  
}

