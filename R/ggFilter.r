## name: ggFilter
## purpose: plot original against recon/filtered ts
## plot individual bands against the global wavelet power spectrum.

ggFilter <- function( wt, filter, facet=TRUE, siglvl=.95, use.current.window=FALSE ) {

  ## allow a null filter to be the entire range.
  if (missing(filter)) filter <- waveletFilter( wt, window=a.period(wt)[c(1,length(a.period(wt)))] )
    
  ntime=length(a.POSIXct(a.input(wt)))
  scale.dum <- if (filter$unit=="period") a.period(wt) else a.scale(wt)

  ##############################
  ## original ts against recon from filtered
  orig.recon.frame <- data.frame( POSIXct=a.POSIXct(a.input(wt)),
                                  original=a.data(a.input(wt)),
                                  reconstruct=colSums(filter$filter) )
  orig.recon.frame.melt <- melt(orig.recon.frame, id.var='POSIXct')
  
  gg.orig.recon <- ggplot( orig.recon.frame.melt, aes(x=POSIXct,y=value, linetype=variable, colour=variable) ) +
                     geom_line(size=1) +
                     scale_linetype(name='') + #a.data.name(a.input(wt))) +
                     scale_colour_hue(name='') + #a.data.name(a.input(wt))) +
                     opts(legend.position='left',
                          title=paste(format(100*sd(orig.recon.frame$reconstruct)/sd(orig.recon.frame$orig),
                            nsmall=0,digits=1),'% variance explained.',sep=''))
  
  ##############################
  ## bands plot
  nbands=dim(filter$filter)[1]
  
  filter.frame <- data.frame( t(filter$filter) )
  names(filter.frame)=filter$window.labels.out
  filter.frame$POSIXct=rep(a.POSIXct(a.input(wt)))
  
  filter.frame.melt <- melt(filter.frame, id.vars='POSIXct')

  gg.filter <- ggplot( filter.frame.melt, aes(x=POSIXct, y=value, color=variable)) +
          geom_line(alpha=if (facet) 1 else .5,size=1) +
          scale_color_hue(name=paste(filter$unit,' bands\n(Years)',sep='')) +
            opts(title=paste("Wavelet band filtering of ",a.data.name(a.input(wt)),sep=''),
                 legend.position='left')
  if (facet) gg.filter <- gg.filter + facet_wrap(~variable, ncol=1)

  ##############################
  ## gws plot
  ## y-axis label lenghts calculated up-front
  breaks= rev(a.scale(wt)[which( (0:a.j(wt)) %% (1/a.dj(wt))==0)])
  lab.breaks=format(breaks, digits=1, nsmall=2)
  y.width <<- max(nchar(lab.breaks))+1
  font.size <- 10
  
  dof = length(a.POSIXct(a.input(wt))) - a.scale(wt)   ## the -scale corrects for padding at edges
  global.frame <-
    data.frame( period= a.period(wt),
               power=  rowMeans(abs(a.transform(wt))^2),
               signif= waveSignif(wt, sigtest=1, gws=FALSE, dof=dof, siglvl=siglvl)$signif
               )

  bands=as.vector(t(filter$window.out)) #strsplit(gsub(']','',paste(as.character(filter$cut.out),collapse=' ')),'[:[,(:]')[[1]][2:(2*nbands+1)]
  bands=rep(bands,each=2) ## min min max max
  scales=c(min(global.frame$power), max(global.frame$power)) # max min min max
  scales=c(rev(scales), scales)
  
  band.frame <- data.frame( x=as.numeric(bands),
                            y=scales ,
                            lvl = factor( rep(1:nbands, each=4) ),
                            grp = factor( rep(1:nbands, each=4) ) )

  TransRevLog2 <<- Trans$new("revlog2", function(x) -(log(x, 2)), function(x) 2^(-x), function(x) bquote(2^.(-x)))  

  gg.global <- ggplot( global.frame ) 
  for (bb in 0:(nbands-1)) ## some transform (involving negtives) in the area geom isnt playing well with the revlog2
    gg.global <- gg.global + geom_area( data=band.frame[(bb*4)+(1:4),], aes(x=x,y=y,fill=lvl,group=grp), alpha=.8 )

  gg.global <-  gg.global +
                           geom_line( aes(x=period, y=power) ) +
                           geom_line( aes(x=period, y=signif), linetype=2 ) +
                           scale_x_continuous(breaks=breaks, labels=lab.breaks,
                                              name="Period (years)",trans="RevLog2") +
                           scale_y_continuous(name=paste('Power (',a.data.units(a.input(wt)),'^2)',sep='')) +
                           coord_flip() +
                           theme_grey(base_size=font.size) +
                           opts(title=paste("Global Wavelet\nSpectrum,\n",
                                             siglvl,' significance.',sep=''), legend.position='none')   

  if (!use.current.window) check.new.dev()  ## open a new device?

  ## set up view ports
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, 4)))
  vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
  
  ## plot into view ports
  print(gg.orig.recon, vp=vplayout( 1, 1:3) )
  print(gg.filter,  vp=vplayout( 2:5 , 1:3   ) )   
  print(gg.global, vp=vplayout( 1:5, 4))    

  invisible( list( orig.recon.ts=gg.orig.recon, filter.ts=gg.filter, global=gg.global) )
}  
