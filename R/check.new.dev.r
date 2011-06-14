check.new.dev <- function( height=8.5, width=11) {
  if ( !(.Device %in% c('pdf','postscript','cairo_pdf','cairo_ps')) ) dev.new( width=width, height=height )
}
