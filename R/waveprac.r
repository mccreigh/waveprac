## The primary internal setup for waveprac

## imports - should I be more specific??
require('stats')
require('ggplot2')
require('RColorBrewer')
require('plyr')
require('multicore')
require('grDevices')

## EOF depends on waveprac (any routines using both should then be sourced in this file.)

## variable definition since lubridate's ruins normal POSIXct functionality
origin=as.POSIXct( 0, origin='1970-01-01', tz='UTC' )

## some more widely used setClassUnions
setClassUnion("numericNULL", c("numeric", "NULL"))
setClassUnion("vectorNULL", c("vector", "NULL"))
setClassUnion("arrayNULL", c("array", "NULL"))

## ####################################################################
## class: timeSeries
setClass("timeSeries",
         representation(data = "vector",
                        data.name = "character",
                        data.units = "character",
                        POSIXct = "POSIXct")
         )

## spaceTime accessor functions
set.accessors( 'timeSeries' )

## ####################################################################
## class: spaceTime
setClass("spaceTime",
         representation(data = "array",
                        data.name = "character",
                        data.units = "character",
                        lon = "vector", 
                        lat = "vector", 
                        POSIXct = "POSIXct")
         )

## spaceTime accessor functions
set.accessors( 'spaceTime' )

## ####################################################################
## class: spaceTimeBands
## note that his S4 class contains a list so the
## generic accessors are mixed when accessing the list
## @ for S4 objects and $ for the list.
setClass("spaceTimeBands",
         representation(data = "list",
                        data.name = "character",
                        data.units = "character",                        
                        lon = "vector", 
                        lat = "vector", 
                        POSIXct = "POSIXct",
                        bands="matrix",
                        s0="numeric",
                        j="numeric",
                        dj="numeric"
                        )
         )
## spaceTimeBands accessor functions
set.accessors( 'spaceTimeBands' )


## these were sourced but are in EOF... 
## gg.ccf.r
##"spectralCEOF.r"
##"ggSpectralCEOF.r"
