This is an R port and extension of

"A Practical Guide to Wavelet Analysis"  
by C. Torrence and and G. P. Compo,  
1998: Bull. Amer. Meteor. Soc., 79, 61-78.  
http://paos.colorado.edu/research/wavelets/  

This package is not ready to build. I'm using devtools to develop it, see: https://github.com/hadley/devtools . Once you download and unpack it, command line: `$ R CMD INSTALL hadley-devtools-1234567`

In R:  
`require(devtools)`  
then see examples.r for loading the development package and data and running some things. 

I'm not sure if DESCRIPTION will force its dependencies. The current package depenedencies are ggplot2, plyr, multicore (if you want to use more than 1 core on a multicore machine).

See DESCRIPTION for logistics and credits.   
See TODO for things that need to be done.

The current version of waveprac is early beta. Much needs to be done.   
I have checked that I reproduce the results of the IDL code of Torrence and Compo for the basic wavelet transforms. However this should be formalized as a test. 

Please feel free to let me know about issues, suggestions, new code, etc... this is github.



