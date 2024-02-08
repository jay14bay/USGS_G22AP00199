# USGS_G22AP00199
Jeff Bayless
jeff.bayless@aecom.com

MATLAB implementation of Bea24 rupture directivity model:
 
Bea24_Example.m: a script which calculates the directivity effect (median and standard deviation) for an example scenario and produces maps of the model components and the median adjustment, making use of the functions below.
GC2.m: Calculates the GC2 coordinates. This function is a conversion to MATLAB of Brian Chiou's R functions (pers. comm.)
Bea24.m: a function which implements the Bea24 directivity model 

HAZ45.2 implementation of Bea24:

Directivity_bea24.f: a fortran subroutine which implements the directivity model.
Modified versions of the following HAZ45.2 programs and subroutines: Directivity.f, Haz_main2.f, gc2.f, cldist.f, and declare1.h
ChangeLog-7Feb2024.txt: a list of the required changes from HAZ45.2 required to implement the model. 