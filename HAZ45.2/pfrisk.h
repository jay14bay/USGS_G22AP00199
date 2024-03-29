c     Set array dimensions

      integer MAX_FLT, MAX_SEG, MAX_INTEN, MAX_PROB, 
     1        MAX_EPS, MAX_Xcost, MAXPARAM, MAX_MAG, 
     2        MAX_DIST, MAX_N1, MAX_WIDTH, MAX_DIST1, 
     3        MAX_GRID, MAX_SYN, MAX_AMPMAG, MAX_AMPPER,
     4        MAX_AMPGM, MAX_PER, MAXDETM_DIST, MAX_DD
      integer MAXFLT_DD, MAXFLT_AS, MAX_NODE,
     1        MAX_ATTEN, MAX_FTYPE, MAX_S7
      integer MAX_NVERT

      parameter (MAX_FLT=300, MAX_SEG=420, MAX_INTEN=21, MAX_PROB=21,
     1           MAX_EPS=9, MAX_Xcost=7, MAXPARAM=80, MAX_MAG=10,
     2           MAX_DIST=32, MAX_N1=120, MAX_WIDTH=9, MAX_DIST1=10000,
     3           MAX_GRID=170000, MAX_SYN=10, MAX_AMPMAG=25, MAX_AMPPER=15, 
     4           MAX_AMPGM=15, MAX_PER=101, MAXDETM_DIST=2000, MAX_DD=5) 
      parameter (MAXFLT_DD=1000, MAXFLT_AS=2500, MAX_NODE=200,
     1           MAX_ATTEN=143, MAX_FTYPE=5, MAX_S7=70000)
      parameter (MAX_NVERT=50)
