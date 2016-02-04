      function radius(a)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning (adapted from R.K. Tripathi)
c | Date  : July 7, 2004
c | Task  : Radius function for Tripathi formula
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer na(23),ia,i
      real    radius,a,rms(23),fact
c
c *************************** Radius function **************************
c
      data na/1,2,3,4,6,7,9,10,11,12,13,14,15,16,
     +  17,18,19,20,22,23,24,25,26/
      data rms/0.85,2.095,1.976,1.671,2.57,2.41,2.519,2.45,2.42,
     +  2.471,2.440,2.58,2.611,2.730,2.662,2.727,2.900,3.040,2.969,2.94,
     +  3.075,3.11,3.06/
      fact=sqrt(5./3.)
      ia=int(a+0.4)
      radius=fact*(0.84*a**(1./3.)+0.55)
      do 10 i=1,23
        if (ia.eq.na(i)) radius=fact*rms(i)
  10  continue
      return
      end
