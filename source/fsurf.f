      function fsurf(epscloc)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Form factor for the surface energy
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function fsurf originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      real epscloc,dum,fsurf
c
c **********************************************************************
c
      if (epscloc) 10,20,30
 10   dum=sqrt(1.+epscloc**2)
      fsurf=(alog(-epscloc+dum)-epscloc*dum)/(-2.*epscloc*dum**(1./3.))
      return
 20   fsurf=1.
      return
 30   dum=sqrt(1.-epscloc**2)
      fsurf=(asin(epscloc)+epscloc*dum)/(2.*epscloc*dum**(1./3.))
      return
      end
