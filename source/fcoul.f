      function fcoul(epscloc)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Form factor for the coulomb self energy
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function fcoul originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      real epscloc,fcoul
c
c **********************************************************************
c
      if (epscloc) 10,20,30
 10   fcoul=(1.+epscloc**2)**(1./3.)/epscloc*atan(epscloc)
      return
 20   fcoul=1.
      return
 30   fcoul=(1.-epscloc**2)**(1./3.)/(2.*epscloc)*
     +  alog((1.+epscloc)/(1.-epscloc))
      return
      end
