      function rhodi(z)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Returns the shape of the dinuclear system.
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function rho originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real rhodi,z,d
c
c **********************************************************************
c
      d=totl-r1-r3
      if (z.lt.-r1.or.z.gt.d+r3) go to 3
      if (z.gt.z1) go to 1
      rhodi=sqrt(r1**2-z**2)
      return
 1    if (z.gt.z3) go to 2
      rhodi=r2+aaa*aaa*cur*(cosh((z-z2)/aaa)-1.)
      return
 2    rhodi=sqrt(r3**2-(z-d)**2)
      return
 3    rhodi=0.
      return
      end
