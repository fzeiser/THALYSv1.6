      function vr1(z)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Volume of the projectile-like section.
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function vr1 originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real vr1,z
c
c **********************************************************************
c
      vr1=pi*((2.*r1**3-z**3)/3.+r1**2*z)
      return
      end
