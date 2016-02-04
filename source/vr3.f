      function vr3(z)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Volume of the target-like section.
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function vr3 originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real vr3,z,d
c
c **********************************************************************
c
      d=totl-r1-r3
      vr3=pi*((2.*r3**3-(d-z)**3)/3.+r3**2*(d-z))
      return
      end
