      function vr2(zz1,zz2,zz3)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Volume of the neck
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function vr2 originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real vr2,zz1,zz2,zz3,z21,z32,b
c
c **********************************************************************
c
      z21=(zz2-zz1)/aaa
      z32=(zz3-zz2)/aaa
      if (z21.gt.30.) z21=30.
      if (z32.gt.30.) z32=30.
      b=aaa*aaa*cur
      vr2=pi*(((r2-b)**2+.5*b**2)*(zz3-zz1)+
     +  aaa*b*(2.*(r2-b)*(sinh(z32)+sinh(z21))+
     +  0.25*b*(sinh(2.*z32)+sinh(2.*z21))))
      return
      end
