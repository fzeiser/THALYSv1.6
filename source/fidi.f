      function fidi(x,i)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : May 8, 2007
c | Task  : fidi=0 determines the exact shape of the dinuclear complex.
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function fidi originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer i
      real    fidi,x(i),d,wu1,wu3,sq1,sq3,z21,z32,vr3,vr2,vr1,al,be,ga,
     +        de,eq,epsc,ds
c
c **********************************************************************
c
      d=totl-r1-r3
      aaa=.8+x(1)*2.
      z2=d*(.5-x(2))
      z1=r1*x(3)
      z3=d-r3*x(4)
      r1=rp*x(5)
      r3=rt*x(6)
      cur=c0*abs(x(7))
      wu1=r1**2-z1**2
      if (wu1.lt.1.e-8) wu1=1.e-8
      wu3=r3**2-(d-z3)**2
      if (wu3.lt.1.e-8) wu3=1.e-8
      sq1=sqrt(wu1)
      sq3=sqrt(wu3)
      z21=(z2-z1)/aaa
      z32=(z3-z2)/aaa
      if (z21.gt.30.) z21=30.
      if (z32.gt.30.) z32=30.
      al=(vtot-vr1(z1)-vr2(z1,z2,z3)-vr3(z3))**2
      be=((rp**3-r1**3)-rpt*(rt**3-r3**3))**2*1.e1
      ga=(sq1-r2-aaa*aaa*cur*(cosh(z21)-1.))**2*1.e5
      de=(sq3-r2-aaa*aaa*cur*(cosh(z32)-1.))**2*1.e5
      eq=(z1/sq1-aaa*cur*sinh(z21))**2*1.e5
      epsc=((d-z3)/sq3-aaa*cur*sinh(z32))**2*1.e5
      ds=x(7)**2*1.e-4
      fidi=al+be+ga+de+eq+epsc+ds
      return
      end
