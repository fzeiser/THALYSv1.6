      subroutine pol1(x1,x2,y1,y2,x,y)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 5, 2004
c | Task  : Interpolation of first order
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      real x1,x2,y1,y2,fac,x,y
c
c ***************************** Interpolation **************************
c
      fac=(x-x1)/(x2-x1)
      y=y1+fac*(y2-y1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
