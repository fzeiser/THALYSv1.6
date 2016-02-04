      subroutine pol2(x1,x2,x3,y1,y2,y3,x,y)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : Polynomial interpolation of second order
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      real    x1,x2,x3,y1,y2,y3,x,y,yy1,yy2,yy3
c
c ******************** Lagrange formula for order 2 ********************
c
c         (x-x2)(x-x3)        (x-x1)(x-x3)        (x-x1)(x-x2)
c    y = --------------.y1 + --------------.y2 + --------------.y3
c        (x1-x2)(x1-x3)      (x2-x1)(x2-x3)      (x3-x1)(x3-x2)
c
      yy1=(x-x2)*(x-x3)/((x1-x2)*(x1-x3))*y1
      yy2=(x-x1)*(x-x3)/((x2-x1)*(x2-x3))*y2
      yy3=(x-x1)*(x-x2)/((x3-x1)*(x3-x2))*y3
      y=yy1+yy2+yy3
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
