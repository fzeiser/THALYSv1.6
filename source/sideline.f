      function sideline(x,y,xs,ys,xe,ye)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 9, 2004
c | Task  : Indicate if (x,y) is on one side (sideline>0.) of the
c |         segment (xe,ye) <--> (xs,ys) or on the other side
c |         (sideline<0.) or on the segment (sideline=0.)
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c x,y  : coordinates of the point to test
c xs,ys: coordinates of the 1st point of the segment
c xe,ye: coordinates of the 2nd point of the segment
c a    : help variable
c b    : help variable
c
      implicit none
      real sideline,x,y,xs,ys,xe,ye,xdlim,ydlim,a,b
c
c ************************** Initialisation ****************************
c
      xdlim=min(abs(xs),abs(xe))/1.e14
      ydlim=min(abs(ys),abs(ye))/1.e14
c
c segment // y-axis2
c
      if (abs(xs-xe).le.xdlim) then
        sideline=x-xs
        if (abs(sideline).le.xdlim) sideline=0.
        return
      endif
c
c segment // x-axis
c
      if (abs(ys-ye).le.ydlim) then
        sideline=y-ys
        if (abs(sideline).le.ydlim) sideline=0.
        return
      endif
c
c normal segment
c
      a=(ys-ye)/(xs-xe)
      b=ys-a*xs
      sideline=y-(a*x+b)
      if (abs(sideline).le.1.e-12) then
        sideline=0.
        return
      endif
      if (sideline.lt.0.) then
        sideline=-1.
        return
      endif
      if (sideline.gt.0.) then
        sideline=1.
        return
      endif
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
