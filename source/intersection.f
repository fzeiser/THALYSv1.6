      subroutine intersection(xs,ys,xe,ye,xl,yl,xu,yu,xr,yr,ir)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 9, 2004
c | Task  : Calculation of the intersection points defined by a segment
c |         (xs,ys) <--> (xe,ye) crossing a bin defined by its summits
c |         (xl,yl),(xl,yu),(xu,yu),(xu,yl).
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c xs,ys       : coordinates of the segment first point
c xe,ye       : coordinates of the segment second point
c xl,yl,xu,yu : coordinates defining a bin
c xr,yr       : array storing the correct coordinates
c x1,y1,x2,y2 : coordinates of intersection points inside the bin
c ir          : counter for xr,yr
c
c a,b     : coefficient of the segment line if segment // neither
c           to the x-axis nor to the y-axis
c xa,ya   : intersection of the segment with the vertical axis linking
c           (xl,yl) with (xl,yu)
c xb,yb   : intersection of the segment with the vertical axis linking
c           (xu,yl) with (xu,yu)
c xc,yc   : intersection of the segment with the horizontal axis linking
c           (xl,yl) with (xu,yl)
c xd,yd   : intersection of the segment with the horizontal axis linking
c           (xl,yu) with (xu,yu)
c belongs : logical function to test if one value is between two others
c
      implicit none
      logical belongs
      integer ir
      real    xs,ys,xe,ye,xl,yl,xu,yu,xr(4),yr(4),xdlim,ydlim,a,b,
     +        xa,xb,xc,xd,ya,yb,yc,yd
c
c **************************** Initialisation **************************
c
      ir=0
      xdlim=min(abs(xs),abs(xe))/1.e14
      ydlim=min(abs(ys),abs(ye))/1.e14
c
c segment // y-axis ==> 2 cases
c
      if (abs(xs-xe).le.xdlim) then
        if (belongs(xs,xu,xl)) then
          ir=2
          xr(1)=xs
          xr(2)=xs
          yr(1)=yu
          yr(2)=yl
        endif
        return
      endif
c
c segment // x-axis
c
      if (abs(ys-ye).le.ydlim) then
        if (belongs(ys,yu,yl)) then
          ir=2
          xr(1)=xu
          xr(2)=xl
          yr(1)=ys
          yr(2)=ys
        endif
        return
      endif
c
c normal segment crossing the bin
c
      a=(ys-ye)/(xs-xe)
      b=ys-a*xs
c intersections with vertical axis
      xa=xl
      ya=a*xl+b
      xb=xu
      yb=a*xu+b
c intersections with horizontal axis
      yc=yl
      xc=(yl-b)/a
      yd=yu
      xd=(yu-b)/a
c only keep intersections inside the bin limits
      ir=0
      if (belongs(ya,yu,yl)) then
        ir=ir+1
        xr(ir)=xa
        yr(ir)=ya
      endif
      if (belongs(yb,yu,yl)) then
        ir=ir+1
        xr(ir)=xb
        yr(ir)=yb
      endif
      if (belongs(xc,xu,xl)) then
        ir=ir+1
        xr(ir)=xc
        yr(ir)=yc
      endif
      if (belongs(xd,xu,xl)) then
        ir=ir+1
        xr(ir)=xd
        yr(ir)=yd
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
