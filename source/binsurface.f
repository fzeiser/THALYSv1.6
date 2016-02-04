      subroutine binsurface(x1,y1,x2,y2,x3,y3,xl,xu,yl,yu,pi,twopi,
     +                      epsx,epsy,surfloc,surfbin)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 12, 2004
c | Task  : Calculation of the fraction of the bidimensional grid
c |         (xl,yl),(xu,yl),(xu,yu),(xl,yu) covered by the triangle
c |         defined by its summits (x1,y1),(x2,y2) and (x3,y3).
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c epsx           : x-accuracy (used to distinguish 2 points)
c epsy           : y-accuracy (used to distinguish 2 points)
c isurf          : number of points limiting the surface in the bin
c xsurf,ysurf    : arrays containing the coordinates of the points
c                  defining the covered surface in the bin.
c surfloc        : surface to calculate
c
      implicit  none
      real      x1,y1,x2,y2,x3,y3
      real      xl,xu
      real      yl,yu
      real      pi,twopi
      real      epsx,epsy
      dimension xsurf(19),ysurf(19)
      real      xsurf,ysurf
      integer   isurf
      real      surfloc
      real      surfbin
c
c ************************** Local variables ***************************
c
c xr,yr           : help array for storing segment intersection points
c                   with the bin sides
c nr              : number of such intersection points
c ir              : loop counter
c xr1,yr1         : help variable
c xpts,ypts       : help array for storing all intersection points
c ipts            : loop counter
c npts            : total number of intersection points
c xdiff,yddiff    : help variables
c xr2,yr2,xr3,yr3 : help variables
c xg,yg           : center of gravity of the intersection points
c
c
      dimension xr(4),yr(4)
      real      xr,yr
      integer   nr,ir
      real      xr1,yr1
      dimension xpts(12),ypts(12)
      real      xpts,ypts
      integer   ipts,npts
      real      xdiff,ydiff
      real      xr2,yr2,xr3,yr3
      real      xg,yg,weight
      real      norm1,norm2,norm12
      dimension angsurf(19)
      real      angsurf
      real      cos12,ang12,deter12
      real      ang1,ang2
c
c ************************** Logical functions *************************
c
      logical   belongs
      logical   invect
      logical   intri
c
c determine the bin summits included in the considered surface
c
      isurf=0
      if (intri(xl,yl,x1,y1,x2,y2,x3,y3)) then
        isurf=isurf+1
        xsurf(isurf)=xl
        ysurf(isurf)=yl
      endif
      if (intri(xl,yu,x1,y1,x2,y2,x3,y3)) then
        isurf=isurf+1
        xsurf(isurf)=xl
        ysurf(isurf)=yu
      endif
      if (intri(xu,yl,x1,y1,x2,y2,x3,y3)) then
        isurf=isurf+1
        xsurf(isurf)=xu
        ysurf(isurf)=yl
      endif
      if (intri(xu,yu,x1,y1,x2,y2,x3,y3)) then
        isurf=isurf+1
        xsurf(isurf)=xu
        ysurf(isurf)=yu
      endif
      if (isurf.eq.4) then
        surfloc=surfbin
        return
      endif
c
c segment coordinates definition
c
      ipts=0
      call intersection(x1,y1,x2,y2,xl,yl,xu,yu,xr,yr,nr)
      do ir=1,nr
        xr1=xr(ir)
        yr1=yr(ir)
        if (invect(xr1,yr1,x1,y1,x2,y2)) then
          ipts=ipts+1
          xpts(ipts)=xr1
          ypts(ipts)=yr1
        endif
      enddo
      call intersection(x2,y2,x3,y3,xl,yl,xu,yu,xr,yr,nr)
      do ir=1,nr
        xr1=xr(ir)
        yr1=yr(ir)
        if (invect(xr1,yr1,x2,y2,x3,y3)) then
          ipts=ipts+1
          xpts(ipts)=xr1
          ypts(ipts)=yr1
        endif
      enddo
      call intersection(x3,y3,x1,y1,xl,yl,xu,yu,xr,yr,nr)
      do ir=1,nr
        xr1=xr(ir)
        yr1=yr(ir)
        if (invect(xr1,yr1,x3,y3,x1,y1)) then
          ipts=ipts+1
          xpts(ipts)=xr1
          ypts(ipts)=yr1
        endif
      enddo
c
c check if the triangle summits are inside that bin
c
      if (belongs(x1,xu,xl).and.belongs(y1,yu,yl)) then
        ipts=ipts+1
        xpts(ipts)=x1
        ypts(ipts)=y1
      endif
      if (belongs(x2,xu,xl).and.belongs(y2,yu,yl)) then
        ipts=ipts+1
        xpts(ipts)=x2
        ypts(ipts)=y2
      endif
      if (belongs(x3,xu,xl).and.belongs(y3,yu,yl)) then
        ipts=ipts+1
        xpts(ipts)=x3
        ypts(ipts)=y3
      endif
      npts=ipts
c
c intersection analysis : eliminate redundancies
c
c eliminate redundant points
c
      if (isurf.eq.0.and.npts.gt.0) then
        isurf=1
        xsurf(1)=xpts(1)
        ysurf(1)=ypts(1)
      endif
      do 10 ipts=1,npts
        xr1=xpts(ipts)
        yr1=ypts(ipts)
        do ir=1,isurf
          xdiff=abs(xr1-xsurf(ir))
          ydiff=abs(yr1-ysurf(ir))
          if ((xdiff.le.epsx).and.(ydiff.le.epsy)) goto 10
        enddo
        isurf=isurf+1
        xsurf(isurf)=xr1
        ysurf(isurf)=yr1
 10   continue
c
c calculation of the surface defined by the intersection points
c
      if (isurf.le.2) then
        surfloc=0.
        return
      endif
      if (isurf.eq.3) then
        xr1=xsurf(1)
        yr1=ysurf(1)
        xr2=xsurf(2)
        yr2=ysurf(2)
        xr3=xsurf(3)
        yr3=ysurf(3)
        surfloc=abs(0.5*((xr3-xr2)*(yr1-yr2)-(xr1-xr2)*(yr3-yr2)))
        return
      endif
c
c if we have more than 4 points defining the covered surface we have
c to order these points using (xg,yg) as origin and xsurf(1),ysurf(1)
c as (1,0) point
c
      xg=0.
      yg=0.
      do ipts=1,isurf
        xg=xg+xsurf(ipts)
        yg=yg+ysurf(ipts)
      enddo
      weight=1./float(isurf)
      xg=xg*weight
      yg=yg*weight
      xr1=xsurf(1)-xg
      yr1=ysurf(1)-yg
      norm1=xr1*xr1+yr1*yr1
      angsurf(1)=0.
      do ipts=2,isurf
        xr2=xsurf(ipts)-xg
        yr2=ysurf(ipts)-yg
        norm2=xr2*xr2+yr2*yr2
        norm12=norm1*norm2
        if (norm12.eq.0.) then
          cos12=1.
        else
          cos12=(xr1*xr2+yr1*yr2)/(sqrt(norm12))
        endif
        if (cos12.le.-1.0) then
          ang12=pi
          goto 20
        endif
        if (cos12.ge.1.0) then
          ang12=0.
          goto 20
        endif
        ang12=acos(cos12)
 20     deter12=xr1*yr2-xr2*yr1
        if (deter12.lt.0.) ang12=twopi-ang12
        angsurf(ipts)=ang12
      enddo
      do 30 ipts=1,isurf-1
        do 30 ir=ipts+1,isurf
          xr1=xsurf(ipts)
          yr1=ysurf(ipts)
          ang1=angsurf(ipts)
          xr2=xsurf(ir)
          yr2=ysurf(ir)
          ang2=angsurf(ir)
          if (ang2.gt.ang1) goto 30
          xsurf(ipts)=xr2
          ysurf(ipts)=yr2
          angsurf(ipts)=ang2
          xsurf(ir)=xr1
          ysurf(ir)=yr1
          angsurf(ir)=ang1
 30   continue
c
c Determine the covered surface by summing all the sub-triangles
c
      xr1=xsurf(1)
      yr1=ysurf(1)
      xr2=xsurf(isurf)
      yr2=ysurf(isurf)
      surfloc=abs(0.5*((xg-xr2)*(yr1-yr2)-(xr1-xr2)*(yg-yr2)))
      do ipts=1,isurf-1
        xr1=xsurf(ipts)
        yr1=ysurf(ipts)
        xr2=xsurf(ipts+1)
        yr2=ysurf(ipts+1)
        surfloc=surfloc+
     +          abs(0.5*((xg-xr2)*(yr1-yr2)-(xr1-xr2)*(yg-yr2)))
      enddo
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
