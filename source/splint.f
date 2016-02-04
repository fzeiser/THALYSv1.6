      subroutine splint(xa,ya,y2a,i,x,y)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : May 12, 2008
c | Task  : Spline fit
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c  (C) Copr. 1986-92 Numerical Recipes Software %#..5#2P.
c
c ****************** Declarations and common blocks ********************
c
      integer i,k,khi,klo
      real    a,b,hsp,x,y,xa(i),y2a(i),ya(i)
c
c **********************************************************************
c
      klo=1
      khi=i
   10 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
          khi=k
        else
          klo=k
        endif
        goto 10
      endif
      hsp=xa(khi)-xa(klo)
      if (hsp.eq.0.) then
        write(*,*) 'bad xa input in splint'
        stop
      endif
      a=(xa(khi)-x)/hsp
      b=(x-xa(klo))/hsp
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*
     +  (hsp**2)/6.
      return
      end
