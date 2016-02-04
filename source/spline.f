      subroutine spline(x,y,nk,yp1,ypn,y2)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Spline fit
c +---------------------------------------------------------------------
c
c *************************** Comments ********************************
c
c  (C) Copr. 1986-92 Numerical Recipes Software %#..5#2P.
c
c ****************** Declarations and common blocks ********************
c
      integer nk
      real yp1,ypn,x(nk),y(nk),y2(nk)
      integer i,k
      real psp,qn,sig,un,u(500)
c
c **********************************************************************
c
      if (yp1.gt.0.99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 10 i=2,nk-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        psp=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/psp
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     +    (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/psp
   10 continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(nk)-x(nk-1)))*(ypn-(y(nk)-y(nk-1))/(x(nk)-x(nk-1)))
      endif
      y2(nk)=(un-qn*u(nk-1))/(qn*y2(nk-1)+1.)
      do 20 k=nk-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
20    continue
      return
      end
