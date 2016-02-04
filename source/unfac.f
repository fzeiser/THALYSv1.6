      subroutine unfac(l,rho,rhoc,amun,vl,ps)
c
c +---------------------------------------------------------------------
c | Author: Gilles Noguere
c | Date  : September 14, 2011
c | Task  : Penetrability factor (vl) and phase shift (ps) from NJOY
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer l
      real    rho,rhoc,amun,vl,ps,r2,r4
c
c **********************************************************************
c
      vl=1.
      ps=1.
      r2=rho*rho
      if (l.eq.0) then
        vl=amun
        ps=rhoc
      else if (l.eq.1) then
        vl=amun*r2/(1.+r2)
        ps=rhoc-atan(rhoc)
      else if (l.eq.2) then
        r4=r2*r2
        vl=amun*r4/(9.+3.*r2+r4)
        ps=rhoc-atan(3.*rhoc/(3.-rhoc*rhoc))
      endif
      return
      end
