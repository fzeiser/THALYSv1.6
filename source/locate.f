      subroutine locate(xx,ib,ie,x,j)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning (adapted from Numerical Recipes)
c | Date  : August 24, 2004
c | Task  : Find value in ordered table
c +---------------------------------------------------------------------
c
c **************************** Declarations ****************************
c
      implicit none
      logical  ascend
      integer  ib,ie,j,jl,jm,ju
      real     x,xx(0:ie)
c
c ******************************* Search *******************************
c
c Find j such that xx(j) <= x < xx(j+1) or xx(j) > x >= xx(j+1)
c
      j=0
      if (ib.gt.ie) return
      jl=ib-1
      ju=ie+1
      ascend=xx(ie).ge.xx(ib)
  10  if (ju-jl.gt.1) then
        jm=(ju+jl)/2
        if (ascend.eqv.(x.ge.xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
        goto 10
      endif
      if (x.eq.xx(ib)) then
        j=ib
      else if (x.eq.xx(ie)) then
        j=ie-1
      else
        j=jl
      endif
      return
      end
