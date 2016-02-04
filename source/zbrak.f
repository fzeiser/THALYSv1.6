      subroutine zbrak(func,x1,x2,n,xb1,xb2,nb)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning (adapted from Numerical Recipes)
c | Date  : June 26, 2007
c | Task  : Bracket the function
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer  n,nb,i,nbb
      real     func,x1,x2,xb1(nb),xb2(nb),x,dx,fp,fc
c
c ************************** Bracket function **************************
c
c The function has a zero between xb1 and xb2.
c
c func: function
c
      nbb=0
      x=x1
      dx=(x2-x1)/n
      fc=0.
      fp=func(x)
      do 10 i=1,n
        x=x+dx
        fc=func(x)
        if (fc*fp.lt.0.) then
          nbb=nbb+1
          xb1(nbb)=x-dx
          xb2(nbb)=x
          if (nbb.eq.nb) goto 20
        endif
        fp=fc
   10 continue
   20 nb=nbb
      return
      end
