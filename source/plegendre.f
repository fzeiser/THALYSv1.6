      function plegendre(l,x)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : July 2, 2004
c | Task  : Calculation of Legendre polynomial
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      implicit none
      integer  l,i
      real     plegendre,x,pl(0:200),pl0,pl1
c
c ************************ Legendre polynomial *************************
c
c x           : cosine of the angle
c plegendre,pl: legendre polynomial
c l           : l-value
c pl0, pl1    : help variables
c
      pl0=1.
      pl1=x
      pl(0)=pl0
      pl(1)=pl1
      do 10 i=2,l
         pl(i)=(x*(2*i-1)*pl1-(i-1)*pl0)/i
         pl0=pl1
         pl1=pl(i)
 10   continue
      plegendre=pl(l)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
