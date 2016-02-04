      function prodp(sx,tjl,numtjl,numtr,numinc)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 10, 2004
c | Task  : Calculation of the product of (1+t*x)**(1/20)
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numtjl,numtr,numinc,i
      real             sx
      double precision tjl(0:5,numtr),tav,prodp
c
c **************** Calculation of the product of 1+t*x *****************
c
c prodp : product function for GOE
c numtjl: number of transmission coefficients
c tav   : average transmission coefficients
c
      prodp=1.
      do 10 i=numinc+1,numtjl
        if (tjl(0,i).eq.0.) then
          tav=0.
        else
          tav=tjl(1,i)/tjl(0,i)
        endif
        prodp=prodp*real((1.+tav*sx)**(tjl(0,i)*0.05))
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
