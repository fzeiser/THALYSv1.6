      subroutine molprepare(tjl,numtjl,st,npmold,xmo,wmo,tav,vnu,
     +  product,numtr,numinc)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : March 2, 2010
c | Task  : Preparation of Moldauer width fluctuation correction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numtjl,npmold,numtr,numinc,i,im
      real             xmo(npmold),wmo(npmold),vnu(numtr),
     +                 product(npmold),x,expo,capt,factor,eps,yy,x1,
     +                 prod,fxmsqrt,fxmold
      double precision tjl(0:5,numtr),st,tav(numtr),tgamma
c
c **************** Preparation of Moldauer integral ********************
c
c tgamma: gamma transmission coefficient
c numtjl: number of transmission coefficients
c tav   : average transmission coefficients
c vnu   : number of degrees of freedom
c tjl   : transmission coefficients
c
c Initialisation and average transmission coefficients
c
      tgamma=0.
      do 10 i=1,numtr
        tav(i)=0.
        vnu(i)=1.
   10 continue
      do 20 i=1,numtjl
        if (tjl(0,i).gt.0.) tav(i)=tjl(1,i)/tjl(0,i)
   20 continue
      tgamma=tjl(1,numtjl+1)
c
c Calculation of number of degrees of freedom
c
c st : denominator of compound nucleus formula
c
      do 30 i=1,numtjl
        vnu(i)=1.78+real(((tav(i)**1.212)-0.78)*exp(-0.228*st))
   30 continue
      vnu(numtjl+1)=1.
c
c Loop over integration points
c
c x,expo,capt: help variables
c
      do 40 im=1,npmold
        x=xmo(im)
c
c Special case for capture
c
        expo=real(tgamma*x/st)
        if (expo.gt.80.) then
          capt=0.
        else
          capt=exp(-expo)
        endif
c
c Calculation of product over tjl
c
c factor,eps,yy,x1   : help variables
c prod,fxmsqrt,fxmold: help variables
c product            : product used in final Moldauer calculation
c
        factor=0.
        do 50 i=numinc+1,numtjl
          eps=real(2.*tav(i)*x/(st*vnu(i)))
          yy=real(-vnu(i)*0.5*tjl(0,i))
          if (eps.le.1.e-30) goto 50
          if (eps.lt.1.e-5) then
            x1=eps-0.5*eps*eps+(eps**3)/3.-0.25*(eps**4)+0.2*(eps**5)
          else
            x1=log(1.+eps)
          endif
          factor=factor+yy*x1
   50   continue
        prod=exp(factor)
        fxmsqrt=wmo(im)*exp(0.5*x)
        fxmold=fxmsqrt*prod*fxmsqrt
        product(im)=fxmold*capt
   40 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
