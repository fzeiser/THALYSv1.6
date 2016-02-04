      subroutine moldauer(numtjl,na,nb,st,npmold,xmo,tav,vnu,product,
     +  res,numtr,ielas)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : August 22, 2004
c | Task  : Moldauer width fluctuation correction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numtjl,na,nb,npmold,numtr,ielas,im
      real             xmo(npmold),vnu(numtr),product(npmold),res,
     +                 dab,factor,factora,factorb,x
      double precision st,tav(numtr)
c
c **************** Calculation of Moldauer integral ********************
c
c numtjl           : number of transmission coefficients
c na,nb            : counters for width fluctuation calculation
c st               : denominator of compound nucleus formula
c npmold,xmo       : variables for Gauss-Legendre integration
c tav              : average transmission coefficients
c vnu              : number of degrees of freedom
c product          : product used in final Moldauer calculation
c res              : width fluctuation factor
c numtr            : number of transmission coefficients
c ielas            : designator for elastic channel
c dab,factor       : help variables
c factora,factorb,x: help variables
c
c Initialization
c
      res=0.
      dab=0.
      if (ielas.eq.1) dab=1.
      factor=1.+2.*dab/vnu(na)
      factora=real(2.*tav(na)/(st*vnu(na)))
      factorb=real(2.*tav(nb)/(st*vnu(nb)))
c
c Loop over integration points
c
      do 10 im=1,npmold
        x=xmo(im)
c
c Final result
c
        if (nb.eq.numtjl+1) then
c
c Special case for capture
c
          res=res+product(im)/(1.+x*factora)
        else
          res=res+product(im)*factor/((1.+x*factora)*(1.+x*factorb))
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
