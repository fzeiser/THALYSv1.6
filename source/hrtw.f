      subroutine hrtw(tjl,na,nb,st,sv,v,w,res,numtr,ielas)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : February 6, 2007
c | Task  : HRTW width fluctuation correction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          na,nb,numtr,ielas
      real             sv,v(numtr),w(numtr),res,dab
      double precision tjl(0:5,numtr),st,tt,res1
c
c ************************** HRTW calculation **************************
c
c We use the model of HRTW (1975), revised in 1980.
c
c tjl             : transmission coefficients
c na,nb           : counters for width fluctuation calculation
c st              : denominator of compound nucleus formula
c res             : width fluctuation factor
c numtr           : number of transmission coefficients
c ielas           : designator for elastic channel
c dab,tt,sv,w,res1: help variables
c
      dab=0.
      if (ielas.eq.1) dab=1.
c
c Final result
c
      tt=tjl(1,na)*tjl(1,nb)/tjl(0,na)
      if (tt.ne.0) then
        res1=v(na)*v(nb)*tjl(0,nb)/sv*(1.+dab*(w(na)-1.))
        res=real(res1*st/tt)
      else
        res=0.
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
