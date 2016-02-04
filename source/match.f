      function match(Eex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Matching function
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          i
      real             match,Eex,dEx,temp,logrhof
      double precision rhof,factor1,factor2,term
c
c *********************** Matching function ****************************
c
c match       : matching function
c Eex         : excitation energy
c pol1        : subroutine for interpolation of first order
c temp,temprho: nuclear temperature
c rhof        : value for total level density
c E0save      : E0 value saved for matching routine
c NLo,EL,NP,EP: matching level numbers and energies
c factor1,2   : help variables
c
      match=0.
      dEx=0.1
      i=max(int(Eex/dEx),1)
      call pol1(i*dEx,(i+1)*dEx,temprho(i),temprho(i+1),Eex,temp)
      if (temp.gt.0.) then
        call pol1(i*dEx,(i+1)*dEx,logrho(i),logrho(i+1),Eex,logrhof)
        rhof=exp(dble(logrhof))
        if (E0save.eq.1.e-20) then
          factor1=exp(-Eex/temp)
          factor2=exp(dble(EP/temp))
          if (EL.ne.0.) factor2=factor2-exp(EL/temp)
          term=real(min(temp*rhof*factor1*factor2,1.d30))
          match=term+NLo-NP
        else
          match=Eex-temp*log(temp*rhof)-E0save
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
