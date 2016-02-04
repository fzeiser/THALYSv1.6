      subroutine lifetime2(ppi,hpi,pnu,hnu)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : August 9, 2004
c | Task  : Calculation of lifetime of two-component exciton state
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ppi,hpi,pnu,hnu,p
      real    term1,term2,term3,term4,term5,term6,term7,term8
c
c **************** Calculation of total strength ***********************
c
c ppi         : proton particle number
c hpi         : proton hole number
c pnu         : neutron particle number
c hnu         : neutron hole number
c p           : particle number
c PP2         : total strength
c ppi0        : initial proton number
c hpi0        : initial proton hole number
c pnu0        : initial neutron number
c hnu0        : initial neutron hole number
c term1,......: help variables
c Gpiplus,... : two-component branching ratios
c maxpar      : maximal particle number
c Lexc        : exchange term
c Spre        : time-integrated strength of two-component exciton state
c tauexc2     : lifetime of two-component exciton state
c
c The strength of the exciton state is calculated
c
      p=ppi+pnu
      if (p.eq.p0) then
        PP2(ppi0,hpi0,pnu0,hnu0)=1.
      else
        term1=0.
        term2=0.
        term3=0.
        term4=0.
        term5=0.
        term6=0.
        term7=0.
        term8=0.
        if (ppi.gt.ppi0.and.hpi.gt.hpi0) then
          term1=PP2(ppi-1,hpi-1,pnu,hnu)*Gpiplus(ppi-1,hpi-1,pnu,hnu)
          term4=PP2(ppi-1,hpi-1,pnu,hnu)*Gnuplus(ppi-1,hpi-1,pnu,hnu)
        endif
        if (pnu.gt.pnu0.and.hnu.gt.hnu0) then
          term2=PP2(ppi,hpi,pnu-1,hnu-1)*Gnuplus(ppi,hpi,pnu-1,hnu-1)
          term6=PP2(ppi,hpi,pnu-1,hnu-1)*Gpiplus(ppi,hpi,pnu-1,hnu-1)
        endif
        if (pnu.lt.maxpar.and.hnu.lt.maxpar)  then
          if (ppi.gt.ppi0.and.hpi.gt.hpi0)
     +      term5=Gnupi(ppi-1,hpi-1,pnu+1,hnu+1)
          if (ppi-1.gt.ppi0.and.hpi-1.gt.hpi0)
     +      term3=PP2(ppi-2,hpi-2,pnu+1,hnu+1)*
     +      Gpiplus(ppi-2,hpi-2,pnu+1,hnu+1)
        endif
        if (ppi.lt.maxpar.and.hpi.lt.maxpar) then
          if (pnu.gt.pnu0.and.hnu.gt.hnu0)
     +      term8=Gpinu(ppi+1,hpi+1,pnu-1,hnu-1)
          if (pnu-1.gt.pnu0.and.hnu-1.ge.hnu0)
     +      term7=PP2(ppi+1,hpi+1,pnu-2,hnu-2)*
     +      Gnuplus(ppi+1,hpi+1,pnu-2,hnu-2)
        endif
        PP2(ppi,hpi,pnu,hnu)=term1+term2+Lexc(ppi,hpi,pnu,hnu)*
     +    ((term3+term4)*term5+(term6+term7)*term8)
      endif
      Spre(ppi,hpi,pnu,hnu)=PP2(ppi,hpi,pnu,hnu)*
     +  tauexc2(ppi,hpi,pnu,hnu)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
