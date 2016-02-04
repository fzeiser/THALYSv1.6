      function preeqpair(Zix,Nix,n,E,pmodel)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : March 7, 2010
c | Task  : Pairing effects in exciton model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,n,pmodel
      real    preeqpair,E,pairCN,gsp,gsn,gs,pair0,Tc,ncrit,pairex
c
c **************************** Fu formula ******************************
c
c preeqpair  : pre-equilibrium pairing energy
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c n          : exciton number
c E          : excitation energy
c pmodel     : model for preequilibrium pairing energy
c pair,pairCN: pairing energy
c flag2comp  : flag for two-component pre-equilibrium model
c gp,gsp     : single-particle proton level density parameter
c gn,gsn     : single-particle neutron level density parameter
c gs         : single-particle level density parameter
c pair0      : ground-state pairing gap
c Tc         : critical temperature
c ncrit      : average quasi-particle number at the Tc
c pairex     : excited state pairing gap
c
      preeqpair=0.
      pairCN=pair(Zix,Nix)
      if (pairCN.le.0.) return
c
c pmodel 1: Fu formula
c
      if (pmodel.eq.1) then
        if (flag2comp) then
          gsp=gp(Zix,Nix)
          gsn=gn(Zix,Nix)
          gs=gsp+gsn
        else
          gs=g(Zix,Nix)
        endif
        pair0=sqrt(pairCN/(0.25*gs))
        Tc=2.*pair0/3.5
        ncrit=2.*gs*Tc*log(2.)
        if (E/pairCN.ge.(0.716+2.44*(n/ncrit)**2.17)) then
          pairex=pair0*(0.996-1.76*(n/ncrit)**1.6/((E/pairCN)**0.68))
        else
          pairex=0.
        endif
        preeqpair=pairCN-0.25*gs*(pairex**2.)
      else
c
c pmodel 2: compoud nucleus value
c
        preeqpair=pairCN
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
