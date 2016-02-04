      subroutine exchange2(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 9, 2004
c | Task  : Calculation of two-component exchange terms
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,ppi,hpi,pnu,hnu,p,h
      real    tpiplus,lambdapiplus,tnuplus,lambdanuplus,tpinu,
     +        lambdapinu,tnupi,lambdanupi,tplus,texchange,wemis,
     +        tauexc2p,tauinv,tauinvp
c
c ************************* Exchange terms *****************************
c
c Zcomp          : charge number index for compound nucleus
c Ncomp          : neutron number index for compound nucleus
c ppi0           : initial proton number
c pnu0           : initial neutron number
c ppi            : proton particle number
c maxpar         : maximal particle number
c hpi            : proton hole number
c hpi0           : initial proton hole number
c pnu            : neutron particle number
c hnu            : neutron hole number
c hnu0           : initial neutron hole number
c p              : particle number
c h              : hole number
c tpiplus,....   : help variables
c lambdapiplus   : proton transition rate for n --> n+2
c lambdanuplus   : neutron transition rate for n --> n+2
c lambdapinu     : proton-neutron transition rate for n --> n
c lambdanupi     : neutron-proton transition rate for n --> n
c n              : exciton number
c emissionrate2  : subroutine for two-component emission rate
c wemistot2,wemis: total two-component emission rate per exciton number
c tauexc2        : lifetime of two-component exciton state
c tauinvp        : help variable
c tauexc2p       : lifetime for creation and emission only
c Lexc           : exchange term
c Gpiplus,....   : two-component branching ratios
c
c The strength of the exciton state is calculated
c
      do 10 ppi=ppi0,maxpar
        hpi=hpi0+ppi-ppi0
        do 20 pnu=pnu0,maxpar
          hnu=hnu0+pnu-pnu0
          p=ppi+pnu
          h=hpi+hnu
          if (p.gt.maxpar.or.h.gt.maxpar) goto 20
          tpiplus=lambdapiplus(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
          tnuplus=lambdanuplus(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
          tpinu=lambdapinu(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
          tnupi=lambdanupi(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
          tplus=tpiplus+tnuplus
          texchange=tpinu+tnupi
          call emissionrate2(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
          wemis=wemistot2(ppi,hpi,pnu,hnu)
          tauinv=tplus+texchange+wemis
          if (tplus.eq.0.) then
            tauexc2(ppi,hpi,pnu,hnu)=0.
            Lexc(ppi,hpi,pnu,hnu)=0.
          else
            tauexc2(ppi,hpi,pnu,hnu)=1./tauinv
            tauinvp=tplus+wemis
            if (tauinvp.eq.0.) then
              tauexc2p=0.
            else
              tauexc2p=1./tauinvp
            endif
            Lexc(ppi,hpi,pnu,hnu)=tauexc2p/tauexc2(ppi,hpi,pnu,hnu)
          endif
          Gpiplus(ppi,hpi,pnu,hnu)=tpiplus*tauexc2(ppi,hpi,pnu,hnu)
          Gnuplus(ppi,hpi,pnu,hnu)=tnuplus*tauexc2(ppi,hpi,pnu,hnu)
          Gpinu(ppi,hpi,pnu,hnu)=tpinu*tauexc2(ppi,hpi,pnu,hnu)
          Gnupi(ppi,hpi,pnu,hnu)=tnupi*tauexc2(ppi,hpi,pnu,hnu)
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
