      subroutine exciton2out
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 5, 2006
c | Task  : Output of two-component exciton model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer p,ppi,hpi,pnu,hnu,h,n,type
      real    lambdapiplus,lambdanuplus,lambdapinu,lambdanupi
c
c ************************ Exciton model *******************************
c
      write(*,'(/" ++++++++++ TWO-COMPONENT EXCITON MODEL ++++++++++")')
c
c 1. Output of matrix element
c
c Ecomp     : total energy of composite system
c M2constant: constant for matrix element in exciton model
c Rpinu,....: ratio for two-component matrix element
c p         : particle number
c p0        : initial particle number
c maxpar    : maximal particle number
c ppi       : proton particle number
c ppi0      : initial proton number
c hpi       : proton hole number
c pnu       : neutron particle number
c pnu0      : initial neutron number
c hnu       : neutron hole number
c h         : hole number
c n         : exciton number
c matrix    : subroutine for matrix element for exciton model
c Ainit     : mass number of initial compound nucleus
c M2pipi    : square of proton-proton matrix element
c M2nunu    : square of neutron-neutron matrix element
c M2pinu    : square of proton-neutron matrix element
c M2nupi    : square of neutron-proton matrix element
c
      write(*,'(/" 1. Matrix element for E= ",f7.3/)') Ecomp
      write(*,'(" Constant for matrix element : ",f7.3)') M2constant
      write(*,'(" p-p ratio for matrix element: ",f7.3)') Rpipi
      write(*,'(" n-n ratio for matrix element: ",f7.3)') Rnunu
      write(*,'(" p-n ratio for matrix element: ",f7.3)') Rpinu
      write(*,'(" n-p ratio for matrix element: ",f7.3/)') Rnupi
      write(*,'(" p(p) h(p) p(n) h(n)     M2pipi      M2nunu",
     +  "      M2pinu      M2nupi"/)')
      do 10 p=p0,maxpar
        do 10 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 10 pnu=pnu0,maxpar
            if (ppi+pnu.eq.p) then
              hnu=pnu-pnu0
              h=hpi+hnu
              n=p+h
              call matrix(Ainit,n)
              write(*,'(1x,4(i2,3x),1p,4e12.5)') ppi,hpi,pnu,hnu,M2pipi,
     +          M2nunu,M2pinu,M2nupi
            endif
   10 continue
c
c 2. Output of emission rates or escape widths
c
c parname   : name of particle
c wemispart2: two-component emission rate per particle and exciton
c wemistot2 : total two-component emission rate per exciton number
c hbar      : Planck's constant / 2.pi in MeV.s
c
      write(*,'(/" 2. Emission rates or escape widths"/)')
      write(*,'(" A. Emission rates ( /sec)"/)')
      write(*,'(" p(p) h(p) p(n) h(n)",4x,7(a8,4x),"Total"/)')
     +  (parname(type),type=0,6)
      do 20 p=p0,maxpar
        do 20 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 20 pnu=pnu0,maxpar
            if (ppi+pnu.eq.p) then
              hnu=pnu-pnu0
              write(*,'(1x,4(i2,3x),1p,8e12.5)') ppi,hpi,pnu,hnu,
     +        (wemispart2(type,ppi,hpi,pnu,hnu),type=0,6),
     +        wemistot2(ppi,hpi,pnu,hnu)
            endif
   20 continue
      write(*,'(/" B. Escape widths (MeV)"/)')
      write(*,'(" p(p) h(p) p(n) h(n)",4x,7(a8,4x),"Total"/)')
     +  (parname(type),type=0,6)
      do 30 p=p0,maxpar
        do 30 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 30 pnu=pnu0,maxpar
            if (ppi+pnu.eq.p) then
              hnu=pnu-pnu0
              write(*,'(1x,4(i2,3x),1p,8e12.5)') ppi,hpi,pnu,hnu,
     +        (wemispart2(type,ppi,hpi,pnu,hnu)*hbar,type=0,6),
     +        wemistot2(ppi,hpi,pnu,hnu)*hbar
            endif
   30 continue
c
c 3. Output of transition rates or damping widths and total widths
c
c lambdapiplus: proton transition rate for n --> n+2
c lambdanuplus: neutron transition rate for n --> n+2
c lambdapinu  : proton-neutron transition rate for n --> n
c lambdanupi  : neutron-proton transition rate for n --> n
c
      write(*,'(/" 3. Internal transition rates or damping widths,",
     +  " total widths"/)')
      write(*,'(" A. Internal transition rates ( /sec)"/)')
      write(*,'(" p(p) h(p) p(n) h(n)     lambdapiplus   ",
     +  "lambdanuplus    lambdapinu     lambdanupi"/)')
      do 40 p=p0,maxpar
        do 40 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 40 pnu=pnu0,maxpar
            if (ppi+pnu.eq.p) then
              hnu=pnu-pnu0
              write(*,'(1x,4(i2,3x),1p,4e15.5)') ppi,hpi,pnu,hnu,
     +          lambdapiplus(0,0,ppi,hpi,pnu,hnu),
     +          lambdanuplus(0,0,ppi,hpi,pnu,hnu),
     +          lambdapinu(0,0,ppi,hpi,pnu,hnu),
     +          lambdanupi(0,0,ppi,hpi,pnu,hnu)
            endif
   40 continue
      write(*,'(/" B. Damping widths (MeV)"/)')
      write(*,'(" p(p) h(p) p(n) h(n)     gammapiplus    ",
     +  "gammanuplus    gammapinu      gammanupi"/)')
      do 50 p=p0,maxpar
        do 50 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 50 pnu=pnu0,maxpar
            if (ppi+pnu.eq.p) then
              hnu=pnu-pnu0
              write(*,'(1x,4(i2,3x),1p,4e15.5)') ppi,hpi,pnu,hnu,
     +          lambdapiplus(0,0,ppi,hpi,pnu,hnu)*hbar,
     +          lambdanuplus(0,0,ppi,hpi,pnu,hnu)*hbar,
     +          lambdapinu(0,0,ppi,hpi,pnu,hnu)*hbar,
     +          lambdanupi(0,0,ppi,hpi,pnu,hnu)*hbar
            endif
   50 continue
      write(*,'(/" C. Total widths (MeV)"/)')
      write(*,'(" p(p) h(p) p(n) h(n)      gammatot"/)')
      do 60 p=p0,maxpar
        do 60 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 60 pnu=pnu0,maxpar
            if (ppi+pnu.eq.p) then
              hnu=pnu-pnu0
              write(*,'(1x,4(i2,3x),1p,e15.5)') ppi,hpi,pnu,hnu,
     +          hbar*(lambdapiplus(0,0,ppi,hpi,pnu,hnu)+
     +          lambdanuplus(0,0,ppi,hpi,pnu,hnu)+
     +          lambdapinu(0,0,ppi,hpi,pnu,hnu)+
     +          lambdanupi(0,0,ppi,hpi,pnu,hnu)+
     +          wemistot2(ppi,hpi,pnu,hnu))
            endif
   60 continue
c
c 4. Output of lifetimes of exciton states
c
c Spre: time-integrated strength of two-component exciton state
c
      write(*,'(/" 4. Lifetimes")')
      write(*,'(" p(p) h(p) p(n) h(n)      Strength"/)')
      do 70 p=p0,maxpar
        do 70 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 70 pnu=pnu0,maxpar
            if (ppi+pnu.eq.p) then
              hnu=pnu-pnu0
              write(*,'(1x,4(i2,3x),1p,e15.5)') ppi,hpi,pnu,hnu,
     +          Spre(ppi,hpi,pnu,hnu)
            endif
   70 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
