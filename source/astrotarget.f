      subroutine astrotarget
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : December 30, 2011
c | Task  : Compound reaction for many target states
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          elas1,elas2,elastic
      integer          Zcomp,Ncomp,Zixtarget,Nixtarget,
     +                 nexastro,l2maxhf,NL,Pbeg,Pend,spin2beg,spin2end,
     +                 Ptarget,spin2target,spintarget,parity,J2,J2cnend,
     +                 J,type,
     +                 parspin2o,Zix,Nix,nexout,Pprimebeg,Pprimeend,
     +                 Irspin2beg,Irspin2end,J2res,Pprime,Irspin2,Ir,
     +                 ielas
      real             Tinc,fluxsum,Tout
      double precision Wabinelastic
      double precision rhoinc,factor1,sumIPE,sumIP,rho,compterm,rhoel
      double precision sumIPas,suminl

c
c ********************** Astrophysical situation ***********************
c
c The astrophysical calculation may involve also the excited states of
c the target located above the ground state of the nucleus or both
c below and above when dealing with a isomeric target. The temperature
c during the evolution of the universe or in stars is such that a
c nucleus can exist in various excited states.
c This would in principle imply that the loops over the quantum number
c must be performed for all possible inelastic channels. However,
c in practice these loops are reduced by a lumping approximation for
c lprime and jprime quantum numbers which means that we only loop
c over excited states number,spin and parity of the target nucleus.
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c Zixtarget: charge number index of target nucleus
c Nixtarget: neutron number index of target nucleus
c k0       : index of incident particle
c parspin  : spin of particle
c flagcheck: flag for output of numerical checks
c flagwidth: flag for width fluctuation calculation
c maxex    : maximum excitation energy bin for residual nucleus
c Ltarget  : excited level of target
c l2maxhf  : 2 * lmaxhf
c lmaxhf   : maximal l-value for transmission coefficients
c Nlast,NL : last discrete level
c parlev   : parity of level
c Pbeg,Pend: begin and end of parity summation
c jdis     : spin of level
c spin2beg : 2 * start of target spin summation
c spin2end : 2 * end of target spin summation
c
c 10: Sum over target states
c
      Zcomp=0
      Ncomp=0
      Zixtarget=Zindex(Zcomp,Ncomp,k0)
      Nixtarget=Nindex(Zcomp,Ncomp,k0)
      if (flagcheck.and.flagwidth)
     +  write(*,'(/"Flux check for astrophysical case"/)')
      do 10 nexastro=0,maxex(Zixtarget,Nixtarget)
        if (nexastro.eq.Ltarget) goto 10
        l2maxhf=2*lmaxhf(k0,nexastro)
        NL=Nlast(Zixtarget,Nixtarget,0)
        if (nexastro.le.NL) then
          Pbeg=parlev(Zixtarget,Nixtarget,nexastro)
          Pend=Pbeg
          spin2beg=int(2.*jdis(Zixtarget,Nixtarget,nexastro))
          spin2end=spin2beg
        else
          Pbeg=-1
          Pend=1
          spin2beg=mod(int(2*jdis(Zixtarget,Nixtarget,0)),2)
          spin2end=2*maxJ(Zixtarget,Nixtarget,nexastro)
        endif
c
c 20: Sum over target parity
c
c Ptarget: target parity
c
        do 20 Ptarget=Pbeg,Pend,2
c
c 30: Sum over target spin
c
c spin2target: 2 * target spin
c spintarget : target spin
c rhoinc     : integrated level density for target
c J2beg      : 2 * start of J summation
c J2cnend    : 2 * end of J summation of CN
c numJ       : maximal J-value
c
          do 30 spin2target=spin2beg,spin2end,2
            spintarget=spin2target/2
            rhoinc=rho0(k0,nexastro,spintarget,Ptarget)
            if (rhoinc.eq.0.) goto 30
            J2cnend=int(2*(lmaxhf(k0,nexastro)+parspin(k0)+spintarget))
            J2cnend=min(J2cnend,numJ)
c
c 110: Sum over Compound Nucleus parity
c
c parity: parity
c pardif: difference between target and compound nucleus parity
c
            do 110 parity=-1,1,2
              pardif=abs(Ptarget-parity)/2
c
c 120: Sum over Compound Nucleus total angular momentum J (J2)
c
c J2          : 2 * J
c flagfission : flag for fission
c nfisbar     : number of fission barrier parameters
c tfission    : subroutine for fission transmission coefficients
c astroprepare: subroutine to prepare information for astrophysical
c               compound nucleus calculations
c Tinc        : transmission coefficients as a function of j and l
c               for the incident channel
c Tastroinc   : transmission coefficient for incident channel
c               (Astrophysical case)
c tnumi       : counter for width fluctuation calculation
c denomhf     : denominator for compound nucleus formula
c widthprepare: subroutine for preparation of width fluctuation
c               corrections


              do 120 J2=J2beg,J2cnend,2
                J=J2/2
                if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0)
     +            call tfission(Zcomp,Ncomp,nexastro,J2,parity)
                call astroprepare(Zcomp,Ncomp,J2,parity,spin2target,
     +            Ptarget,nexastro)
                if (denomhf.eq.0.) goto 120
                if (flagwidth) then
                  Tinc=Tastroinc(1,nexastro,spintarget,Ptarget)
                  rhoinc=Tastroinc(0,nexastro,spintarget,Ptarget)
                  tnumi=1
                  if (Tinc.lt.transeps) goto 120
                  call widthprepare
                else
                  if (feed.lt.transeps) goto 120
                endif
c
c 210: Sum over outgoing channels
c
c fluxsum      : check for conservation of flux per P,J,j,l
c wmode        : designator for width fluctuation model
c tnumo,tnuminc: counters for width fluctuation calculation
c tnum         : total number of transmission coefficients
c widthfluc    : subroutine for width fluctuation correction
c factor1      : help variable
c tfis         : fission transmission coefficients
c xsbinary     : cross section from initial compound to residual nucleus
c CNfactor     : factor for compound nucleus cross section:
c                pi/[ k**2 (2s+1)(2I+1) ]
c
                fluxsum=0.
c
c 1. Fission channel
c
                if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0) then
                  if (flagwidth) then
                    tnumo=tnum
                    if (wmode.ge.1) then
                      call widthfluc(0)
                      factor1=real(Tinc*tfis(J,parity)/denomhf*Wab)
                      fluxsum=fluxsum+factor1
                    else
                      factor1=real(Tinc*tfis(J,parity)/denomhf)
                    endif
                  else
                    factor1=real(feed*tfis(J,parity)/denomhf)
                  endif
                  xsbinary(-1)=xsbinary(-1)+CNfactor*(J2+1.)*factor1
                endif
c
c 2. Gamma and particle channels
c
c parskip   : logical to skip outgoing particle
c parspin2o : 2 * particle spin for outgoing channel
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c
                tnumo=tnum+1
                do 210 type=0,6
                  if (type.eq.1) tnumo=tnuminc
                  if (parskip(type)) goto 210
                  parspin2o=int(2.*parspin(type))
                  Zix=Zindex(Zcomp,Ncomp,type)
                  Nix=Nindex(Zcomp,Ncomp,type)
                  NL=Nlast(Zix,Nix,0)
c
c 220: Sum over outgoing excitation energies
c
c sumIPE    : compound contribution summed over residual spin and parity
c             and energy
c elastic   : designator for elastic channel
c Pprimebeg : start of residual parity summation
c parlev    : parity of level
c Pprimeend : end of residual parity summation
c Irspin2beg: 2 * start of residual spin summation
c jdis      : spin of level
c Irspin2end: 2 * end of residual spin summation
c J2res     : help variable
c maxJ      : maximal J-value
c sumIP     : compound contribution summed over residual spin and parity
c
                  sumIPE=0.
                  do 220 nexout=0,maxex(Zix,Nix)
                    l2maxhf=2*lmaxhf(type,nexout)
                    elas1=(type.eq.k0.and.nexout.eq.nexastro)
                    if (nexout.le.NL) then
                      Pprimebeg=parlev(Zix,Nix,nexout)
                      Pprimeend=Pprimebeg
                      Irspin2beg=int(2.*jdis(Zix,Nix,nexout))
                      Irspin2end=Irspin2beg
                    else
                      Pprimebeg=-1
                      Pprimeend=1
                      J2res=J2+parspin2o
                      Irspin2beg=mod(J2res,2)
                      Irspin2end=J2res+l2maxhf
                      Irspin2end=min(Irspin2end,2*maxJ(Zix,Nix,nexout))
                    endif
                    sumIP=0.
                    sumIPas=0.
c
c 230: Sum over residual parity
c
                    do 230 Pprime=Pprimebeg,Pprimeend,2
c
c 240: Sum over residual spin
c
c Irspin2  : 2 * residual spin
c Ir       : residual spin
c rho      : integrated level density
c Tout     : transmission coefficient
c Tastroout: transmission coefficient for outgoing channel
c            (Astrophysical case)
c
                      do 240 Irspin2=Irspin2beg,Irspin2end,2
                        Ir=Irspin2/2
                        elas2=(spintarget.eq.Ir.and.Ptarget.eq.Pprime)
                        elastic=(elas1.and.elas2)
                        if (flagwidth) then
                          rho=Tastroout(0,type,nexout,Ir,Pprime)
                          if (rho.lt.1.0d-20) goto 240
                          if (rho.ne.0.)
     +                      Tout=Tastroout(1,type,nexout,Ir,Pprime)/rho
                          if (type.ge.1) tnumo=tnumo+1
                          if (wmode.ge.1) then
                            ielas=0
                            call widthfluc(ielas)
                            Wabinelastic=Wab
                            if (elastic) then
                              ielas=1
                              call widthfluc(ielas)
                              rhoel=rho
                            endif
                          endif
                          factor1=real(Tinc*rho*Tout/denomhf*Wab)
c
c We avoid double counting for the elastic channel
c take into account the rho*(rho-1) inelastic contribution of the lumped
c  (l,j) incident channels
c For more explanations, please call 911
c
c enumhf    : enumerator for compound nucleus formula
c compterm  : partial contribution to compound nucleus term
c xspop     : population cross section
c xspopex   : population cross section summed over spin and parity
c xscompcont: compound cross section for continuum
c contrib   : contribution to emission spectrum
c xspopnuc  : population cross section per nucleus
c
c
                          if (ielas.eq.1) then
                             factor1=factor1+real((rho-1.)*Tinc*rho*Tout
     +                         /denomhf*Wabinelastic)
                             if (rho.gt.0.) factor1=factor1/rho
                             suminl=real((rho-1.)*Tinc*rho*Tout
     +                         /denomhf*Wabinelastic)
                             if (rho.gt.0.) suminl=suminl/rho
                          endif
                          fluxsum=fluxsum+factor1
                        else
                          factor1=feed*
     +                      enumhf(type,nexout,Ir,Pprime)/denomhf
                          if (flagastro.and.elastic) then
                            sumIPas=sumIPas+CNfactor*(J2+1)*
     +                        enumhf(type,nexout,Ir,Pprime)**2/denomhf
                          endif
                        endif
                        compterm=CNfactor*(J2+1.)*factor1
                        xspop(Zix,Nix,nexout,Ir,Pprime)=
     +                    xspop(Zix,Nix,nexout,Ir,Pprime)+compterm
                        sumIP=sumIP+compterm
  240                 continue
  230               continue
                    xspopex(Zix,Nix,nexout)=
     +                xspopex(Zix,Nix,nexout)+sumIP
                    if (nexout.gt.NL) then
                      xscompcont(type)=xscompcont(type)+sumIP
                      contrib(type,nexout)=contrib(type,nexout)+sumIP
                    endif
c
c Compound elastic scattering is excluded from the residual production
c cross sections
c
                    if (.not.flagastro) then
                      if (.not.elastic) sumIPE=sumIPE+sumIP
                    else
                      if (flagwidth) then
                        if (.not.elastic) then
                          sumIPE=sumIPE+sumIP
                        else
                          sumIPE=sumIPE+suminl
                        endif
                      else
                        if (elastic) then
                          sumIPE=sumIPE+(sumIP-sumIPas)
                        else
                          sumIPE=sumIPE+sumIP
                        endif
                      endif
                    endif
  220             continue
                  xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)+sumIPE
                  xsbinary(type)=xsbinary(type)+sumIPE
  210           continue
c
c ****** Check of flux conservation of transmission coefficients *******
c
c cparity: parity (character)
c
c This check is included to test the stability of the width fluctuation
c calculation.
c
                if (flagcheck.and.flagwidth) then
                  if (fluxsum.eq.0.) fluxsum=Tinc
                  write(*,'(" Parity=",a1,"  J=",f4.1," Tinc",1p,e12.5,
     +              " Sum over outgoing channels=",1p,e12.5,"  Ratio=",
     +              f8.5,"  Rho=",2g14.6)') cparity(parity),0.5*J2,
     +              Tinc,fluxsum,Tinc/fluxsum,rhoinc,rhoel
                endif
  120         continue
  110       continue
   30     continue
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
