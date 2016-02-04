      subroutine multipreeq2(Zcomp,Ncomp,nex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 9, 2013
c | Task  : Two-component multiple preequilibrium model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell
      integer          Zcomp,Ncomp,nex,ipp,ihp,ipn,ihn,ip,ih,p,ppi,hpi,
     +                 pnu,hnu,h,type,Zix,Nix,zejec,nejec,nexout,nen,
     +                 parity,J
      real             sumfeed,gsp,gsn,damp,ignatyuk,omegaph,phdens2,
     +                 dEx,Eex,Eo,term(2,0:numex),omegap1h,
     +                 EoplusS,omega1p,proba,Tswave,Pescape,Exm,Exmin,
     +                 factor(2,0:numex,0:numparx,0:numparx),Jterm
      double precision summpe,feedph,sumterm,sumph,sumtype(2)
c
c mpreeqmode=1: Multiple emission exciton model
c mpreeqmode=2: Transmission coefficients for particle emission
c
c Multiple preequilibrium emission model 2 is adopted from Chadwick
c and Young, Phys. Rev. C50 (1994) p. 996.
c
c ************ Check presence of multiple pre-equilibrium **************
c
c Zcomp   : charge number index for compound nucleus
c Ncomp   : neutron number index for compound nucleus
c sumfeed : help variable
c ipp     : proton particle number of mother nucleus
c maxpar  : maximal particle number
c ihp     : proton hole number of mother nucleus
c ipn     : neutron particle number of mother nucleus
c ihn     : neutron hole number of mother nucleus
c xspopph2: population cross section per two-component particle-hole
c           configuration
c
c There must be excited particles and holes present for multiple
c pre-equilibrium emission to occur.
c
      if (Zcomp.eq.0.and.Ncomp.eq.0) return
      sumfeed=0.
      do ipn=0,numparx
        do ipp=0,numparx
          do nexout=0,numex
            do type=1,2
              factor(type,nexout,ipn,ipp)=0.
            enddo
          enddo
        enddo
      enddo
      do 10 ipp=0,maxpar
        do 10 ihp=0,maxpar
          do 10 ipn=0,maxpar
            do 10 ihn=0,maxpar
              sumfeed=sumfeed+xspopph2(Zcomp,Ncomp,nex,ipp,ihp,ipn,ihn)
   10 continue
      if (sumfeed.le.1.e-10) return
c
c ***************** Multiple preequilibrium emission *******************
c
c surfwell  : flag for surface effects in finite well
c Ecomp     : total energy of composite system
c Exinc     : excitation energy of entrance bin
c gp,gsp    : single-particle proton level density parameter
c gn,gsn    : single-particle neutron level density parameter
c flaggshell: flag for energy dependence of single particle level
c             density parameter g
c damp      : shell damping factor
c ignatyuk  : function for energy dependent level density parameter a
c alev      : level density parameter
c summpe    : multiple preequilibrium emission for excitation energy bin
c
      surfwell=.false.
      Ecomp=Exinc
      gsp=gp(Zcomp,Ncomp)
      gsn=gn(Zcomp,Ncomp)
      if (flaggshell) then
        damp=ignatyuk(Zcomp,Ncomp,Ecomp,0)/alev(Zcomp,Ncomp)
        gsp=gsp*damp
        gsn=gsn*damp
      endif
      summpe=0.
c
c 110: Loops over all possible particle-hole excitations of the mother
c      bin.
c
c ip        : particle number
c ih        : hole number
c feedph    : feeding term from previous particle-hole calculation
c sumph     : multiple preequilibrium flux for particle-hole pair
c mpreeqmode: designator for multiple pre-equilibrium model
c omegaph   : particle-hole state density for compound system
c phdens2   : function for two-component particle-hole state density
c Efermi    : depth of Fermi well
c p0        : initial particle number
c ppi0      : initial proton number
c hpi0      : initial proton hole number
c pnu0      : initial neutron number
c hnu0      : initial neutron hole number
c exchange2 : subroutine for calculation of two-component exchange terms
c p         : particle number
c ppi       : proton particle number
c hpi       : proton hole number
c pnu       : neutron particle number
c hnu       : neutron hole number
c h         : hole number
c lifetime2 : subroutine for calculation of lifetime of two-component
c             exciton state
c
      do 110 ipp=0,maxpar
        do 110 ihp=0,maxpar
          do 110 ipn=0,maxpar
            ip=ipp+ipn
            if (ip.eq.0.or.ip.gt.maxpar) goto 110
            do 120 ihn=0,maxpar
              ih=ihp+ihn
              if (ih.eq.0.or.ih.gt.maxpar) goto 120
              feedph=xspopph2(Zcomp,Ncomp,nex,ipp,ihp,ipn,ihn)
              if (feedph.le.1.e-10) goto 120
              sumph=0.
              if (mpreeqmode.eq.2) then
                omegaph=phdens2(Zcomp,Ncomp,ipp,ihp,ipn,ihn,gsp,gsn,
     +            Exinc,Efermi,surfwell)
              else
                p0=ip
                ppi0=ipp
                hpi0=ihp
                pnu0=ipn
                hnu0=ihn
                call exchange2(Zcomp,Ncomp)
                do 130 p=p0,maxpar
                  do 130 ppi=ppi0,maxpar
                    hpi=hpi0+ppi-ppi0
                    do 130 pnu=pnu0,maxpar
                      hnu=hnu0+pnu-pnu0
                      h=hpi+hnu
                      if (ppi+pnu.eq.p.and.h.le.maxpar)
     +                  call lifetime2(ppi,hpi,pnu,hnu)
  130           continue
              endif
c
c 140: Loop over neutron and proton type of next emitted particle
c
c parskip    : logical to skip outgoing particle
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c numZph     : maximal number of protons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c numNph     : maximal number of neutrons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c zejec,parZ : charge number of leading particle
c nejec,parN : neutron number of leading particle
c sumtype    : help variable
c deltaEx,dEx: excitation energy bin for population arrays
c
          sumterm=0.
          do 140 type=1,2
            sumtype(type)=0.
            do nexout=0,numex
              term(type,nexout)=0.
            enddo
            if (parskip(type)) goto 140
            Zix=Zindex(Zcomp,Ncomp,type)
            Nix=Nindex(Zcomp,Ncomp,type)
            if (Zix.gt.numZph.or.Nix.gt.numNph) goto 140
            zejec=parZ(type)
            nejec=parN(type)
            if (mpreeqmode.eq.2) then
              if (ipp-zejec.lt.0) goto 140
              if (ipn-nejec.lt.0) goto 140
            endif
c
c 150: Loop over residual excitation energy bins
c
c Nlast     : last discrete level
c nexmax    : maximum excitation energy bin for residual nucleus
c Ex,Eex    : excitation energy
c Eo        : outgoing energy
c S         : separation energy per particle
c locate    : subroutine to find value in ordered table
c egrid     : energies of basic energy grid in MeV
c ebegin    : first energy point of energy grid
c eend      : last energy point for energy grid of each particle
c term      : help variables
c omegap1h  : particle-hole state density for residual system
c EoplusS   : outgoing energy+S
c omega1p   : particle-hole state density for continuum particle
c proba     : probability to find continuum particle
c Tswave    : transmission coefficient for s-wave
c Tjl       : transmission coefficients as a function of particle
c             type, energy, spin and l-value
c Pescape   : escape probability
c Exm       : maximal attainable energy
c dExinc    : excitation energy bin for mother nucleus
c Exmin     : lower bound of energy bin
c factor    : help variable
c Spre      : time-integrated strength of two-component exciton state
c wemission : two-component emission rate
c mcontrib  : contribution to emission spectrum
c mpecontrib: contribution to multiple pre-equilibrium emission spectrum
c xspopex   : population cross section summed over spin and parity
c
            do 150 nexout=Nlast(Zix,Nix,0)+1,nexmax(type)
              dEx=deltaEx(Zix,Nix,nexout)
              Eex=Ex(Zix,Nix,nexout)
              Eo=Exinc-Eex-S(Zcomp,Ncomp,type)
              call locate(egrid,ebegin(type),eend(type),Eo,nen)
              if (mpreeqmode.eq.2) then
                if (omegaph.gt.0.) then
                  gsp=gp(Zix,Nix)
                  gsn=gn(Zix,Nix)
                  if (flaggshell) then
                    damp=ignatyuk(Zix,Nix,Eex,0)/alev(Zix,Nix)
                    gsp=gsp*damp
                    gsn=gsn*damp
                  endif
                  omegap1h=phdens2(Zix,Nix,ipp-zejec,ihp,ipn-nejec,ihn,
     +              gsp,gsn,Eex,Efermi,surfwell)
                  EoplusS=Exinc-Eex
                  omega1p=phdens2(Zix,Nix,zejec,0,nejec,0,gsp,gsn,
     +              EoplusS,Efermi,surfwell)
                  proba=omega1p*omegap1h/omegaph/(ipp+ipn)
                  Tswave=Tjl(type,nen,1,0)
                  Pescape=proba*Tswave
                  if (nexout.eq.nexmax(type)) then
                    Exm=Exinc+0.5*dExinc-S(Zcomp,Ncomp,type)
                    Exmin=Ex(Zix,Nix,nexout)-0.5*dEx
                    dEx=Exm-Exmin
                  endif
                  term(type,nexout)=feedph*Pescape*dEx
                  sumterm=sumterm+term(type,nexout)
                endif
              else
                do 160 ppi=ppi0,maxpar
                  if (ppi-zejec.lt.0) goto 160
                  hpi=hpi0+ppi-ppi0
                  do 170 pnu=pnu0,maxpar
                    hnu=hnu0+pnu-pnu0
                    if (pnu-nejec.lt.0) goto 170
                    h=hpi+hnu
                    if (h.gt.maxpar) goto 170
                    factor(type,nexout,ppi,pnu)=feedph*
     +                Spre(ppi,hpi,pnu,hnu)*
     +                wemission2(type,ppi,hpi,pnu,hnu,nen)*dEx
                    term(type,nexout)=term(type,nexout)+
     +                factor(type,nexout,ppi,pnu)
  170             continue
  160           continue
                sumterm=sumterm+term(type,nexout)
              endif
              sumtype(type)=sumtype(type)+term(type,nexout)
  150       continue
  140     continue
c
c Normalization
c
          if (sumterm.gt.feedph) then
            do 210 type=1,2
              do 220 nexout=Nlast(Zix,Nix,0)+1,nexmax(type)
                term(type,nexout)=term(type,nexout)*feedph/sumterm
                if (mpreeqmode.eq.1) then
                  do 230 ppi=ppi0,maxpar
                    do 240 pnu=pnu0,maxpar
                      factor(type,nexout,ppi,pnu)=
     +                  factor(type,nexout,ppi,pnu)*feedph/sumterm
  240               continue
  230             continue
                endif
  220         continue
              sumtype(type)=sumtype(type)*feedph/sumterm
  210       continue
          endif
c
c Feed new population bins
c
c preeqpopex: pre-equilibrium population cross section summed over
c             spin and parity
c
          do 310 type=1,2
            if (parskip(type)) goto 310
            Zix=Zindex(Zcomp,Ncomp,type)
            Nix=Nindex(Zcomp,Ncomp,type)
            if (Zix.gt.numZph.or.Nix.gt.numNph) goto 310
            zejec=parZ(type)
            nejec=parN(type)
            if (mpreeqmode.eq.2) then
              if (ipp-zejec.lt.0) goto 310
              if (ipn-nejec.lt.0) goto 310
            endif
            do 320 nexout=Nlast(Zix,Nix,0)+1,nexmax(type)
              if (mpreeqmode.eq.2) then
                xspopph2(Zix,Nix,nexout,ipp-zejec,ihp,ipn-nejec,ihn)=
     +            xspopph2(Zix,Nix,nexout,ipp-zejec,ihp,ipn-nejec,ihn)
     +            +term(type,nexout)
              else
                do 330 ppi=ppi0,maxpar
                  if (ppi-zejec.lt.0) goto 330
                  hpi=hpi0+ppi-ppi0
                  do 340 pnu=pnu0,maxpar
                    hnu=hnu0+pnu-pnu0
                    if (pnu-nejec.lt.0) goto 340
                    h=hpi+hnu
                    if (h.gt.maxpar) goto 340
                    xspopph2(Zix,Nix,nexout,ppi-zejec,hpi,pnu-nejec,hnu)
     +                =xspopph2(Zix,Nix,nexout,ppi-zejec,hpi,pnu-nejec,
     +                hnu)+factor(type,nexout,ppi,pnu)
  340             continue
  330           continue
              endif
              mcontrib(type,nex,nexout)=mcontrib(type,nex,nexout)+
     +          term(type,nexout)
              mpecontrib(type,nex,nexout)=mpecontrib(type,nex,nexout)+
     +          term(type,nexout)
              xspopex(Zix,Nix,nexout)=xspopex(Zix,Nix,nexout)+
     +          term(type,nexout)
              preeqpopex(Zix,Nix,nexout)=preeqpopex(Zix,Nix,nexout)+
     +          term(type,nexout)
c
c 350: Adopt spin distribution from primary pre-equilibrium.
c
c parity  : parity
c maxJ    : maximal J-value
c xspop   : population cross section
c preeqpop: pre-equilibrium population cross section
c RnJ     : spin distribution for particle-hole states
c RnJsum  : (2J+1)*sum over spin distributions
c
              do 350 parity=-1,1,2
                do 350 J=0,maxJ(Zix,Nix,nexout)
                  Jterm=0.5*(2*J+1)*RnJ(2,J)/RnJsum(2)*term(type,nexout)
                  xspop(Zix,Nix,nexout,J,parity)=
     +              xspop(Zix,Nix,nexout,J,parity)+Jterm
                  preeqpop(Zix,Nix,nexout,J,parity)=
     +              preeqpop(Zix,Nix,nexout,J,parity)+Jterm
  350         continue
  320       continue
c
c Add total multiple pre-equilibrium contributions
c
c xspopnuc  : population cross section per nucleus
c xspoppreeq: preequilibrium population cross section per nucleus
c xspartial : emitted cross section flux per energy bin
c xsfeed    : cross section from compound to residual nucleus
c xsmpe     : multiple-preequilibrium cross section per energy bin
c xsmpetot  : total multiple-preequilibrium cross section
c mulpreZN  : logical for multiple pre-equilibrium per nucleus
c
            sumph=sumph+sumtype(type)
            summpe=summpe+sumtype(type)
            xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)+sumtype(type)
            xspoppreeq(Zix,Nix)=xspoppreeq(Zix,Nix)+sumtype(type)
            xspartial(type,nex)=xspartial(type,nex)+sumtype(type)
            xsfeed(Zcomp,Ncomp,type)=xsfeed(Zcomp,Ncomp,type)+
     +        sumtype(type)
            xsmpe(type,nex)=xsmpe(type,nex)+sumtype(type)
            xsmpetot(type)=xsmpetot(type)+sumtype(type)
            if (sumtype(type).ne.0.) mulpreZN(Zix,Nix)=.true.
  310     continue
c
c Flux that is not emitted during a particular stage, is always
c transferred to the next stage.
c
          if (mpreeqmode.eq.2) then
            if (ip.le.maxpar-1.and.ih.le.maxpar-1) then
              if (ipp.le.maxpar-1.and.ihp.le.maxpar-1)
     +          xspopph2(Zcomp,Ncomp,nex,ipp+1,ihp+1,ipn,ihn)=
     +          xspopph2(Zcomp,Ncomp,nex,ipp+1,ihp+1,ipn,ihn)+
     +          0.5*(feedph-sumph)
              if (ipn.le.maxpar-1.and.ihn.le.maxpar-1)
     +          xspopph2(Zcomp,Ncomp,nex,ipp,ihp,ipn+1,ihn+1)=
     +          xspopph2(Zcomp,Ncomp,nex,ipp,ihp,ipn+1,ihn+1)+
     +          0.5*(feedph-sumph)
            endif
          endif
  120   continue
  110 continue
c
c ************************ Normalization *******************************
c
c Dmulti: depletion factor for multiple preequilibrium
c
      Dmulti(nex)=summpe/xspopex(Zcomp,Ncomp,nex)
      xspopex(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)-summpe
      preeqpopex(Zcomp,Ncomp,nex)=preeqpopex(Zcomp,Ncomp,nex)-summpe
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
