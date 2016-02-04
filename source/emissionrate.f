      subroutine emissionrate(Zcomp,Ncomp,p,h)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 14, 2010
c | Task  : Emission rates for exciton model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell,surfgam
      integer Zcomp,Ncomp,p,h,n,type,Zix,Nix,aejec,nen,nen1,pres,nres
      real    gs,gsg,ignatyuk,edepth,U,preeqpair,phcomp,phcompg,
     +        wemissum(0:numparx,0:numparx,0:numen),Qfac,Ewell,phdens,
     +        Eout,xs,factor,Eres,phres1,phres2,g2E,branchplus,
     +        branchzero,phratio,Ures,phres,Emax,dE
c
c *************************** Emission rates ***************************
c
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c p           : particle number
c h           : hole number
c n           : exciton number
c g,gs,gsg    : single-particle level density parameter
c Atarget     : mass number of target nucleus
c flaggshell  : flag for energy dependence of single particle level
c               density parameter g
c ignatyuk    : function for energy dependent level density parameter a
c Ecomp       : total energy of composite system
c alev        : level density parameter
c surfwell    : flag for surface effects in finite well
c flagsurface : flag for surface effects in exciton model
c primary     : flag to designate primary (binary) reaction
c edepth      : depth of potential well
c Esurf       : well depth for surface interaction
c Efermi      : depth of Fermi well
c U,Ures      : excitation energy minus pairing energy
c preeqpair   : pre-equilibrium pairing energy
c pairmodel   : model for preequilibrium pairing energy
c phcomp      : particle-hole state density for compound system
c phdens      : function for particle-hole state density
c phcompg     : particle-hole state density for compound system (gamma)
c wemistot    : total emission rate per exciton number
c wemispart   : emission rate per particle and exciton number
c wemission   : emission rate per particle, exciton number and energy
c wemissum    : emission rate per exciton number and energy
c parskip     : logical to skip outgoing particle
c Zindex,Zix  : charge number index for residual nucleus
c Nindex,Nix  : neutron number index for residual nucleus
c aejec,parA  : mass number of leading particle
c Qfactor,Qfac: Q-factor for neutron/proton distinction
c Ewell       : depth of potential well
c ebegin      : first energy point of energy grid
c eend        : last energy point of energy grid
c egrid,Eout  : energies of basic energy grid in MeV
c factor,xs   : help variables
c wfac        : factor for emission rate
c xsreac      : reaction cross section
c Eres        : total energy of residual system
c S           : separation energy per particle
c
c The emission rates are derived from detailed balance, see
c the manual.
c
      n=p+h
      gs=g(Zcomp,Ncomp)
      gsg=Atarget/13.
      if (flaggshell) gs=gs*ignatyuk(Zcomp,Ncomp,Ecomp,0)/
     +  alev(Zcomp,Ncomp)
      surfwell=flagsurface.and.h.eq.1.and.primary
      if (surfwell) then
        edepth=Esurf
      else
        edepth=Efermi
      endif
      U=Ecomp-preeqpair(Zcomp,Ncomp,n,Ecomp,pairmodel)
      phcomp=phdens(Zcomp,Ncomp,p,h,gs,U,edepth,surfwell)
      phcompg=phdens(Zcomp,Ncomp,p,h,gsg,U,edepth,surfwell)
      wemistot(p,h)=0.
      do 10 nen=0,numen
        wemissum(p,h,nen)=0.
   10 continue
      do 20 type=0,6
        wemispart(type,p,h)=0.
        do 30 nen=0,numen
          wemission(type,p,h,nen)=0.
   30   continue
        if (parskip(type)) goto 20
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        gs=g(Zix,Nix)
        aejec=parA(type)
        if (primary.and.type.gt.0) then
          Qfac=Qfactor(type,p)
        else
          Qfac=1.
        endif
        if (type.gt.2) then
          Ewell=Efermi
        else
          Ewell=edepth
        endif
        if (phcomp.le.1.e-10) goto 20
        if (type.eq.0.and.phcompg.le.1.e-10) goto 20
        do 40 nen=ebegin(type),eend(type)
          Eout=egrid(nen)
          xs=xsreac(type,nen)
          factor=wfac(type)*xs*Eout
          Eres=Ecomp-S(Zcomp,Ncomp,type)-Eout
c
c Avoid anomalies at the drip line
c
          if (type.le.2) Eres=min(Eres,Ecomp)
c
c Check if outgoing energy exceeds maximal possible energy
c
          if (Eres.lt.0.) then
            nen1=nen-1
            goto 50
          endif
          if (flaggshell) gs=g(Zix,Nix)*ignatyuk(Zix,Nix,Eres,0)/
     +      alev(Zix,Nix)
c
c 1. Gamma emission rates
c
c surfgam   : flag for surface effects for photons (always false)
c phres1-2  : particle-hole state density for residual system
c g2E       : help variable
c branchplus: branching ratio for n-2 --> n
c branchzero: branching ratio for n --> n
c phratio   : ratio between residual and compound particle-hole density
c Rgamma    : adjustable parameter for pre-equilibrium gamma decay
c
c Gamma emission is only included for primary pre-equilibrium decay
c and n <= 7.
c
          if (type.eq.0) then
            if (primary.and.n.le.7) then
              factor=factor*Eout
              U=max(Eres-preeqpair(Zcomp,Ncomp,n,Eres,pairmodel),
     +          preeqpair(Zcomp,Ncomp,n,Eres,pairmodel))
              surfgam=.false.
              phres1=phdens(Zcomp,Ncomp,p-1,h-1,gsg,U,Efermi,surfgam)
              phres2=phdens(Zcomp,Ncomp,p,h,gsg,U,Efermi,surfgam)
              g2E=gsg*gsg*Eout
              if (n.ge.2) then
                branchplus=g2E/(gsg*(n-2)+g2E)
              else
                branchplus=0.
              endif
              branchzero=gsg*n/(gsg*n+g2E)
              phratio=(branchplus*phres1+branchzero*phres2)/phcompg
              wemission(type,p,h,nen)=Rgamma*factor*phratio
              wemissum(p,h,nen)=wemission(type,p,h,nen)
            endif
          else
c
c 2. Particle emission rates
c
c pres : help variable
c nres : exciton number for residual system
c phres: particle-hole state density for residual system
c
            pres=p-aejec
            nres=n-aejec
            if (pres.lt.0.or.h.eq.0) goto 40
            Ures=max(Eres-preeqpair(Zix,Nix,nres,Eres,pairmodel),
     +        preeqpair(Zix,Nix,nres,Eres,pairmodel))
            phres=phdens(Zix,Nix,pres,h,gs,Ures,Ewell,surfwell)
            phratio=phres/phcomp
            wemission(type,p,h,nen)=factor*phratio*Qfac
            wemissum(p,h,nen)=wemissum(p,h,nen)+wemission(type,p,h,nen)
          endif
c
c *** Integration of emission rates over all energies and particles ****
c
c deltaE: energy bin around outgoing energies
c
          wemispart(type,p,h)=wemispart(type,p,h)+
     +      wemission(type,p,h,nen)*deltaE(nen)
   40   continue
        goto 60
c
c Correction in integration for last outgoing energy.
c
c Emax: maximal outgoing energy
c dE  : extra part for energy integration
c
   50   Eout=egrid(nen1)
        Emax=Ecomp-S(Zcomp,Ncomp,type)
        dE=Emax-(Eout+0.5*deltaE(nen1))
        wemispart(type,p,h)=wemispart(type,p,h)+
     +    wemission(type,p,h,nen1)*dE
   60   wemistot(p,h)=wemistot(p,h)+wemispart(type,p,h)
   20 continue
c
c Prevent divergence of pre-equilibrium gamma cross sections in case
c of absence of particle competition.
c
c
      if (n.gt.1) then
        if (wemistot(p,h).eq.wemispart(0,p,h)) wemispart(0,p,h)=0.
        do 110 nen=0,numen
          if (wemissum(p,h,nen).eq.wemission(0,p,h,nen)) then
            wemissum(p,h,nen)=0.
            wemission(0,p,h,nen)=0.
          endif
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
