      subroutine emissionrate2(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 14, 2010
c | Task  : Two-component emission rates for exciton model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell,surfgam
      integer Zcomp,Ncomp,ppi,hpi,pnu,hnu,n,h,type,Zix,Nix,zejec,nejec,
     +        nen,nen1,ppires,pnures,nres
      real    gsp,gsn,gs,damp,ignatyuk,edepth,U,preeqpair,phcomp,
     +        wemissum2(0:numparx,0:numparx,0:numparx,0:numparx,
     +        0:numen),phdens2,Ewell,Eout,xs,factor,Eres,phres1,phres2,
     +        g2E,branchplus,branchzero,phratio,Ures,phres,Emax,dE
c
c *************************** Emission rates ***************************
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c ppi        : proton particle number
c hpi        : proton hole number
c pnu        : neutron particle number
c hnu        : neutron hole number
c n          : exciton number
c h          : hole number
c gp,gsp     : single-particle proton level density parameter
c gn,gsn     : single-particle neutron level density parameter
c gs         : single-particle level density parameter
c flaggshell : flag for energy dependence of single particle level
c              density parameter g
c damp       : shell damping factor
c ignatyuk   : function for energy dependent level density parameter a
c Ecomp      : total energy of composite system
c alev       : level density parameter
c surfwell   : flag for surface effects in finite well
c flagsurface: flag for surface effects in exciton model
c primary    : flag to designate primary (binary) reaction
c edepth     : depth of potential well
c Esurf      : well depth for surface interaction
c Efermi     : depth of Fermi well
c U,Ures     : excitation energy minus pairing energy
c preeqpair  : pre-equilibrium pairing energy
c pairmodel  : model for preequilibrium pairing energy
c phcomp     : particle-hole state density for compound system
c phdens2    : function for two-component particle-hole state density
c wemistot2  : total two-component emission rate per exciton number
c wemispart2 : two-component emission rate per particle and exciton
c              number
c wemission2 : two-component emission rate per particle, exciton
c              number and energy
c wemissum2  : two-component emission rate per exciton number and energy
c parskip    : logical to skip outgoing particle
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c zejec,parZ : charge number of leading particle
c nejec,parN : neutron number of leading particle
c Ewell      : depth of potential well
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c egrid,Eout : energies of basic energy grid in MeV
c xs,factor  : help variables
c xsreac     : reaction cross section
c wfac       : factor for emission rate
c Eres       : total energy of residual system
c S          : separation energy per particle
c
      n=ppi+hpi+pnu+hnu
      h=hpi+hnu
      gsp=gp(Zcomp,Ncomp)
      gsn=gn(Zcomp,Ncomp)
      gs=gsp+gsn
      if (flaggshell) then
        damp=ignatyuk(Zcomp,Ncomp,Ecomp,0)/alev(Zcomp,Ncomp)
        gsp=gsp*damp
        gsn=gsn*damp
        gs=gsp+gsn
      endif
      surfwell=flagsurface.and.h.eq.1.and.primary
      if (surfwell) then
        edepth=Esurf
      else
        edepth=Efermi
      endif
      U=Ecomp-preeqpair(Zcomp,Ncomp,n,Ecomp,pairmodel)
      phcomp=phdens2(Zcomp,Ncomp,ppi,hpi,pnu,hnu,gsp,gsn,U,edepth,
     +  surfwell)
      wemistot2(ppi,hpi,pnu,hnu)=0.
      do 10 nen=0,numen
        wemissum2(ppi,hpi,pnu,hnu,nen)=0.
   10 continue
      do 20 type=0,6
        wemispart2(type,ppi,hpi,pnu,hnu)=0.
        do 30 nen=0,numen
          wemission2(type,ppi,hpi,pnu,hnu,nen)=0.
   30   continue
        if (parskip(type)) goto 20
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        gsp=gp(Zix,Nix)
        gsn=gn(Zix,Nix)
        zejec=parZ(type)
        nejec=parN(type)
        if (type.gt.2) then
          Ewell=Efermi
        else
          Ewell=edepth
        endif
        if (phcomp.le.1.e-10) goto 20
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
          if (flaggshell) then
            damp=ignatyuk(Zix,Nix,Eres,0)/alev(Zix,Nix)
            gsp=gp(Zix,Nix)*damp
            gsn=gn(Zix,Nix)*damp
            gs=gsp+gsn
          endif
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
              phres1=0.5*(phdens2(Zcomp,Ncomp,ppi-1,hpi-1,pnu,hnu,
     +          gsp,gsn,U,Efermi,surfgam)+
     +          phdens2(Zcomp,Ncomp,ppi,hpi,pnu-1,hnu-1,
     +            gsp,gsn,U,Efermi,surfgam))
              phres2=phdens2(Zcomp,Ncomp,ppi,hpi,pnu,hnu,
     +          gsp,gsn,U,Efermi,surfgam)
              g2E=gs*gs*Eout
              if (n.ge.2) then
                branchplus=g2E/(gs*(n-2)+g2E)
              else
                branchplus=0.
              endif
              branchzero=gs*n/(gs*n+g2E)
              phratio=(branchplus*phres1+branchzero*phres2)/phcomp
              wemission2(type,ppi,hpi,pnu,hnu,nen)=Rgamma*factor*phratio
              wemissum2(ppi,hpi,pnu,hnu,nen)=
     +          wemission2(type,ppi,hpi,pnu,hnu,nen)
            endif
          else
c
c 2. Particle emission rates
c
c ppires,pnures: help variables
c nres         : exciton number for residual system
c phres        : particle-hole state density for residual system
c
            ppires=ppi-zejec
            pnures=pnu-nejec
            nres=n-zejec-nejec
            if (ppires.lt.0.or.pnures.lt.0.or.h.eq.0) goto 40
            Ures=max(Eres-preeqpair(Zix,Nix,nres,Eres,pairmodel),
     +        preeqpair(Zix,Nix,nres,Eres,pairmodel))
            phres=phdens2(Zix,Nix,ppires,hpi,pnures,hnu,
     +        gsp,gsn,Ures,Ewell,surfwell)
            phratio=phres/phcomp
            wemission2(type,ppi,hpi,pnu,hnu,nen)=factor*phratio
            wemissum2(ppi,hpi,pnu,hnu,nen)=
     +        wemissum2(ppi,hpi,pnu,hnu,nen)+
     +        wemission2(type,ppi,hpi,pnu,hnu,nen)
          endif
c
c *** Integration of emission rates over all energies and particles ****
c
c deltaE: energy bin around outgoing energies
c
          wemispart2(type,ppi,hpi,pnu,hnu)=
     +      wemispart2(type,ppi,hpi,pnu,hnu)+
     +      wemission2(type,ppi,hpi,pnu,hnu,nen)*deltaE(nen)
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
        wemispart2(type,ppi,hpi,pnu,hnu)=
     +    wemispart2(type,ppi,hpi,pnu,hnu)+
     +    wemission2(type,ppi,hpi,pnu,hnu,nen1)*dE
   60   wemistot2(ppi,hpi,pnu,hnu)=wemistot2(ppi,hpi,pnu,hnu)+
     +    wemispart2(type,ppi,hpi,pnu,hnu)
   20 continue
c
c Prevent divergence of pre-equilibrium gamma cross sections in case
c of absence of particle competition.
c
c
      if (n.gt.1) then
        if (wemistot2(ppi,hpi,pnu,hnu).eq.wemispart2(0,ppi,hpi,pnu,hnu))
     +    then
          wemistot2(ppi,hpi,pnu,hnu)=0.
          wemispart2(0,ppi,hpi,pnu,hnu)=0.
        endif
        do 110 nen=0,numen
          if (wemissum2(ppi,hpi,pnu,hnu,nen).eq.
     +      wemission2(0,ppi,hpi,pnu,hnu,nen))
     +      wemission2(0,ppi,hpi,pnu,hnu,nen)=0.
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
