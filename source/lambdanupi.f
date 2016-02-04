      function lambdanupi(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 1, 2008
c | Task  : Neutron-proton transition rates for n --> n
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          surfwell
      integer          Zcomp,Ncomp,ppi,hpi,pnu,hnu,h,p,n,A,nexcbins,i,
     +                 Zix,Nix,k,nen
      real             lambdanupi,edepth,gsp,gsn,damp,ignatyuk,U,
     +                 preeqpair,factor1,factor2,factor3,factor4,
     +                 factor23,Bfactor,L1,L2,dEx,uu,lambdanupi1p,eopt,
     +                 Weff,phdens2,densh,densp,ratio,finitewell
      double precision lnupi,sumnupi1p,termnupi1p,phtot
c
c *************************** Transition rates *************************
c
c lambdanupi   : proton-neutron transition rate for n --> n
c Zcomp        : charge number index for compound nucleus
c Ncomp        : neutron number index for compound nucleus
c ppi          : proton particle number
c hpi          : proton hole number
c pnu          : neutron particle number
c hnu          : neutron hole number
c h            : hole number
c p            : particle number
c n            : exciton number
c A            : mass number of compound nucleus
c Ainit        : mass number of initial compound nucleus
c matrix       : subroutine for matrix element for exciton model
c surfwell     : flag for surface effects in finite well
c flagsurface  : flag for surface effects in exciton model
c primary      : flag to designate primary (binary) reaction
c edepth       : depth of potential well
c Efermi       : depth of Fermi well
c Esurf        : well depth for surface interaction
c gp,gsp       : single-particle proton level density parameter
c gn,gsn       : single-particle neutron level density parameter
c flaggshell   : flag for energy dependence of single particle level
c                density parameter g
c damp         : shell damping factor
c ignatyuk     : function for energy dependent level density parameter a
c Ecomp        : total energy of composite system
c alev         : level density parameter
c U            : excitation energy minus pairing energy
c preeqpair    : pre-equilibrium pairing energy
c pairmodel    : model for preequilibrium pairing energy
c Bfactor,lnupi: help variables
c twopihbar    : 2*pi/hbar
c factor1-4    : help variables
c Apauli2      : two-component Pauli blocking correction factor
c M2nupi       : square of neutron-proton matrix element
c finitewell   : correction function for finite well depth
c
c A. Analytical solution: The transition rates are taken from Kalbach,
c    PRC33, 818 (1986).
c
      lambdanupi=0.
      h=hpi+hnu
      p=ppi+pnu
      n=p+h
      if (n.eq.0) return
      A=Ainit-Zcomp-Ncomp
      call matrix(A,n)
      surfwell=flagsurface.and.h.eq.1.and.primary
      if (surfwell) then
        edepth=Esurf
      else
        edepth=Efermi
      endif
      gsp=gp(Zcomp,Ncomp)
      gsn=gn(Zcomp,Ncomp)
      if (flaggshell) then
        damp=ignatyuk(Zcomp,Ncomp,Ecomp,0)/alev(Zcomp,Ncomp)
        gsp=gsp*damp
        gsn=gsn*damp
      endif
      U=Ecomp-preeqpair(Zcomp,Ncomp,n,Ecomp,pairmodel)
      if ((preeqmode.ne.2.and.preeqmode.ne.3).or.n.eq.1) then
        factor1=twopihbar*pnu*hnu*M2nupi/n*gsp*gsp
        Bfactor=Apauli2(ppi,hpi,pnu,hnu)
        if (pnu.gt.0.and.hnu.gt.0)
     +    Bfactor=max(Bfactor,Apauli2(ppi+1,hpi+1,pnu-1,hnu-1))
        factor2=U-Bfactor
        factor3=U-Apauli2(ppi,hpi,pnu,hnu)
        if (factor2.le.0..or.factor3.le.0.) return
        factor23=factor2/factor3
        if (factor23.lt.0.01) return
        factor4=2.*(U-Bfactor)
        if (pnu.gt.0.and.hnu.gt.0) factor4=factor4+n*
     +    abs(Apauli2(ppi,hpi,pnu,hnu)-Apauli2(ppi+1,hpi+1,pnu-1,hnu-1))
        lnupi=factor1*(factor23**(n-1))*factor4
        lambdanupi=lnupi*finitewell(p,h,U,edepth,surfwell)
      else
c
c B. Numerical solution: Transition rates based on either matrix
c    element (preeqmode=2) or optical model (preeqmode=3).
c
c L1,L2       : integration limits
c nexcbins    : number of integration bins
c nbins       : number of continuum excitation energy bins
c dEx         : integration bin width
c sumnupi1p   : help variable
c uu          : residual excitation energy
c lambdanupi1p: collision probability for neutron-proton particle
c phdens2     : two-component particle-hole state density
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c eopt        : energy for optical model calculation
c S           : separation energy per particle
c Weff        : effective imaginary well depth
c Wompfac     : adjustable constant for OMP based transition rates
c wvol        : absorption part of the optical potential averaged over
c               the volume
c hbar        : Planck's constant / 2.pi in MeV.s
c densh,densp : help variables
c ratio       : state density ratio for hole scattering
c termnupi1p,.: help variables
c phtot       : total particle-hole state density
c
        L1=Apauli2(ppi,hpi,pnu,hnu)-Apauli2(ppi,hpi,pnu-1,hnu-1)
        L2=U-Apauli2(ppi,hpi,pnu-1,hnu-1)
        if (primary) then
          nexcbins=max(nbins/2,2)
        else
          nexcbins=max(nbins/4,2)
        endif
        dEx=(L2-L1)/nexcbins
        sumnupi1p=0.
        do 10 i=1,nexcbins
          uu=L1+(i-0.5)*dEx
          if (flaggshell) then
            damp=ignatyuk(Zcomp,Ncomp,uu,0)/alev(Zcomp,Ncomp)
            gsp=gp(Zcomp,Ncomp)*damp
            gsn=gn(Zcomp,Ncomp)*damp
          endif
          if (preeqmode.eq.2) then
            lambdanupi1p=twopihbar*M2nupi*
     +        phdens2(Zcomp,Ncomp,1,1,0,0,gsp,gsn,uu,edepth,surfwell)
          else
            Zix=0
            Nix=1
            k=1
            eopt=max(uu-S(Zix,Nix,k),-20.)
            nen=min(10*numen,int(eopt*10.))
            Weff=Wompfac(0)*wvol(k,nen)
            densh=phdens2(Zcomp,Ncomp,1,1,0,0,gsp,gsn,uu,
     +        edepth,surfwell)
            densp=phdens2(Zcomp,Ncomp,1,1,1,0,gsp,gsn,uu,
     +        edepth,surfwell)
            if (densp.gt.1.) then
              ratio=densh/densp
            else
              ratio=1.
            endif
            lambdanupi1p=2.*Weff/hbar*ratio
          endif
          termnupi1p=lambdanupi1p*
     +      phdens2(Zcomp,Ncomp,0,0,1,1,gsp,gsn,uu,edepth,surfwell)*dEx
          if (flaggshell) then
            damp=ignatyuk(Zcomp,Ncomp,uu,0)/alev(Zcomp,Ncomp)
            gsp=gp(Zcomp,Ncomp)*damp
            gsn=gn(Zcomp,Ncomp)*damp
          endif
          sumnupi1p=sumnupi1p+termnupi1p*
     +      phdens2(Zcomp,Ncomp,ppi,hpi,pnu-1,hnu-1,gsp,gsn,U-uu,
     +      edepth,surfwell)
   10   continue
        if (flaggshell) then
          damp=ignatyuk(Zcomp,Ncomp,U,0)/alev(Zcomp,Ncomp)
          gsp=gp(Zcomp,Ncomp)*damp
          gsn=gn(Zcomp,Ncomp)*damp
        endif
        phtot=phdens2(Zcomp,Ncomp,ppi,hpi,pnu,hnu,gsp,gsn,U,
     +    edepth,surfwell)
        if (phtot.gt.0.) then
          lnupi=sumnupi1p/phtot
          lambdanupi=real(lnupi)
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
