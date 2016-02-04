      function lambdaplus(Zcomp,Ncomp,p,h)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 1, 2008
c | Task  : Transition rates for n --> n+2
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          surfwell
      integer          Zcomp,Ncomp,p,h,n,A,nexcbins,Zix,Nix,i,j,k,nen
      real             lambdaplus,edepth,gs,ignatyuk,U,preeqpair,
     +                 factor1,term1,term2,term12,finitewell,L1p,L2p,
     +                 L1h,L2h,dExpp,dExhh,uup,uuh,lambda1p,lambda1h,
     +                 Zratio,Nratio,uu,eopt,Weff(2),ww,densh,densp,
     +                 ratio,phdens
      double precision lplus,sum1p,sum1h,term1p,term1h,phtot
c
c *************************** Transition rates *************************
c
c lambdaplus : transition rate for n --> n+2
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c p          : particle number
c h          : hole number
c n          : exciton number
c A          : mass number of compound nucleus
c Ainit      : mass number of initial compound nucleus
c matrix     : subroutine for matrix element for exciton model
c surfwell   : flag for surface effects in finite well
c flagsurface: flag for surface effects in exciton model
c primary    : flag to designate primary (binary) reaction
c edepth     : depth of potential well
c Esurf      : well depth for surface interaction
c Efermi     : depth of Fermi well
c gs,g       : single-particle level density parameter
c flaggshell : flag for energy dependence of single particle level
c              density parameter g
c ignatyuk   : function for energy dependent level density parameter a
c Ecomp      : total energy of composite system
c alev       : level density parameter
c U          : excitation energy minus pairing energy
c preeqpair  : pre-equilibrium pairing energy
c pairmodel  : model for preequilibrium pairing energy
c
      lambdaplus=0.
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
      gs=g(Zcomp,Ncomp)
      if (flaggshell) gs=gs*ignatyuk(Zcomp,Ncomp,Ecomp,0)/
     +  alev(Zcomp,Ncomp)
      U=Ecomp-preeqpair(Zcomp,Ncomp,n,Ecomp,pairmodel)
c
c A. Analytical solution
c    First we calculate the transition density for the infinite well as
c    given by Oblozinsky et al (Nuc Phys A226, 347 (1974)) and then
c    multiply it by the finite well correction function.
c
c preeqmode    : designator for pre-equilibrium model
c factor1,lplus: help variables
c twopihbar    : 2*pi/hbar
c M2           : square of matrix element
c term1,term2  : help variables
c Apauli       : Pauli blocking correction energy
c finitewell   : correction function for finite well depth
c
      if ((preeqmode.ne.2.and.preeqmode.ne.3).or.n.eq.1) then
        factor1=twopihbar*M2*0.5/(n+1.)*gs**3
        term1=U-Apauli(p+1,h+1)
        term2=U-Apauli(p,h)
        if (term1.le.0.or.term2.le.0.) return
        term12=term1/term2
        if (term12.lt.0.01) return
        lplus=factor1*term1**2*(term12**(n-1))
        lambdaplus=lplus*finitewell(p+1,h+1,U,edepth,surfwell)
      else
c
c B. Numerical solution: Transition rates based on either matrix
c    element (preeqmode=2) or optical model (preeqmode=3).
c
c L1p,...      : integration limits
c nexcbins     : number of integration bins
c nbins        : number of continuum excitation energy bins
c dExpp....    : integration bin width
c sum1p,sum1h  : help variables
c uup,uuh      : residual excitation energy
c lambda1p     : collision probability for particle
c lambda1h     : collision probability for hole
c phdens       : particle-hole state density
c Zratio,Nratio: help variables
c Zinit        : charge number of initial compound nucleus
c Ninit        : neutron number of initial compound nucleus
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c eopt         : energy for optical model calculation
c S            : separation energy per particle
c Weff         : effective imaginary well depth
c Wompfac      : adjustable constant for OMP based transition rates
c wvol         : absorption part of the optical potential averaged over
c                the volume
c densh,densp  : help variables
c term1p,term1h: help variables
c ratio        : state density ratio for hole scattering
c hbar         : Planck's constant / 2.pi in MeV.s
c phtot        : total particle-hole state density
c
        L1p=Apauli(p+1,h+1)-Apauli(p-1,h)
        L2p=U-Apauli(p-1,h)
        L1h=Apauli(p+1,h+1)-Apauli(p,h-1)
        L2h=U-Apauli(p,h-1)
        if (primary) then
          nexcbins=max(nbins/2,2)
        else
          nexcbins=max(nbins/4,2)
        endif
        dExpp=(L2p-L1p)/nexcbins
        dExhh=(L2h-L1h)/nexcbins
        sum1p=0.
        sum1h=0.
        do 10 i=1,nexcbins
          uup=L1p+(i-0.5)*dExpp
          uuh=L1h+(i-0.5)*dExhh
          if (flaggshell) gs=g(Zcomp,Ncomp)*
     +      ignatyuk(Zcomp,Ncomp,uup,0)/alev(Zcomp,Ncomp)
          if (preeqmode.eq.2) then
            lambda1p=twopihbar*M2*phdens(Zcomp,Ncomp,2,1,gs,uup,
     +        edepth,surfwell)
            lambda1h=twopihbar*M2*phdens(Zcomp,Ncomp,1,2,gs,uuh,
     +        edepth,surfwell)
          else
            Zratio=real(Zinit)/Ainit
            Nratio=real(Ninit)/Ainit
            do 20 j=1,2
              if (j.eq.1) then
                uu=uup
              else
                uu=uuh
              endif
              do 30 k=1,2
                if (k.eq.1) then
                  Zix=0
                  Nix=1
                else
                  Zix=1
                  Nix=0
                endif
                eopt=max(uu-S(Zix,Nix,k),-20.)
                nen=min(10*numen,int(eopt*10.))
                Weff(k)=0.5*Wompfac(0)*wvol(k,nen)
   30         continue
              ww=Zratio*Weff(2)+Nratio*Weff(1)
              if (j.eq.1) then
                lambda1p=2.*ww/hbar
              else
                densh=phdens(Zcomp,Ncomp,1,2,gs,uu,edepth,surfwell)
                densp=phdens(Zcomp,Ncomp,2,1,gs,uu,edepth,surfwell)
                if (densp.gt.1.) then
                  ratio=densh/densp
                else
                  ratio=1.
                endif
                lambda1h=2.*ww/hbar*ratio
              endif
   20       continue
          endif
          term1p=lambda1p*gs*dExpp
          term1h=lambda1h*gs*dExhh
          if (flaggshell) gs=g(Zcomp,Ncomp)*
     +      ignatyuk(Zcomp,Ncomp,U-uup,0)/alev(Zcomp,Ncomp)
          sum1p=sum1p+term1p*phdens(Zcomp,Ncomp,p-1,h,gs,U-uup,
     +      edepth,surfwell)
          sum1h=sum1h+term1h*phdens(Zcomp,Ncomp,p,h-1,gs,U-uuh,
     +      edepth,surfwell)
   10   continue
        if (flaggshell) gs=g(Zcomp,Ncomp)*
     +    ignatyuk(Zcomp,Ncomp,U,0)/alev(Zcomp,Ncomp)
        phtot=phdens(Zcomp,Ncomp,p,h,gs,U,edepth,surfwell)
        if (phtot.gt.0.) then
          lplus=(sum1p+sum1h)/phtot
          lambdaplus=real(lplus)
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
