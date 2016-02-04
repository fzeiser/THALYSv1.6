      function lambdapiplus(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 1, 2008
c | Task  : Proton transition rates for n --> n+2
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          surfwell
      integer          Zcomp,Ncomp,ppi,hpi,pnu,hnu,h,p,n,A,nexcbins,i,j,
     +                 Zix,Nix,k,nen
      real             lambdapiplus,edepth,gsp,gsn,damp,ignatyuk,U,
     +                 preeqpair,fac1,factor1,factor2,factor3,term1,
     +                 term2,term12,L1pip,L2pip,L1pih,L2pih,L1nup,L2nup,
     +                 L1nuh,L2nuh,dExpip,dExpih,dExnup,dExnuh,uupip,
     +                 uupih,uunup,uunuh,lambdapipi1p,lambdapipi1h,
     +                 lambdanupi1p,lambdanupi1h,uu,eopt,Weff,phdens2,
     +                 densh,densp,ratio,finitewell
      double precision lplus,sumpipi1p,sumpipi1h,sumnupi1p,sumnupi1h,
     +                 termpipi1p,termpipi1h,termnupi1p,termnupi1h,phtot
c
c *************************** Transition rates *************************
c
c lambdapiplus: proton transition rate for n --> n+2
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c ppi         : proton particle number
c hpi         : proton hole number
c pnu         : neutron particle number
c hnu         : neutron hole number
c h           : hole number
c p           : particle number
c n           : exciton number
c A           : mass number of compound nucleus
c Ainit       : mass number of initial compound nucleus
c matrix      : subroutine for matrix element for exciton model
c surfwell    : flag for surface effects in finite well
c flagsurface : flag for surface effects in exciton model
c primary     : flag to designate primary (binary) reaction
c edepth      : depth of potential well
c Efermi      : depth of Fermi well
c Esurf       : well depth for surface interaction
c gp,gsp      : single-particle proton level density parameter
c gn,gsn      : single-particle neutron level density parameter
c flaggshell  : flag for energy dependence of single particle level
c               density parameter g
c damp        : shell damping factor
c ignatyuk    : function for energy dependent level density parameter a
c Ecomp       : total energy of composite system
c alev        : level density parameter
c U           : excitation energy minus pairing energy
c preeqpair   : pre-equilibrium pairing energy
c pairmodel   : model for preequilibrium pairing energy
c fac1,lplus  : help variable
c factor1-3   : help variables
c twopihbar   : 2*pi/hbar
c Apauli2     : two-component Pauli blocking correction factor
c term1-2     : help variables
c M2pipi      : square of proton-proton matrix element
c M2pinu      : square of proton-neutron matrix element
c finitewell  : correction function for finite well depth
c
c A. Analytical solution: The transition rates are taken from Kalbach,
c    PRC33, 818 (1986).
c
      lambdapiplus=0.
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
        fac1=2.*n*(n+1.)
        factor1=twopihbar*gsp*gsp/fac1
        term1=U-Apauli2(ppi+1,hpi+1,pnu,hnu)
        term2=U-Apauli2(ppi,hpi,pnu,hnu)
        if (term1.le.0.or.term2.le.0.) return
        term12=term1/term2
        if (term12.lt.0.01) return
        factor2=term1**2*(term12**(n-1))
        factor3=(ppi+hpi)*gsp*M2pipi+2.*(pnu+hnu)*gsn*M2pinu
        lplus=factor1*factor2*factor3
        lambdapiplus=lplus*finitewell(p+1,h+1,U,edepth,surfwell)
      else
c
c B. Numerical solution: Transition rates based on either matrix
c    element (preeqmode=2) or optical model (preeqmode=3).
c
c L1pip,...   : integration limits
c nexcbins    : number of integration bins
c nbins       : number of continuum excitation energy bins
c dExpip,...  : integration bin width
c sumpipi1p...: help variables
c uupip,..... : residual excitation energy
c lambdapipi1p: collision probability for proton-proton particle
c lambdapipi1h: collision probability for proton-proton hole
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
c termpipi1p,.: help variables
c phtot       : total particle-hole state density
c
        L1pip=Apauli2(ppi+1,hpi+1,pnu,hnu)-Apauli2(ppi-1,hpi,pnu,hnu)
        L2pip=U-Apauli2(ppi-1,hpi,pnu,hnu)
        L1pih=Apauli2(ppi+1,hpi+1,pnu,hnu)-Apauli2(ppi,hpi-1,pnu,hnu)
        L2pih=U-Apauli2(ppi,hpi-1,pnu,hnu)
        L1nup=Apauli2(ppi+1,hpi+1,pnu,hnu)-Apauli2(ppi,hpi,pnu-1,hnu)
        L2nup=U-Apauli2(ppi,hpi,pnu-1,hnu)
        L1nuh=Apauli2(ppi+1,hpi+1,pnu,hnu)-Apauli2(ppi,hpi,pnu,hnu-1)
        L2nuh=U-Apauli2(ppi,hpi,pnu,hnu-1)
        if (primary) then
          nexcbins=max(nbins/2,2)
        else
          nexcbins=max(nbins/4,2)
        endif
        dExpip=(L2pip-L1pip)/nexcbins
        dExpih=(L2pih-L1pih)/nexcbins
        dExnup=(L2nup-L1nup)/nexcbins
        dExnuh=(L2nuh-L1nuh)/nexcbins
        sumpipi1p=0.
        sumpipi1h=0.
        sumnupi1p=0.
        sumnupi1h=0.
        do 10 i=1,nexcbins
          uupip=L1pip+(i-0.5)*dExpip
          uupih=L1pih+(i-0.5)*dExpih
          uunup=L1nup+(i-0.5)*dExnup
          uunuh=L1nuh+(i-0.5)*dExnuh
          if (flaggshell) then
            damp=ignatyuk(Zcomp,Ncomp,uupip,0)/alev(Zcomp,Ncomp)
            gsp=gp(Zcomp,Ncomp)*damp
            gsn=gn(Zcomp,Ncomp)*damp
          endif
          if (preeqmode.eq.2) then
            lambdapipi1p=twopihbar*M2pipi*
     +        phdens2(Zcomp,Ncomp,2,1,0,0,gsp,gsn,uupip,edepth,surfwell)
            lambdapipi1h=twopihbar*M2pipi*
     +        phdens2(Zcomp,Ncomp,1,2,0,0,gsp,gsn,uupih,edepth,surfwell)
            lambdanupi1p=twopihbar*M2nupi*
     +        phdens2(Zcomp,Ncomp,1,1,1,0,gsp,gsn,uunup,edepth,surfwell)
            lambdanupi1h=twopihbar*M2nupi*
     +        phdens2(Zcomp,Ncomp,1,1,0,1,gsp,gsn,uunuh,edepth,surfwell)
          else
            do 20 j=1,4
              if (j.eq.1) uu=uupip
              if (j.eq.2) uu=uupih
              if (j.eq.3) uu=uunup
              if (j.eq.4) uu=uunuh
              if (j.gt.2) then
                Zix=0
                Nix=1
                k=1
              else
                Zix=1
                Nix=0
                k=2
              endif
              eopt=max(uu-S(Zix,Nix,k),-20.)
              nen=min(10*numen,int(eopt*10.))
              if (j.le.2) then
                Weff=Wompfac(1)*wvol(k,nen)
              else
                Weff=Wompfac(2)*wvol(k,nen)
              endif
              if (j.eq.1) lambdapipi1p=2.*Weff/hbar
              if (j.eq.2) then
                densh=phdens2(Zcomp,Ncomp,1,2,0,0,gsp,gsn,uu,
     +            edepth,surfwell)
                densp=phdens2(Zcomp,Ncomp,2,1,0,0,gsp,gsn,uu,
     +            edepth,surfwell)
                if (densp.gt.1.) then
                  ratio=densh/densp
                else
                  ratio=1.
                endif
                lambdapipi1h=2.*Weff/hbar*ratio
              endif
              if (j.eq.3) lambdanupi1p=2.*Weff/hbar
              if (j.eq.4) then
                densh=phdens2(Zcomp,Ncomp,1,1,0,1,gsp,gsn,uu,
     +            edepth,surfwell)
                densp=phdens2(Zcomp,Ncomp,1,1,1,0,gsp,gsn,uu,
     +            edepth,surfwell)
                if (densp.gt.1.) then
                  ratio=densh/densp
                else
                  ratio=1.
                endif
                lambdanupi1h=2.*Weff/hbar*ratio
              endif
   20       continue
          endif
          termpipi1p=lambdapipi1p*gsp*dExpip
          termpipi1h=lambdapipi1h*gsp*dExpih
          termnupi1p=lambdanupi1p*gsn*dExnup
          termnupi1h=lambdanupi1h*gsn*dExnuh
          if (flaggshell) then
            damp=ignatyuk(Zcomp,Ncomp,U-uupip,0)/alev(Zcomp,Ncomp)
            gsp=gp(Zcomp,Ncomp)*damp
            gsn=gn(Zcomp,Ncomp)*damp
          endif
          sumpipi1p=sumpipi1p+termpipi1p*
     +      phdens2(Zcomp,Ncomp,ppi-1,hpi,pnu,hnu,gsp,gsn,U-uupip,
     +      edepth,surfwell)
          sumpipi1h=sumpipi1h+termpipi1h*
     +      phdens2(Zcomp,Ncomp,ppi,hpi-1,pnu,hnu,gsp,gsn,U-uupih,
     +      edepth,surfwell)
          sumnupi1p=sumnupi1p+termnupi1p*
     +      phdens2(Zcomp,Ncomp,ppi,hpi,pnu-1,hnu,gsp,gsn,U-uunup,
     +      edepth,surfwell)
          sumnupi1h=sumnupi1h+termnupi1h*
     +      phdens2(Zcomp,Ncomp,ppi,hpi,pnu,hnu-1,gsp,gsn,U-uunuh,
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
          lplus=(sumpipi1p+sumpipi1h+sumnupi1p+sumnupi1h)/phtot
          lambdapiplus=real(lplus)
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
