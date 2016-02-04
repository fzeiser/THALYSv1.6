      function lambdanuplus(Zcomp,Ncomp,ppi,hpi,pnu,hnu)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 1, 2008
c | Task  : Neutron transition rates for n --> n+2
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          surfwell
      integer          Zcomp,Ncomp,ppi,hpi,pnu,hnu,h,p,n,A,nexcbins,i,
     +                 j,Zix,Nix,k,nen
      real             lambdanuplus,edepth,gsp,gsn,damp,ignatyuk,U,
     +                 preeqpair,fac1,factor1,factor2,factor3,term1,
     +                 term2,term12,L1nup,L2nup,L1nuh,L2nuh,L1pip,L2pip,
     +                 L1pih,L2pih,dExnup,dExnuh,dExpip,dExpih,uunup,
     +                 uunuh,uupip,uupih,lambdanunu1p,lambdanunu1h,
     +                 lambdapinu1p,lambdapinu1h,uu,eopt,Weff,phdens2,
     +                 densh,densp,ratio,finitewell
      double precision lplus,sumnunu1p,sumnunu1h,sumpinu1p,sumpinu1h,
     +                 termnunu1p,termnunu1h,termpinu1p,termpinu1h,phtot
c
c *************************** Transition rates *************************
c
c lambdanuplus: neutron transition rate for n --> n+2
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
c twopihbar   : 2*pi/hbar
c factor1-3   : help variables
c Apauli2     : two-component Pauli blocking correction factor
c term1-2     : help variables
c M2nunu      : square of neutron-neutron matrix element
c M2nupi      : square of neutron-proton matrix element
c finitewell  : correction function for finite well depth
c
c A. Analytical solution: The transition rates are taken from Kalbach,
c    PRC33, 818 (1986).
c
      lambdanuplus=0.
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
        factor1=twopihbar*gsn*gsn/fac1
        term1=U-Apauli2(ppi,hpi,pnu+1,hnu+1)
        term2=U-Apauli2(ppi,hpi,pnu,hnu)
        if (term1.le.0.or.term2.le.0.) return
        term12=term1/term2
        if (term12.lt.0.01) return
        factor2=term1**2*(term12**(n-1))
        factor3=(pnu+hnu)*gsn*M2nunu+2.*(ppi+hpi)*gsp*M2nupi
        lplus=factor1*factor2*factor3
        lambdanuplus=lplus*finitewell(p+1,h+1,U,edepth,surfwell)
      else
c
c B. Numerical solution: Transition rates based on either matrix
c    element (preeqmode=2) or optical model (preeqmode=3).
c
c L1nup,...   : integration limits
c nexcbins    : number of integration bins
c nbins       : number of continuum excitation energy bins
c dExnup,...  : integration bin width
c sumnunu1p...: help variables
c uunup,..... : residual excitation energy
c lambdanunu1p: collision probability for neutron-neutron particle
c lambdanunu1h: collision probability for neutron-neutron hole
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
c termpipip,..: help variables
c phtot       : total particle-hole state density
c
        L1nup=Apauli2(ppi,hpi,pnu+1,hnu+1)-Apauli2(ppi,hpi,pnu-1,hnu)
        L2nup=U-Apauli2(ppi,hpi,pnu-1,hnu)
        L1nuh=Apauli2(ppi,hpi,pnu+1,hnu+1)-Apauli2(ppi,hpi,pnu,hnu-1)
        L2nuh=U-Apauli2(ppi,hpi,pnu,hnu-1)
        L1pip=Apauli2(ppi,hpi,pnu+1,hnu+1)-Apauli2(ppi-1,hpi,pnu,hnu)
        L2pip=U-Apauli2(ppi-1,hpi,pnu,hnu)
        L1pih=Apauli2(ppi,hpi,pnu+1,hnu+1)-Apauli2(ppi,hpi-1,pnu,hnu)
        L2pih=U-Apauli2(ppi,hpi-1,pnu,hnu)
        if (primary) then
          nexcbins=max(nbins/2,2)
        else
          nexcbins=max(nbins/4,2)
        endif
        dExnup=(L2nup-L1nup)/nexcbins
        dExnuh=(L2nuh-L1nuh)/nexcbins
        dExpip=(L2pip-L1pip)/nexcbins
        dExpih=(L2pih-L1pih)/nexcbins
        sumnunu1p=0.
        sumnunu1h=0.
        sumpinu1p=0.
        sumpinu1h=0.
        do 10 i=1,nexcbins
          uunup=L1nup+(i-0.5)*dExnup
          uunuh=L1nuh+(i-0.5)*dExnuh
          uupip=L1pip+(i-0.5)*dExpip
          uupih=L1pih+(i-0.5)*dExpih
          if (flaggshell) then
            damp=ignatyuk(Zcomp,Ncomp,uunup,0)/alev(Zcomp,Ncomp)
            gsp=gp(Zcomp,Ncomp)*damp
            gsn=gn(Zcomp,Ncomp)*damp
          endif
          if (preeqmode.eq.2) then
            lambdanunu1p=twopihbar*M2nunu*
     +        phdens2(Zcomp,Ncomp,0,0,2,1,gsp,gsn,uunup,edepth,surfwell)
            lambdanunu1h=twopihbar*M2nunu*
     +        phdens2(Zcomp,Ncomp,0,0,1,2,gsp,gsn,uunuh,edepth,surfwell)
            lambdapinu1p=twopihbar*M2pinu*
     +        phdens2(Zcomp,Ncomp,1,0,1,1,gsp,gsn,uupip,edepth,surfwell)
            lambdapinu1h=twopihbar*M2pinu*
     +        phdens2(Zcomp,Ncomp,0,1,1,1,gsp,gsn,uupih,edepth,surfwell)
          else
            do 20 j=1,4
              if (j.eq.1) uu=uunup
              if (j.eq.2) uu=uunuh
              if (j.eq.3) uu=uupip
              if (j.eq.4) uu=uupih
              if (j.gt.2) then
                Zix=1
                Nix=0
                k=2
              else
                Zix=0
                Nix=1
                k=1
              endif
              eopt=max(uu-S(Zix,Nix,k),-20.)
              nen=min(10*numen,int(eopt*10.))
              if (j.le.2) then
                Weff=Wompfac(1)*wvol(k,nen)
              else
                Weff=Wompfac(2)*wvol(k,nen)
              endif
              if (j.eq.1) lambdanunu1p=2.*Weff/hbar
              if (j.eq.2) then
                densh=phdens2(Zcomp,Ncomp,0,0,1,2,gsp,gsn,uu,
     +            edepth,surfwell)
                densp=phdens2(Zcomp,Ncomp,0,0,2,1,gsp,gsn,uu,
     +            edepth,surfwell)
                if (densp.gt.1.) then
                  ratio=densh/densp
                else
                  ratio=1.
                endif
                lambdanunu1h=2.*Weff/hbar*ratio
              endif
              if (j.eq.3) lambdapinu1p=2.*Weff/hbar
              if (j.eq.4) then
                densh=phdens2(Zcomp,Ncomp,0,1,1,1,gsp,gsn,uu,
     +            edepth,surfwell)
                densp=phdens2(Zcomp,Ncomp,1,0,1,1,gsp,gsn,uu,
     +            edepth,surfwell)
                if (densp.gt.1.) then
                  ratio=densh/densp
                else
                  ratio=1.
                endif
                lambdapinu1h=2.*Weff/hbar*ratio
              endif
   20       continue
          endif
          termnunu1p=lambdanunu1p*gsn*dExnup
          termnunu1h=lambdanunu1h*gsn*dExnuh
          termpinu1p=lambdapinu1p*gsp*dExpip
          termpinu1h=lambdapinu1h*gsp*dExpih
          if (flaggshell) then
            damp=ignatyuk(Zcomp,Ncomp,U-uunup,0)/alev(Zcomp,Ncomp)
            gsp=gp(Zcomp,Ncomp)*damp
            gsn=gn(Zcomp,Ncomp)*damp
          endif
          sumnunu1p=sumnunu1p+termnunu1p*
     +      phdens2(Zcomp,Ncomp,ppi,hpi,pnu-1,hnu,gsp,gsn,U-uunup,
     +            edepth,surfwell)
          sumnunu1h=sumnunu1h+termnunu1h*
     +      phdens2(Zcomp,Ncomp,ppi,hpi,pnu,hnu-1,gsp,gsn,U-uunuh,
     +            edepth,surfwell)
          sumpinu1p=sumpinu1p+termpinu1p*
     +      phdens2(Zcomp,Ncomp,ppi-1,hpi,pnu,hnu,gsp,gsn,U-uupip,
     +            edepth,surfwell)
          sumpinu1h=sumpinu1h+termpinu1h*
     +      phdens2(Zcomp,Ncomp,ppi,hpi-1,pnu,hnu,gsp,gsn,U-uupih,
     +            edepth,surfwell)
   10   continue
        if (flaggshell) then
          damp=ignatyuk(Zcomp,Ncomp,U,0)/alev(Zcomp,Ncomp)
          gsp=gp(Zcomp,Ncomp)*damp
          gsn=gn(Zcomp,Ncomp)*damp
        endif
        phtot=phdens2(Zcomp,Ncomp,ppi,hpi,pnu,hnu,gsp,gsn,U,
     +    edepth,surfwell)
        if (phtot.gt.0.) then
          lplus=(sumnunu1p+sumnunu1h+sumpinu1p+sumpinu1h)/phtot
          lambdanuplus=real(lplus)
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
