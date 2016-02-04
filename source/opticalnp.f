      subroutine opticalnp(Zix,Nix,k,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 11, 2013
c | Task  : Optical potential for neutrons and protons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,k,nen,Z,A,i,mw,md,omptype,nrange,nr
      real         eopt,elow,eup,eint,vloc(19),f,v1loc,v2loc,v3loc,
     +             v4loc,Vc,Vcoul,w1loc,w2loc,w3loc,w4loc,d1loc,d2loc,
     +             d3loc,vso1loc,vso2loc,wso1loc,wso2loc,fjoin,V0term,
     +             Vinf,Vterm,en1(numrange),en2(numrange),Dr(numrange),
     +             sr(numrange),factor
c
c ************************ Calculate parameters ************************
c
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c eopt         : incident energy
c ompadjust    : subroutine for local optical model parameter adjustment
c optmodfile   : file with optical model parameters
c optmod       : file with optical model parameters
c numNph       : maximal number of neutrons away from the initial
c                compound nucleus for multiple pre-equilibrium emission
c numZph       : maximal number of protons away from the initial
c                compound nucleus for multiple pre-equilibrium emission
c eomp         : energies on optical model file
c omplines     : number of lines on optical model file
c elow,eup,eint: help variables
c v1adjust..   : adjustable factors for OMP (default 1.)
c vloc         : interpolated optical model parameters
c vomp         : optical model parameters from file
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c Fv1,....     : adjustable factors for OMP (default 1.)
c
c 1. In case of an optical model file, we interpolate between the
c    tabulated values.
c
      call ompadjust(eopt,k)
      optmodfile='                                                     '
      if (Zix.le.numZph.and.Nix.le.numNph) optmodfile=optmod(Zix,Nix,k)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(Zix,Nix,k,1).or.
     +    eopt.gt.eomp(Zix,Nix,k,omplines(Zix,Nix,k))) goto 100
        do 10 nen=1,omplines(Zix,Nix,k)-1
          elow=eomp(Zix,Nix,k,nen)
          eup=eomp(Zix,Nix,k,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 20 i=1,19
              vloc(i)=vomp(Zix,Nix,k,nen,i)+
     +          eint*(vomp(Zix,Nix,k,nen+1,i)-vomp(Zix,Nix,k,nen,i))
   20       continue
            v=Fv1*vloc(1)
            rv=Frv*vloc(2)
            av=Fav*vloc(3)
            w=Fw1*vloc(4)
            rw=Frw*vloc(5)
            aw=Faw*vloc(6)
            vd=Fd1*vloc(7)
            rvd=Frvd*vloc(8)
            avd=Favd*vloc(9)
            wd=Fd1*vloc(10)
            rwd=Frwd*vloc(11)
            awd=Fawd*vloc(12)
            vso=Fvso1*vloc(13)
            rvso=Frvso*vloc(14)
            avso=Favso*vloc(15)
            wso=Fwso1*vloc(16)
            rwso=Frwso*vloc(17)
            awso=Fawso*vloc(18)
            rc=Frc*vloc(19)
            goto 200
          endif
   10   continue
      endif
c
c 2. The general energy-dependent form of the optical potential
c    using parameters per nucleus or the global optical model,
c    both from subroutine omppar.
c
c ZZ,Z         : charge number of residual nucleus
c AA,A         : mass number of residual nucleus
c ompglobal    : flag for use of global optical model
c soukhovitskii: subroutine for global optical model parameters for
c                actinides by Soukhovitskii et al.
c f            : E-Ef
c v1loc.....   : help variables
c enincmax     : maximum incident energy
c ef           : Fermi energy
c v1,v2,v3     : components for V
c w1,w2        : components for W
c d1,d2,d3     : components for Wd
c mw,md        : powers for W and Wd
c vso1,vso2    : components for Vso
c wso1,wso2    : components for Wso
c
  100 Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      if (Z.ge.90.and.ompglobal(Zix,Nix,k)) then
        call soukhovitskii(k,Z,A,eopt)
      else
        f=eopt-ef(Zix,Nix,k)
        rc=Frc*rc0(Zix,Nix,k)
        v1loc=Fv1*v1(Zix,Nix,k)
        v2loc=Fv2*v2(Zix,Nix,k)
        v3loc=Fv3*v3(Zix,Nix,k)
        v4loc=Fv4*7.e-9
c
c Coulomb term for protons
c
        if (k.eq.2.and.ompglobal(Zix,Nix,k)) then
          Vc=1.73/rc*Z/(A**onethird)
          vcoul=Vc*v1loc*(v2loc-2.*v3loc*f+3.*v4loc*f**2)
        else
          vcoul=0.
        endif
c
c Extension up to 1 GeV.
c
c Ejoin     : joining energy for high energy OMP
c fjoin     : help variable
c Vinfadjust: adjustable factor for high energy limit of
c             real central potential
c Vinf      : high energy limit for real central potential
c V0        : V at zero MeV
c Vjoin     : V at joining energy
c Wjoin     : W at joining energy
c
        if (eopt.gt.Ejoin(k)) fjoin=Ejoin(k)-ef(Zix,Nix,k)
        if (eopt.le.Ejoin(k)) then
          v=v1loc*(1.-v2loc*f+v3loc*f**2-v4loc*f**3)+vcoul
        else
          Vinf=-Vinfadjust(k)*30.
          v=Vinf
          V0term=V0(k)-Vinf
          if (V0term.gt.0.) then
            Vterm=(Vjoin(k)-Vinf)/V0term
            if (Vterm.gt.0.) v=Vinf+V0term*exp(f/fjoin*log(Vterm))
ctest     if (k.eq.2.and.ompglobal(Zix,Nix,k))
c    +      vcoul=-Vc*(V0(k)-Vinf)*logterm/fjoin*exp(f/fjoin*logterm)
          endif
        endif
        rv=Frv*rv0(Zix,Nix,k)
        av=Fav*av0(Zix,Nix,k)
        mw=2
        w1loc=Fw1*w1(Zix,Nix,k)
        w2loc=Fw2*w2(Zix,Nix,k)
        if (eopt.le.Ejoin(k)) then
          w=w1loc*f**mw/(f**mw+w2loc**mw)
        else
          w3loc=Fw3*w3(Zix,Nix,k)
          w4loc=Fw4*w4(Zix,Nix,k)
          w=Wjoin(k)-w3loc*fjoin**4/
     +      (fjoin**4+w4loc**4)+w3loc*f**4/(f**4+w4loc**4)
        endif
        rw=Frw*rv0(Zix,Nix,k)
        aw=Faw*av0(Zix,Nix,k)
        vd=0.
        rvd=Frvd*rvd0(Zix,Nix,k)
        avd=Favd*avd0(Zix,Nix,k)
        md=2
        d1loc=Fd1*d1(Zix,Nix,k)
        d2loc=Fd2*d2(Zix,Nix,k)
        d3loc=Fd3*d3(Zix,Nix,k)
        wd=d1loc*f**md*exp(-d2loc*f)/(f**md+d3loc**md)
        rwd=Frwd*rvd0(Zix,Nix,k)
        awd=Fawd*avd0(Zix,Nix,k)
        vso1loc=Fvso1*vso1(Zix,Nix,k)
        vso2loc=Fvso2*vso2(Zix,Nix,k)
        vso=vso1loc*exp(-vso2loc*f)
        rvso=Frvso*rvso0(Zix,Nix,k)
        avso=Favso*avso0(Zix,Nix,k)
        wso1loc=Fwso1*wso1(Zix,Nix,k)
        wso2loc=Fwso2*wso2(Zix,Nix,k)
        wso=wso1loc*f**2/(f**2+wso2loc**2)
        rwso=Frwso*rvso0(Zix,Nix,k)
        awso=Fawso*avso0(Zix,Nix,k)
      endif
c
c Possible additional energy-dependent adjustment of the geometry
c
c ompadjustF       : logical for local OMP adjustment
c ompadjustN,nrange: number of energy ranges for local OMP adjustment
c ompadjustE1,en1  : start energy of local OMP adjustment
c ompadjustE2,en2  : end energy of local OMP adjustment
c ompadjustD,Dr    : depth of local OMP adjustment
c ompadjusts,sr    : variance of local OMP adjustment
c adjustF          : subroutine for local parameter adjustment
c factor           : Woods-Saxon multiplication factor
c
  200 if (ompadjustF(k)) then
        do 210 omptype=1,13
          nrange=ompadjustN(k,omptype)
          do 220 nr=1,nrange
            en1(nr)=ompadjustE1(k,omptype,nr)
            en2(nr)=ompadjustE2(k,omptype,nr)
            Dr(nr)=ompadjustD(k,omptype,nr)
            sr(nr)=ompadjusts(k,omptype,nr)
  220     continue
          call adjustF(eopt,nrange,en1,en2,Dr,sr,factor)
          if (omptype.eq.1) rv=factor*rv
          if (omptype.eq.2) av=factor*av
          if (omptype.eq.3) rw=factor*rw
          if (omptype.eq.4) aw=factor*aw
          if (omptype.eq.5) rvd=factor*rvd
          if (omptype.eq.6) avd=factor*avd
          if (omptype.eq.7) rwd=factor*rwd
          if (omptype.eq.8) awd=factor*awd
          if (omptype.eq.9) rvso=factor*rvso
          if (omptype.eq.10) avso=factor*avso
          if (omptype.eq.11) rwso=factor*rwso
          if (omptype.eq.12) awso=factor*awso
          if (omptype.eq.13) rc=factor*rc
  210   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
