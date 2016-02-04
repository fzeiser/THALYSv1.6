      subroutine opticalcomp(Zix,Nix,k,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 11, 2013
c | Task  : Optical potential for composite particles
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,k,nen,i,iz,in,ia,omptype,nrange,nr
      real         eopt,elow,eup,eint,vloc(19),e,vn,rvn,avn,wn,rwn,awn,
     +             vdn,rvdn,avdn,wdn,vp,rvp,avp,wp,rwp,awp,vdp,rvdp,
     +             avdp,wdp,vson,rvson,avson,wson,vsop,rvsop,avsop,wsop,
     +             en1(numrange),en2(numrange),Dr(numrange),
     +             sr(numrange),factor
c
c ******************* Optical model input file *************************
c
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c eopt         : incident energy
c numNph       : maximal number of neutrons away from the initial
c                compound nucleus for multiple pre-equilibrium emission
c numZph       : maximal number of protons away from the initial
c                compound nucleus for multiple pre-equilibrium emission
c optmodfile   : file with optical model parameters
c optmod       : file with optical model parameters
c eomp         : energies on optical model file
c omplines     : number of lines on optical model file
c elow,eup,eint: help variables
c vloc         : interpolated optical model parameters
c vomp         : optical model parameters from file
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c Fv1,...      : factor
c
c 1. In case of an optical model file, we interpolate between the
c    tabulated values
c
      optmodfile='                                                     '
      if (Zix.le.numZph.and.Nix.le.numNph) optmodfile=optmod(Zix,Nix,k)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(Zix,Nix,k,1).or.
     +    eopt.gt.eomp(Zix,Nix,k,omplines(Zix,Nix,k))) goto 100
        call ompadjust(eopt,k)
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
c    using parameters per nucleus or a global optical model,
c    both from subroutine omppar.
c
c We use the Watanabe method (S. Watanabe, Nucl. Phys. 8, 484 (1958),
c see also D.G. Madland, Proceedings of a Specialists' Meeting on
c preequilibrium nuclear reactions, Semmering, Austria, February 10-12
c 1988, p. 103) to make a composite particle potential out of the
c proton and neutron potential. The simplified formula is:
c                V(d,E)=z*V(n,E/a)+n*V(p,E/a)
c with z, n and a, the proton, neutron and mass number of the particle,
c and similarly for Wd and W. VSO and Wso are as in the nucleon
c potentials.
c We take a similar weighting for the geometry parameters.
c
c 1. Neutron potential.
c
c e        : incident energy (help variable)
c parA     : mass number of particle
c opticalnp: subroutine for optical potential for neutrons and protons
c vn,wn....: optical model parameters for neutrons
c
  100 e=eopt/parA(k)
      call opticalnp(Zix,Nix,1,e)
      vn=v
      rvn=rv
      avn=av
      wn=w
      rwn=rw
      awn=aw
      vdn=vd
      rvdn=rvd
      avdn=avd
      wdn=wd
c
c 2. Proton potential.
c
c vp,wp....: optical model parameters for protons
c
      call opticalnp(Zix,Nix,2,e)
      vp=v
      rvp=rv
      avp=av
      wp=w
      rwp=rw
      awp=aw
      vdp=vd
      rvdp=rvd
      avdp=avd
      wdp=wd
c
c 3. Another 2 calls to opticalnp for eopt=E to determine Vso and Wso.
c
c vson,wson: Vso and Wso for neutrons
c vsop,wsop: Vso and Wso for protons
c
      call opticalnp(Zix,Nix,1,eopt)
      vson=vso
      rvson=rvso
      avson=avso
      wson=wso
      call opticalnp(Zix,Nix,2,eopt)
      vsop=vso
      rvsop=rvso
      avsop=avso
      wsop=wso
c
c 4. Final potential depths: construct V, W, Wd, Vso and Wso.
c
c parZ  : charge number of particle
c parN  : neutron number of particle
c Fv1,..: adjustable factors for OMP (default 1.)
c
      call ompadjust(eopt,k)
      iz=parZ(k)
      in=parN(k)
      ia=parA(k)
      v=Fv1*(in*vn+iz*vp)
      rv=Frv*(in*rvn+iz*rvp)/ia
      av=Fav*(in*avn+iz*avp)/ia
      w=Fw1*(in*wn+iz*wp)
      rw=Frw*(in*rwn+iz*rwp)/ia
      aw=Faw*(in*awn+iz*awp)/ia
      vd=Fd1*(in*vdn+iz*vdp)
      rvd=Frvd*(in*rvdn+iz*rvdp)/ia
      avd=Favd*(in*avdn+iz*avdp)/ia
      wd=Fd1*(in*wdn+iz*wdp)
      rwd=Frwd*(in*rvdn+iz*rvdp)/ia
      awd=Fawd*(in*avdn+iz*avdp)/ia
      rvso=Frvso*(in*rvson+iz*rvsop)/ia
      avso=Favso*(in*avson+iz*avsop)/ia
      rwso=Frwso*(in*rvson+iz*rvsop)/ia
      awso=Fawso*(in*avson+iz*avsop)/ia
      rc=Frc*rc
      if (k.eq.3) then
        vso=Fvso1*(vson+vsop)/2.
        wso=Fwso1*(wson+wsop)/2.
      endif
      if (k.eq.4.or.k.eq.5) then
        vso=Fvso1*(vson+vsop)/6.
        wso=Fwso1*(wson+wsop)/6.
      endif
      if (k.eq.6) then
        vso=0.
        wso=0.
      endif
c
c Alternative options for complex particle OMP
c
c deuteronomp: deuteron optical model
c alphaomp   : alpha optical model
c
      if (k.eq.3.and.deuteronomp.ge.2)
     +  call opticaldeut(Zix,Nix,eopt)
      if (k.eq.6.and.alphaomp.eq.2) call opticalalpha
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
