      subroutine densitypar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : December 13, 2013
c | Task  : Level density parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist,inpalev,inpdeltaW,inpalimit,inpgammald
      character*4      denchar
      character*22     denformat
      character*90     denfile
      integer          Zix,Nix,Z,N,A,ldmod,ia,Nlow0,Ntop0,ibar,imax,
     +                 imin,i,oddZ,oddN
      real             ald0,pshift0,scutoffsys,sigsum,denom,rj,sd,ald,
     +                 Spair,expo,fU,difprev,factor,argum
      double precision mldm,mliquid1,mliquid2
c
c *************************** Initialization ***************************
c
c ldmodel 1: Gilbert and Cameron
c ldmodel 2: Back-shifted Fermi gas
c ldmodel 3: Superfluid model
c ldmodel 4: Statistical HFB model (Goriely)
c ldmodel 5: Combinatorial HFB model (Hilaire and Goriely)
c ldmodel 6: Combinatorial HFB model - Gogny force (Hilaire and Goriely)
c
c Zix : charge number index for residual nucleus
c Nix : neutron number index for residual nucleus
c ZZ,Z: charge number of residual nucleus
c NN,N: neutron number of residual nucleus
c AA,A: mass number of residual nucleus
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      ldmod=ldmodel(Zix,Nix)
c
c
c *********** Read values from level density parameter file ************
c
c denchar   : string for level density file
c denfile   : level density parameter file
c path      : directory containing structure files to be read
c flagcol   : flag for collective enhancement of level density
c ia        : mass number from level density table
c Nlow0,..  : help variables
c ldparexist: flag for existence of tabulated level density parameters
c Nlow      : lowest discrete level for temperature matching
c Ntop      : highest discrete level for temperature matching
c flagasys  : flag for all level density parameters a from systematics
c alev      : level density parameter
c aadjust...: adjustable factors for level density parameters
c             (default 1.)
c Pshift    : adjustable pairing shift
c ctable    : constant to adjust tabulated level densities
c ptable    : constant to adjust tabulated level densities
c
c Level density parameters from the table can always be overruled
c by values given in the input file. With flagasys, all experimental
c level density parameters a from the table can be overruled by the
c systematics.
c We allow a maximum of Ntop=50
c
      denchar='z   '
      write(denchar(2:4),'(i3.3)') Z
      if (ldmod.eq.1)
     +  denfile=path(1:lenpath)//'density/ground/ctm/'//denchar
      if (ldmod.eq.2)
     +  denfile=path(1:lenpath)//'density/ground/bfm/'//denchar
      if (ldmod.eq.3)
     +  denfile=path(1:lenpath)//'density/ground/gsm/'//denchar
      if (ldmod.eq.4)
     +  denfile=path(1:lenpath)//'density/ground/goriely/'//denchar
      if (ldmod.eq.5)
     +  denfile=path(1:lenpath)//'density/ground/hilaire/'//denchar
      if (ldmod.eq.6)
     +  denfile=path(1:lenpath)//'density/ground/hilaireD1M/'//denchar
      inquire (file=denfile,exist=lexist)
      if (.not.lexist) goto 30
      if (flagcol(Zix,Nix).and.ldmod.le.3) then
        denformat='(4x,i4,32x,2i4,2f12.5)'
      else
        denformat='(4x,3i4,2f12.5)'
      endif
      open (unit=2,status='old',file=denfile)
   10 read(2,fmt=denformat,end=30) ia,Nlow0,Ntop0,ald0,pshift0
      if (A.ne.ia) goto 10
      ldparexist(Zix,Nix)=.true.
      if (Nlow(Zix,Nix,0).eq.-1) Nlow(Zix,Nix,0)=Nlow0
      if (Ntop(Zix,Nix,0).eq.-1) Ntop(Zix,Nix,0)=min(Ntop0,50)
      if (.not.flagasys) then
        if (ldmod.le.3) then
          if (alev(Zix,Nix).eq.0.) alev(Zix,Nix)=aadjust(Zix,Nix)*ald0
          do 20 ibar=0,nfisbar(Zix,Nix)
            if (Pshift(Zix,Nix,ibar).eq.1.e-20)
     +        Pshift(Zix,Nix,ibar)=pshift0
   20     continue
        else
          if (ctable(Zix,Nix,0).eq.1.e-20) ctable(Zix,Nix,0)=ald0
          if (ptable(Zix,Nix,0).eq.1.e-20) ptable(Zix,Nix,0)=pshift0
        endif
      endif
   30 close (unit=2)
c
c Matching levels
c
c ibar     : fission barrier
c nfisbar  : number of fission barrier parameters
c Nlast    : last discrete level
c nlev     : number of excited levels for nucleus
c nfistrrot: number of rotational transition states for barrier
c
      do 40 ibar=0,nfisbar(Zix,Nix)
        if (ibar.eq.0) then
          Nlast(Zix,Nix,ibar)=nlev(Zix,Nix)
        else
          Nlast(Zix,Nix,ibar)=max(nfistrrot(Zix,Nix,ibar),1)
        endif
        if (Ntop(Zix,Nix,ibar).eq.-1) Ntop(Zix,Nix,ibar)=
     +    Nlast(Zix,Nix,ibar)
        if (Nlow(Zix,Nix,ibar).eq.-1) Nlow(Zix,Nix,ibar)=2
        if (Ntop(Zix,Nix,ibar).le.2) Nlow(Zix,Nix,ibar)=0
   40 continue
c
c Determine spin cut-off parameter for discrete level region
c
c scutoffsys  : spin cutoff factor for discrete level from systematics
c scutoffdisc : spin cutoff factor for discrete level region
c imin,imax   : help variables
c Ediscrete   : energy of middle of discrete level region
c edis        : energy of level
c efistrrot   : energy of rotational transition states
c sigsum,denom: help variables
c jdis        : spin of the level
c jfistrrot   : spin of rotational transition states
c
c First assign the systematics value, then overrule in case of enough
c discrete level info.
c
      scutoffsys=(0.83*(A**0.26))**2
      do 50 ibar=0,nfisbar(Zix,Nix)
        scutoffdisc(Zix,Nix,ibar)=scutoffsys
        if (ldparexist(Zix,Nix)) then
          imax=Ntop(Zix,Nix,ibar)
          if (ibar.eq.0) then
            imin=Nlow(Zix,Nix,0)
            Ediscrete(Zix,Nix,0)=0.5*(edis(Zix,Nix,imin)+
     +        edis(Zix,Nix,imax))
          else
            imin=1
            Ediscrete(Zix,Nix,ibar)=0.5*(efistrrot(Zix,Nix,ibar,imin)+
     +        efistrrot(Zix,Nix,ibar,imax))
          endif
          sigsum=0.
          denom=0.
          do 60 i=imin,imax
            if (ibar.eq.0) then
              rj=jdis(Zix,Nix,i)
            else
              rj=jfistrrot(Zix,Nix,ibar,i)
            endif
            sigsum=sigsum+rj*(rj+1)*(2*rj+1)
            denom=denom+2*rj+1
   60     continue
          sd=0.
          if (denom.ne.0.) sd=sigsum/(3.*denom)
          if (sd.gt.scutoffsys/3..and.sd.lt.scutoffsys*3.)
     +      scutoffdisc(Zix,Nix,ibar)=sd
        endif
   50 continue
c
c Check input of various level density parameters
c
c inpalev    : logical to determine existence of input value for a
c inpdeltaW  : logical to determine existence of input value for deltaW
c deltaW     : shell correction in nuclear mass
c nucmass    : mass of nucleus
c shellmodel : model for shell correction energy
c mldm       : liquid drop mass
c mliquid1   : function for liquid drop mass (Myers-Swiatecki)
c mliquid2   : function for liquid drop mass (Goriely)
c amu        : atomic mass unit in MeV
c inpalimit  : logical to determine existence of input value for alimit
c alimit     : asymptotic level density parameter
c alphald    : alpha-constant for asymptotic level density parameter
c betald     : beta-constant for asymptotic level density parameter
c twothird   : 2/3
c inpgammald : logical to determine existence of input value for gammald
c gammald    : gamma-constant for asymptotic level density parameter
c gammashell1: gamma-constant for asymptotic level density parameter
c gammashell2: gamma-constant for asymptotic level density parameter
c onethird   : 1/3
c
c shellmodel 1: Myers-Swiatecki
c shellmodel 2: Goriely
c
      if (alev(Zix,Nix).eq.0.) then
        inpalev=.false.
      else
        inpalev=.true.
      endif
      inpdeltaW=.true.
      if (deltaW(Zix,Nix,0).eq.0.) then
        inpdeltaW=.false.
        if (shellmodel.eq.1) then
          mldm=mliquid1(Z,A)
        else
          mldm=mliquid2(Z,A)
        endif
        deltaW(Zix,Nix,0)=real((nucmass(Zix,Nix)-mldm)*amu)
      endif
      inpalimit=.true.
      if (alimit(Zix,Nix).eq.0.) then
        inpalimit=.false.
        alimit(Zix,Nix)=alphald(Zix,Nix)*A+betald(Zix,Nix)*(A**twothird)
      endif
      inpgammald=.true.
      if (gammald(Zix,Nix).eq.-1.) then
        inpgammald=.false.
        gammald(Zix,Nix)=gammashell1(Zix,Nix)/(A**onethird)+gammashell2
      endif
c
c The Ignatyuk formula implies that alev, deltaW, gammald and alimit can
c not all be given as input. In that case we re-determine alev.
c
      if (inpalev.and.inpdeltaW.and.inpalimit.and.inpgammald) then
        inpalev=.false.
        alev(Zix,Nix)=0.
      endif
c
c Pairing corrections
c
c oddZ,oddN     : help variables
c delta0        : systematical pairing energy
c pairconstant  : constant for pairing energy systematics
c pair          : total pairing correction
c Pshiftconstant: global constant for pairing shift
c delta         : energy shift
c
      oddZ=mod(Z,2)
      oddN=mod(N,2)
      delta0(Zix,Nix)=pairconstant/sqrt(real(A))
c
c Defaults
c
      if (pair(Zix,Nix).eq.1.e-20) then
        if (ldmod.eq.3) then
          pair(Zix,Nix)=(oddZ+oddN)*delta0(Zix,Nix)
        else
          if (ldmod.eq.2) then
            pair(Zix,Nix)=(1.-oddZ-oddN)*delta0(Zix,Nix)
          else
            pair(Zix,Nix)=(2.-oddZ-oddN)*delta0(Zix,Nix)
          endif
        endif
      endif
      do 70 ibar=0,nfisbar(Zix,Nix)
        if (Pshift(Zix,Nix,ibar).eq.1.e-20)
     +    Pshift(Zix,Nix,ibar)=Pshiftconstant(Zix,Nix)
 70   continue
c
c Generalized superfluid model. The critical functions are
c calculated here.
c
c Tcrit    : critical temperature
c factor   : help variable
c aldcrit  : critical level density parameter
c S        : separation energy per particle
c fU,factor: help variables
c Econd    : condensation energy
c Ucrit    : critical U
c Scrit    : critical entropy
c Dcrit    : critical determinant
c
      if (ldmod.eq.3) then
        Tcrit(Zix,Nix)=0.567*delta0(Zix,Nix)
        ald=alimit(Zix,Nix)
        difprev=0.
   80   factor=(1.-exp(-gammald(Zix,Nix)*ald*Tcrit(Zix,Nix)**2))/
     +    (ald*(Tcrit(Zix,Nix)**2))
        aldcrit(Zix,Nix)=alimit(Zix,Nix)*(1.+deltaW(Zix,Nix,0)*factor)
        if (abs(aldcrit(Zix,Nix)-ald).gt.0.001.and.
     +    abs(aldcrit(Zix,Nix)-ald).ne.difprev) then
          difprev=abs(aldcrit(Zix,Nix)-ald)
          ald=aldcrit(Zix,Nix)
          if (ald.gt.1.) goto 80
        endif
        if (aldcrit(Zix,Nix).lt.alimit(Zix,Nix)/3.) then
          expo=min(-gammald(Zix,Nix)*S(Zix,Nix,1),80.)
          fU=1.-exp(expo)
          factor=1.+fU*deltaW(Zix,Nix,0)/S(Zix,Nix,1)
          aldcrit(Zix,Nix)=max(alimit(Zix,Nix)*factor,1.)
        endif
        Econd(Zix,Nix)=1.5/pi2*aldcrit(Zix,Nix)*delta0(Zix,Nix)**2
        Ucrit(Zix,Nix)=aldcrit(Zix,Nix)*Tcrit(Zix,Nix)**2+Econd(Zix,Nix)
        Scrit(Zix,Nix)=2.*aldcrit(Zix,Nix)*Tcrit(Zix,Nix)
        Dcrit(Zix,Nix)=144./pi*(aldcrit(Zix,Nix)**3)*(Tcrit(Zix,Nix)**5)
        do 90 ibar=0,nfisbar(Zix,Nix)
          delta(Zix,Nix,ibar)=Econd(Zix,Nix)-pair(Zix,Nix)-
     +      Pshift(Zix,Nix,ibar)
  90    continue
      else
c
c Constant temperature and back-shifted Fermi gas model
c
        do 100 ibar=0,nfisbar(Zix,Nix)
          delta(Zix,Nix,ibar)=pair(Zix,Nix)+Pshift(Zix,Nix,ibar)
 100    continue
      endif
c
c 1. If no experimental level density parameter is available,
c    i.e. as determined from the neutron resonance spacing, use the
c    Ignatyuk formula to derive the level density parameter at the
c    separation energy.
c
c Spair    : help variable
c
      Spair=S(Zix,Nix,1)-delta(Zix,Nix,0)
      Spair=max(Spair,1.)
      if (.not.inpalev) then
        fU=1.-exp(-gammald(Zix,Nix)*Spair)
        factor=1.+fU*deltaW(Zix,Nix,0)/Spair
        alev(Zix,Nix)=aadjust(Zix,Nix)*alimit(Zix,Nix)*factor
        alev(Zix,Nix)=max(alev(Zix,Nix),1.)
      else
c
c 2. If an experimental level density parameter is available, then we
c    impose the extra boundary boundary condition that it should be
c    equal to the energy dependent level density parameter at the
c    neutron separation energy. There are various possibilities.
c    If alimit is not given as input, we re-adjust it.
c
        if (ldmod.ne.3) then
          if (.not.inpalimit) then
            fU=1.-exp(-gammald(Zix,Nix)*Spair)
            factor=1.+fU*deltaW(Zix,Nix,0)/Spair
            alimit(Zix,Nix)=alev(Zix,Nix)/factor
          else
c
c If both alev and alimit are explicitly given as input, we
c re-adjust deltaW, provided it is not given as input.
c
            if (.not.inpdeltaW) then
              fU=1.-exp(-gammald(Zix,Nix)*Spair)
              factor=alev(Zix,Nix)/alimit(Zix,Nix)-1.
              deltaW(Zix,Nix,0)=Spair*factor/fU
            else
c
c Determine gammald if alev, alimit and deltaW are given by input.
c
c argum: help variable
c
              argum=1.-Spair/deltaW(Zix,Nix,0)*
     +          (alev(Zix,Nix)/alimit(Zix,Nix)-1.)
              if (argum.gt.0..and.argum.lt.1.) then
                gammald(Zix,Nix)=-1./Spair*log(argum)
              else
c
c If gammald can not be solved or is unphysical (this may happen for
c certain parameter combinations) we re-adjust the shell correction.
c
                fU=1.-exp(-gammald(Zix,Nix)*Spair)
                factor=alev(Zix,Nix)/alimit(Zix,Nix)-1.
                deltaW(Zix,Nix,0)=Spair*factor/fU
              endif
            endif
          endif
        endif
      endif
c
c ************************** Fission ***********************************
c
c Determine deltaW on the fission barrier.
c Note that this overrules input parameters for deltaW.
c
c flagfission : flag for fission
c flagcolldamp: flag for damping of collective effects in effective
c               level density (without explicit collective enhancement)
c               Only used for Bruyeres-le-Chatel (Pascal Romain) fission
c               model
c axtype      : type of axiality of barrier
c                  1: axial symmetry
c                  2: left-right asymmetry
c                  3: triaxial and left-right symmetry
c                  4: triaxial no left-right symmetry
c                  5: no symmetry
c
      if (flagfission) then
        do 110 ibar=1,nfisbar(Zix,Nix)
          if (deltaW(Zix,Nix,ibar).eq.0.) then
            if (flagcolldamp) then
              deltaW(Zix,Nix,ibar)=abs(deltaW(Zix,Nix,0))*twothird
            else
              if (ibar.eq.1) then
                if (axtype(Zix,Nix,1).eq.1) then
                  deltaW(Zix,Nix,ibar)=1.5
                else
                  deltaW(Zix,Nix,ibar)=2.5
                endif
              else
                deltaW(Zix,Nix,ibar)=0.6
              endif
            endif
          endif
  110   continue
        if (nfisbar(Zix,Nix).eq.1.and.fbarrier(Zix,Nix,1).eq.0.)
     +    deltaW(Zix,Nix,1)=deltaW(Zix,Nix,2)
      endif
c
c ************** Single-particle level density parameter g *************
c
c g  : single-particle level density parameter
c pi2: pi**2
c Kph: constant for single-particle level density parameter (g=A/Kph)
c gp : single-particle proton level density parameter
c gn : single-particle neutron level density parameter
c
c One component
c
      if (g(Zix,Nix).eq.0.) g(Zix,Nix)=A/Kph
c
c Two component
c
      if (gp(Zix,Nix).eq.0.) gp(Zix,Nix)=Z/Kph
      if (gn(Zix,Nix).eq.0.) gn(Zix,Nix)=N/Kph
      gn(Zix,Nix)=gnadjust(Zix,Nix)*gn(Zix,Nix)
      gp(Zix,Nix)=gpadjust(Zix,Nix)*gp(Zix,Nix)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
