      subroutine densitymatch(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : April 4, 2012
c | Task  : Level density matching solution
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,ibar,i,A,nEx,j,Nstart
      real             ald,Tm,Exm,E0m,ignatyuk,logrholoc(-1:1),
     +                 Exend,dEx,Eex,U,Krot,Kvib,Kcoll,P,logrhomatch
      double precision fermi,rhomatch
c
c ************* Determine level density matching parameters ************
c
c For non-fissile nuclei, there are no fission barriers and only the
c loop for ibar=0 is performed, i.e. for level densities relative to the
c ground state.
c
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c nfisbar   : number of fission barrier parameters
c Ntop,NP   : highest discrete level for temperature matching
c Nlow,NLo  : lowest discrete level for temperature matching
c EL,EP     : matching level energies
c edis      : energy of level
c efistrrot : energy of rotational transition states
c
      do 10 ibar=0,nfisbar(Zix,Nix)
c
c 1. Determine matching levels and energies.
c
        NP=Ntop(Zix,Nix,ibar)
        NLo=Nlow(Zix,Nix,ibar)
c
c A. Parameters for the ground state level density
c
        if (ibar.eq.0) then
          EL=edis(Zix,Nix,NLo)
          EP=edis(Zix,Nix,NP)
        else
c
c B. Parameters for the fission level density
c
          EL=efistrrot(Zix,Nix,ibar,NLo)
          EP=efistrrot(Zix,Nix,ibar,NP)
        endif
        if (ldmodel(Zix,Nix).eq.2.or.ldmodel(Zix,Nix).eq.3.or.
     +    ldexist(Zix,Nix,ibar)) goto 10
c
c 2. Solve level density matching problem for Gilbert-Cameron model
c
c pair,P    : total pairing correction
c AA,A      : mass number of residual nucleus
c Exend     : end of possible energy region
c nEx       : number of energy points for T matching
c nummatchT : maximum number of energy points for T matching
c logrho    : logarithm of level density
c temprho   : temperature
c Eex       : excitation energy
c ald       : level density parameter
c ignatyuk  : function for energy dependent level density parameter a
c colenhance: subroutine for collective enhancement
c Kvib      : vibrational enhancement factor
c Krot      : rotational enhancement factor
c Kcoll     : total collective enhancement
c fermi     : function for Fermi gas level density formula
c Nstart    : energy starting point from which T starts to behave
c             smoothly (for interpolation)
c
c Calculate logarithm of level density and its derivative (temperature)
c on excitation energy grid.
c
        P=pair(Zix,Nix)
        A=AA(Zix,Nix,0)
        Exend=20.+300./A
        dEx=0.1
        nEx=int(Exend/dEx)
        do 20 i=1,nummatchT
          logrho(i)=0.
          temprho(i)=0.
   20   continue
        do 30 i=nEx,1,-1
          do 40 j=-1,1
            Eex=dEx*(i+0.5*j)
            U=Eex-P
            if (U.gt.0.) then
              ald=ignatyuk(Zix,Nix,Eex,ibar)
              call colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib,Kcoll)
              logrholoc(j)=real(log(Kcoll*
     +          fermi(Zix,Nix,ald,Eex,P,ibar)))
            else
              logrholoc(j)=0.
            endif
   40     continue
          logrho(i)=logrholoc(0)
          if (logrholoc(1).ne.logrholoc(-1))
     +      temprho(i)=dEx/(logrholoc(1)-logrholoc(-1))
          if (temprho(i).le.0.1) temprho(i)=temprho(i+1)
   30   continue
        Nstart=1
        do 50 i=nEx,1,-1
          if (i.lt.nEx.and.temprho(i).ge.temprho(i+1)) then
            Nstart=i+1
            goto 60
          endif
   50   continue
c
c T,Tm       : nuclear temperature
c Exmatch,Exm: matching point for Ex
c E0,E0m     : constant of temperature formula
c E0save     : E0 value saved for matching routine
c flagcol    : flag for collective enhancement of level density
c Tmemp      : empirical estimate for temperature
c gammald    : gamma-constant for asymptotic level density parameter
c deltaW     : shell correction in nuclear mass
c Exmemp     : empirical estimate for matching point for Ex
c ldparexist : flag for existence of tabulated level density
c              parameters
c flagctmglob: flag for global CTM model (no discrete level info)
c matching   : subroutine to determine matching between temperature
c              and Fermi-gas region
c pol1       : subroutine for interpolation of first order
c rhomatch   : level density at matching point
c Tadjust....: adjustable factors for level density parameters
c
   60   Tm=T(Zix,Nix,ibar)
        Exm=Exmatchadjust(Zix,Nix,ibar)*Exmatch(Zix,Nix,ibar)
        E0m=E0(Zix,Nix,ibar)
        E0save=E0m
c
c Empirical estimates needed in case of trouble or no discrete levels.
c In those cases, the empirical formula for T is the starting point.
c
        if (flagcol(Zix,Nix)) then
          Tmemp=-0.22+
     +      9.4/sqrt(max(A*(1.+gammald(Zix,Nix)*deltaW(Zix,Nix,0)),1.))
          Exmemp=2.67+253./A+P
        else
          Tmemp=-0.25+
     +      10.2/sqrt(max(A*(1.+gammald(Zix,Nix)*deltaW(Zix,Nix,0)),1.))
          Exmemp=2.33+253./A+P
        endif
        Tmemp=max(Tmemp,0.1)
        Exmemp=max(Exmemp,0.1)
c
c Normal case: CTM parameters are derived from matching automatically,
c i.e. T an Exmatch are not given in the input file.
c
        if (Tm.eq.0..and.Exm.eq.0.) then
c
c A. Discrete levels given
c
          if (ldparexist(Zix,Nix).and..not.flagctmglob) then
            call matching(Zix,Nix,Exm,ibar)
c
c If Exm was set to 0 in matching.f, it means that no reasonable value
c was found. In that case we first use an empirical value for T.
c
            if (Exm.gt.0.) then
              i=int(Exm/dEx)
              call pol1(i*dEx,(i+1)*dEx,temprho(i),temprho(i+1),Exm,Tm)
            else
              Tm=Tmemp
              call locate(temprho,Nstart,nEx-1,Tm,i)
              if (i.gt.0.and.i.le.nEx-1) call pol1(temprho(i),
     +          temprho(i+1),i*dEx,(i+1)*dEx,Tm,Exm)
            endif
          else
c
c B. No discrete levels given and/or global CTM model
c
            Tm=Tmemp
            call locate(temprho,Nstart,nEx-1,Tm,i)
            if (i.gt.0)
     +        call pol1(temprho(i),temprho(i+1),i*dEx,(i+1)*dEx,Tm,Exm)
          endif
c
c If Exm is still unrealistic, we first use an empirical value for T.
c
          ald=ignatyuk(Zix,Nix,Exm,ibar)
          if (Exm.le.max(2.25/ald+P,0.)+0.11) Exm=0.
          if (Exm.gt.3.*Exmemp) Exm=0.
        endif
c
c Special case: either T or Exmatch is given in the input
c If the Exmatch given in the input is unphysical, we choose an
c empirical value.
c
        if (Exm.eq.0.) then
          if (Tm.eq.0.) Tm=Tmemp
          call locate(temprho,Nstart,nEx-1,Tm,i)
          if (i.gt.0)
     +      call pol1(temprho(i),temprho(i+1),i*dEx,(i+1)*dEx,Tm,Exm)
          ald=ignatyuk(Zix,Nix,Exm,ibar)
          if (Exm.le.max(2.25/ald+P,0.)+0.11) Exm=Exmemp
          if (Exm.eq.0.) Exm=Exmemp
          if (Exm.gt.3.*Exmemp) Exm=Exmemp
        endif
        if (Tm.eq.0.) then
          ald=ignatyuk(Zix,Nix,Exm,ibar)
          if (Exm.le.max(2.25/ald+P,0.)+0.11) Exm=Exmemp
          if (Exm.gt.3.*Exmemp) Exm=Exmemp
          i=int(Exm/dEx)
          if (i.gt.0)
     +      call pol1(i*dEx,(i+1)*dEx,temprho(i),temprho(i+1),Exm,Tm)
        endif
        if (E0m.eq.1.e-20) then
          i=int(Exm/dEx)
          if (i.gt.0)
     +      call pol1(i*dEx,(i+1)*dEx,logrho(i),logrho(i+1),Exm,
     +      logrhomatch)
          rhomatch=exp(dble(logrhomatch))
          E0m=Exm-Tm*real(log(dble(Tm)*rhomatch))
        endif
        if (Tm.eq.0.) Tm=Tmemp
c
c Possible iteration after input defined adjustment
c
        if (T(Zix,Nix,ibar).eq.0..and.Tadjust(Zix,Nix,ibar).ne.1.) then
          T(Zix,Nix,ibar)=Tadjust(Zix,Nix,ibar)*Tm
          goto 60
        else
          T(Zix,Nix,ibar)=Tm
        endif
        if (E0(Zix,Nix,ibar).eq.1.e-20.and.E0adjust(Zix,Nix,ibar).ne.1.)
     +    then
          E0(Zix,Nix,ibar)=E0adjust(Zix,Nix,ibar)*E0m
          goto 60
        else
          E0(Zix,Nix,ibar)=E0m
        endif
        Exmatch(Zix,Nix,ibar)=Exm
   10 continue
c
c Set theoretical value of D0
c
c dtheory: subroutine for theoretical calculation of average neutron
c          spacings
c D0theo : mean s-wave resonance spacing
c Dl     : mean s-wave resonance spacing per l value
c
      call dtheory(Zix,Nix,0.)
      D0theo(Zix,Nix)=Dl(0)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
