      subroutine colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib,Kcoll)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn, Arjan Koning, Stephane Hilaire,
c |         Stephane Goriely and Pascal Romain
c | Date  : April 4, 2012
c | Task  : Collective enhancement
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      character*80 key
      integer      Zix,Nix,ibar,A,l
      real         Eex,ald,Krot,Kvib,Kcoll,Krot0,Kvib0,U,aldgs,ignatyuk,
     +             Temp,Kr,factor,rotfactor,aldint,damprot,avib,Pvib,
     +             Tvib,deltaS,deltaU,Cvib,term,omegavib(3),gammavib,
     +             nvib,damper,expo0,expo,spincutgs,spincutbf,spincut,
     +             aldf
c
c ************************ Collective enhancement **********************
c
c Zix      : charge number index for residual nucleus
c Nix      : neutron number index for residual nucleus
c Eex      : excitation energy
c ald,aldgs: level density parameter
c ibar     : fission barrier number, zero for states on ground state
c Krot     : rotational enhancement factor
c Krot0    : rotational enhancement factor (undamped)
c Kvib     : vibrational enhancement
c Kvib0    : vibrational enhancement factor (undamped)
c Kcoll    : total collective enhancement
c flagcol  : flag for collective enhancement of level density
c ldmodel  : level density model
c Exmatch  : matching energy for Ex
c U        : excitation energy minus pairing energy
c delta    : energy shift
c pair     : total pairing correction
c ignatyuk : function for energy dependent level density parameter a
c Temp     : nuclear temperature
c adjust   : subroutine for energy-dependent parameter adjustment
c factor   : multiplication factor
c ldadjust : logical for energy-dependent level density adjustment
c
      Krot=1.
      Krot0=1.
      Kvib=1.
      Kvib0=1.
      Kcoll=1.
      if (.not.flagcol(Zix,Nix)) return
      if (ldmodel(Zix,Nix).eq.1.and.Eex.lt.Exmatch(Zix,Nix,ibar)) return
      U=Eex-delta(Zix,Nix,ibar)
      if (U.le.0.) return
      if (ibar.eq.0) then
        aldgs=ald
      else
        aldgs=ignatyuk(Zix,Nix,Eex,0)
      endif
      Temp=sqrt(U/aldgs)
      if (ldadjust(Zix,Nix)) then
        key='krotconstant'
        call adjust(Eex,key,Zix,Nix,0,ibar,factor)
        Kr=factor*Krotconstant(Zix,Nix,ibar)
      else
        Kr=Krotconstant(Zix,Nix,ibar)
      endif
c
c 1. Specific collective enhancement for Bruyeres-le-Chatel
c    (Pascal Romain) fission model
c
c fismodel    : fission model
c flagcolldamp: flag for damping of collective effects in effective
c               level density (without explicit collective enhancement)
c               Only used for Bruyeres-le-Chatel (Pascal Romain) fission
c               model
c axtype      : type of axiality of barrier
c               1: axial symmetry
c               2: left-right asymmetry
c               3: triaxial and left-right symmetry
c               4: triaxial no left-right symmetry
c               5: no symmetry
c Krotconstant: normalization constant for rotational enhancement
c damprot     : damping factor for rotational enhancement
c rotfactor   : rotational enhancement factor
c aldint      : level density parameter
c
      if (fismodel.eq.1.and.flagcolldamp) then
        if (ibar.eq.0) return
c
c Tri-axial barrier
c
        if (axtype(Zix,Nix,ibar).eq.2) then
          aldint=ald*8./13.5
          rotfactor=Kr*(U/aldint)**(0.25)
        else
c
c Axial and other barriers
c
          rotfactor=max(Kr,1.)
        endif
c
c Damping of enhancement
c
        damprot=1./(1.+exp(-0.5*(U-18.)))
        Krot=rotfactor*(1.-damprot)+damprot
      else
c
c 2. Default calculation
c
c Calculation of Kvib
c
c AA,A     : mass number of residual nucleus
c avib     : level density parameter for vibrational model
c pair,Pvib: total pairing correction
c Tvib     : temperature for vibrational model
c kvibmodel: model for vibrational enhancement
c expo     : exponent
c twothird : 2/3
c deltaS   : entropy change
c deltaU   : excitation energy change
c Cvib     : constant for vibrational enhancement
c onethird : 1/3
c term     : help variable
c deltaW   : shell correction in nuclear mass
c omegavib : energy of vibrational excitation
c l        : multipolarity
c nvib     : occupation number
c pi2      : pi**2
c
        A=AA(Zix,Nix,0)
c
c Kvibmodel 1:  Liquid drop
c
        avib=A/13.
        Pvib=pair(Zix,Nix)
        if (Eex.gt.Pvib) then
          Tvib=sqrt((Eex-Pvib)/avib)
        else
          Tvib=0.
        endif
        if (Tvib.gt.0.) then
          if (kvibmodel.eq.1) then
            expo=min(0.0555*(A**twothird)*(Tvib**(4./3.)),80.)
            Kvib0=exp(expo)
          else
c
c Kvibmodel 2: Bose gas
c
c In this case, Kvib is automatically damped at high energies
c
            deltaS=0.
            deltaU=0.
            Cvib=0.0075*(A**onethird)
            term=A**(-5./6.)/(1.+0.05*deltaW(Zix,Nix,ibar))
            omegavib(2)=65.*term
            omegavib(3)=100.*term
            do 10 l=2,3
              gammavib=Cvib*(omegavib(l)**2+4.*pi2*Tvib**2)
              expo0=gammavib/(2.*omegavib(l))
              expo=omegavib(l)/Tvib
              if (expo0.le.80.and.expo.gt.0..and.expo.le.80.) then
                nvib=exp(-expo0)/(exp(expo)-1.)
                deltaS=deltaS+(2*l+1)*
     +            ((1.+nvib)*log(1.+nvib)-nvib*log(nvib))
                deltaU=deltaU+(2*l+1)*omegavib(l)*nvib
              endif
   10       continue
            Kvib0=exp(deltaS-deltaU/Tvib)
          endif
        endif
c
c Calculation of damping function and Krot
c
c damper      : energy damping function
c Ufermi      : energy of Fermi distribution for damping of ground-state
c             : rotational effects
c cfermi      : width of Fermi distribution for damping of ground-state
c               rotational effects
c spincutgs   : spin-cutoff parameter squared (perpendicular projection)
c               for ground state
c Irigid      : rigid body value of moment of inertia
c
c Ground state
c
        if (ibar.eq.0) then
          damper=0.
          expo=(U-Ufermi)/cfermi
          if (expo.le.80.) damper=1./(1.+exp(expo))
          spincutgs=Kr*Irigid(Zix,Nix,0)*Temp
          Krot0=max(spincutgs,1.)
c
c Fission barrier
c
c Ufermibf : energy of Fermi distribution for damping of barrier
c          : rotational effects
c cfermibf : width of Fermi distribution for damping of barrier
c          : rotational effects
c spincutbf: spin-cutoff parameter squared (perpendicular projection)
c            for fission barrier
c twopi    : 2.*pi
c spincut  : spin cutoff factor
c
        else
          damper=0.
          expo=(U-Ufermibf)/cfermibf
          if (expo.le.80.) damper=1./(1.+exp(expo))
          spincutbf=Krotconstant(Zix,Nix,ibar)*Irigid(Zix,Nix,ibar)*
     +      Temp
          if (axtype(Zix,Nix,ibar).eq.1) Krot0=spincutbf
          if (axtype(Zix,Nix,ibar).eq.2) Krot0=2.*spincutbf
          if (axtype(Zix,Nix,ibar).ge.3) then
            aldf=ignatyuk(Zix,Nix,Eex,ibar)
            term=spincutbf*sqrt(spincut(Zix,Nix,aldf,Eex,ibar)*
     +        (1.-twothird*abs(beta2(Zix,Nix,ibar))))
            if (axtype(Zix,Nix,ibar).eq.3) Krot0=0.5*sqrt(twopi)*term
            if (axtype(Zix,Nix,ibar).eq.4) Krot0=sqrt(twopi)*term
            if (axtype(Zix,Nix,ibar).eq.5) Krot0=2.*sqrt(twopi)*term
          endif
          Krot0=max(Krot0,1.)
        endif
        Krot=1.+(Krot0-1.)*damper
        if (kvibmodel.eq.1) then
          Kvib=1.+(Kvib0-1.)*damper
        else
          Kvib=Kvib0
        endif
      endif
      Kcoll=max(Krot*Kvib,1.)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
