      subroutine energies
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : February 20, 2012
c | Task  : Energies
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          type,nen,i,Zix,Nix,NL,minbins,maxbins
      real             Elast,b
      double precision ma1,ma2,ma,eps
c
c ****************************** Energies ******************************
c
c Center-of-mass energy and wave number
c
c 1. Relativistic case: Ecm = sqrt[ (m+M)**2 + 2ME ] -(m+M)
c                           = 2ME / { sqrt [ (m+M)**2 + 2ME ] + (m+M) }
c
c For incident photons, we always take the relativistic case.
c
c flagrel: flag for relativistic kinematics
c k0     : index of incident particle
c ma1,ma2: help variables
c ma,eps : help variables
c parmass: mass of particle in a.m.u.
c tarmass: mass of target nucleus
c Einc   : incident energy in MeV
c amu    : atomic mass unit in MeV
c eninccm: center-of-mass incident energy in MeV
c wavenum: wave number
c
      if (flagrel.or.k0.eq.0) then
        ma1=parmass(k0)
        ma2=tarmass
        ma=ma1+ma2
        eps=2.*ma2*Einc/amu
        eninccm=real(amu*eps/(sqrt(ma**2+eps)+ma))
        wavenum=real(amu*ma2/hbarc*sqrt(Einc/amu*(Einc/amu+2.*ma1)/
     +    ((ma1+ma2)**2+2.*ma2*Einc/amu)))
      else
c
c 2. Non-relativistic case: Ecm = E*M/(m+M)
c
c specmass: specific mass for target nucleus and all particles
c parZ    : charge number of particle
c parN    : neutron number of particle
c redumass: reduced mass
c hbarc   : hbar.c in MeV.fm
c Etotal  : total energy of compound system (target + projectile)
c S       : separation energy per particle
c targetE : energy of target
c
        eninccm=real(Einc*specmass(parZ(k0),parN(k0),k0))
        wavenum=sqrt(real(2.*amu*redumass(parZ(k0),parN(k0),k0)*
     +    eninccm))/hbarc
      endif
      Etotal=eninccm+S(0,0,k0)+targetE
c
c ***************** Set upper limit for energy grid ********************
c
c eendhigh: last energy point for energy grid for any particle
c eend    : last energy point for energy grid of each particle
c maxen   : total number of energies
c parskip : logical to skip outgoing particle
c egrid   : energies of basic energy grid in MeV
c ebegin  : first energy point of energy grid
c
      eendhigh=0
      do 10 type=0,6
        eend(type)=maxen-1
        if (parskip(type)) goto 10
        do 20 nen=0,maxen
          if (egrid(nen).gt.Etotal-S(0,0,type)) then
            eend(type)=nen
            goto 30
          endif
   20   continue
   30   if (eend(type).gt.ebegin(type))
     +    eend(type)=max(eend(type),ebegin(type)+3)
        eendhigh=max(eendhigh,eend(type))
   10 continue
c
c **** Set limit for cross section spectra ******
c
c speceps: limit for cross section spectra
c xseps  : limit for cross sections
c
      speceps=xseps/eninccm
c
c ********* Set outgoing energies for discrete level scattering ********
c
c nendisc: last discrete bin
c Zindex : charge number index for residual nucleus
c Nindex : neutron number index for residual nucleus
c Nlast  : last discrete level
c numlev2: maximum number of levels
c edis   : energy of level
c Elast  : help variable
c eoutdis: outgoing energy of discrete state reaction
c locate : subroutine to find value in ordered table
c Ebottom: bottom of outgoing energy bin
c
c The variable nendisc is used to determine the last discrete level
c of the population of a residual nucleus. This will be needed for
c the calculation of boundary effects in spectra.
c
      do 110 type=0,6
        nendisc(type)=1
        if (parskip(type)) goto 110
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        NL=Nlast(Zix,Nix,0)
        do 120 i=0,numlev2
          eoutdis(type,i)=Etotal-S(0,0,type)-edis(Zix,Nix,i)
  120   continue
        if (ebegin(type).ge.eend(type)) goto 110
        Elast=eoutdis(type,NL)
        if (Elast.gt.0.) then
          call locate(Ebottom,1,eend(type),Elast,nen)
          nendisc(type)=nen
        endif
  110 continue
c
c ************ Set pre-equilibrium and compound nucleus flags **********
c
c ewfc       : off-set incident energy for width fluctuation calculation
c flagwidth  : flag for width fluctuation calculation
c eurr       : off-set incident energy for URR calculation
c flagurr    : flag for output of unresolved resonance parameters
c flagang    : flag for output of angular distributions
c flagcompang: flag for compound angular distribution calculation
c flagcpang  : flag for compound angular distribution calculation for
c              incident charged particles
c epreeq     : on-set incident energy for preequilibrium calculation
c flagpreeq  : flag for pre-equilibrium calculation
c flaggiant  : flag for collective contribution from giant resonances
c flaggiant0 : flag for collective contribution from giant resonances
c emulpre    : on-set incident energy for multiple preequilibrium
c flagmulpre : flag for multiple pre-equilibrium calculation
c flagendf   : flag for information for ENDF-6 file
c eadd       : on-set incident energy for addition of discrete states
c              to spectra
c eaddel     : on-set incident energy for addition of elastic peak
c              to spectra
c flagadd    : flag for addition of discrete states to spectra
c flagaddel  : flag for addition of elastic peak to spectra
c numZ       : maximal number of protons away from initial compound
c              nucleus
c numN       : maximal number of neutrons away from initial compound
c              nucleus
c mulpreZN   : logical for multiple pre-equilibrium per nucleus
c
      if (Einc.le.ewfc)  then
        flagwidth=.true.
      else
        flagwidth=.false.
      endif
      if (Einc.le.eurr)  then
        flagurr=.true.
      else
        flagurr=.false.
      endif
      if ((k0.eq.1.or.flagcpang).and.flagang.and.Einc.le.50.)  then
        flagcompang=.true.
      else
        flagcompang=.false.
      endif
      if (Einc.lt.epreeq) then
        flagpreeq=.false.
        flaggiant=.false.
      else
        flagpreeq=.true.
        if (flaggiant0) then
          flaggiant=.true.
        else
          flaggiant=.false.
        endif
      endif
      if (Einc.lt.emulpre) then
        flagmulpre=.false.
      else
        flagmulpre=.true.
      endif
      if (Einc.lt.eadd) then
        flagadd=.false.
      else
        flagadd=.true.
      endif
      if (Einc.lt.eaddel) then
        flagaddel=.false.
      else
        flagaddel=.true.
      endif
      do 210 Zix=0,numZ
        do 210 Nix=0,numN
          mulpreZN(Zix,Nix)=.false.
  210 continue
c
c Progressive bin width for excitation energy bins
c
c minbins: minimum number of bins
c maxbins: maximum number of bins
c b      : parameter for slope
c nbins0 : number of continuum excitation energy bins
c nbins  : number of continuum excitation energy bins
c
      if (nbins0.eq.0) then
        minbins=30
        maxbins=numbins
        b=60.
        nbins=minbins+(maxbins-minbins)*Einc*Einc/(Einc*Einc+b*b)
      else
        nbins=nbins0
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
