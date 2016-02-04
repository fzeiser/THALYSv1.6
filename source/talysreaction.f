      subroutine talysreaction
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 2, 2013
c | Task  : Reaction models
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ************************* Reaction calculation ***********************
c
c Basic cross sections (optical model) and initialisation
c
c flagreaction: flag for calculation of nuclear reactions
c flagompall  : flag for new optical model calculation for all residual
c               nuclei
c basicxs     : subroutine for basic cross sections and transmission
c               coefficients
c parinclude  : logical to include outgoing particle
c gamma       : subroutine for gamma cross section and transmission
c               coefficients
c enincmax    : maximum incident energy
c epreeq      : on-set incident energy for preequilibrium calculation
c flagracap   : flag for radiative capture model
c preeqinit   : subroutine for initialization of general
c               pre-equilibrium parameters
c excitoninit : subroutine for initialization of exciton model
c               parameters
c racapinit   : subroutine for initialization of radiative capture model
c flagcomp    : flag for compound nucleus calculation
c compoundinit: subroutine for initialization of compound model
c               parameters
c flagastro   : flag for calculation of astrophysics reaction rate
c astroinit   : subroutine for initialization of astrophysics quantities
c
      if (.not.flagreaction) goto 100
      if (.not.flagompall) call basicxs(0,0)
      if (parinclude(0)) call gamma(0,0)
      if (enincmax.ge.epreeq.or.flagracap) then
        call preeqinit
        call excitoninit
      endif
      if (flagracap) call racapinit
      if (flagcomp) call compoundinit
      if (flagastro) call astroinit
c
c Loop over incident energies
c
c Initialisation
c
c nin        : counter for incident energies
c numinc     : number of incident energies
c eninc,Einc : incident energy in MeV
c energies   : subroutine for energies
c reacinitial: subroutine for initialization of arrays for various
c              cross sections
c eninclow   : minimal incident energy for nuclear model calculations
c
      do 10 nin=1,numinc
        if (flagompall) call basicxs(0,0)
        Einc=eninc(nin)
        call energies
        call reacinitial
        if (Einc.lt.eninclow) goto 10
c
c Optical model
c
c incident   : subroutine for main settings and basic cross sections for
c              incident energy
c exgrid     : subroutine to set excitation energy grid
c flagrecoil : flag for calculation of recoils
c recoilinit : subroutine for calculation of initial recoil velocity
c              and direction
c lmaxinc    : maximal l-value for transmission coefficients
c xsreacinc  : reaction cross section for incident channel
c xseps      : limit for cross sections
c flaginitpop: flag for initial population distribution
c
        call incident
        call exgrid(0,0)
        if (flagrecoil) call recoilinit
c
c In certain cases, there will be no nuclear reaction calculation
c
        if ((lmaxinc.eq.-1.or.xsreacinc.lt.xseps).and..not.flaginitpop)
     +    goto 20
c
c Direct reactions
c
c direct   : subroutine for calculation of direct inelastic cross
c            sections
c flagracap: flag for radiative capture model
c racap    : subroutine for radiative capture model
c
        call direct
        if (flagracap) call racap
c
c Pre-equilibrium reactions
c
c flagpreeq : flag for pre-equilibrium calculation
c preeq     : subroutine for preequilibrium reactions
c population: subroutine for processing of pre-equilibrium spectra
c             into population bins
c
        if (flagpreeq) then
          call preeq
          call population
        endif
c
c Binary compound reactions
c
c compnorm  : subroutine for normalization of compound nucleus
c             cross section
c comptarget: subroutine for compound reaction for initial compound
c             nucleus
c
        if (flagcomp) then
          call compnorm
          call comptarget
        endif
c
c Collecting binary reaction results
c
c binary : subroutine for binary reaction results
c flagang: flag for output of angular distributions
c flagddx: flag for output of double-differential cross sections
c angdis : subroutine for calculation of angular distributions for
c          discrete states
c
        call binary
        if (flagang.or.flagddx.or.flagrecoil) call angdis
c
c Multiple emission
c
c multiple: subroutine for multiple emission
c
        call multiple
c
c Exclusive channels
c
c flagchannels: flag for exclusive channels calculation
c channels    : subroutine for exclusive reaction channels
c
        if (flagchannels) call channels
c
c Collecting total cross sections, spectra, angular distributions, etc.
c
c totalxs      : subroutine for total cross sections
c flagspec     : flag for output of spectra
c spectra      : subroutine for creation of particle spectra
c flagfission  : flag for fission
c flagmassdis  : flag for calculation of fission fragment mass yields
c massdis      : subroutine for fission fragment yields
c massdisout   : subroutine for output of fission fragment yields
c residual     : subroutine for residual production cross sections
c totalrecoil  : subroutine for recoil results
c flagrescue   : flag for final rescue: normalization to data
c normalization: subroutine to normalize cross sections to experimental
c                or evaluated data
c numinclow    : number of incident energies below Elow
c thermal      : subroutine for estimate of thermal cross sections
c flagurr      : flag for output of unresolved resonance parameters
c urr          : subroutine for unresolved resonance range parameters
c output       : subroutine for output
c
        call totalxs
        if (flagspec.or.flagddx) call spectra
        if (flagfission.and.flagmassdis) then
          call massdis
          call massdisout
        endif
        call residual
        if (flagrecoil) call totalrecoil
        if (flagrescue) call normalization
        if (nin.eq.numinclow+1.and.numinclow.gt.0) call thermal
   20   if (flagurr) call urr
        if (.not.flagastro) call output
   10 continue
c
c Final output
c
c astro       : subroutine for astrophysical reaction rates
c finalout    : subroutine for output of final results
c flagintegral: flag for calculation of effective cross section using
c               integral spectrum
c integral    : subroutine to calculate effective cross section for
c               integral spectrum
c flagendf    : flag for information for ENDF-6 file
c endf        : subroutine for cross sections and information for
c               ENDF-6 file
c flagprod    : flag for isotope production
c isoprod     : subroutine for isotope production
c flagmain    : flag for main output
c timer       : subroutine for output of execution time
c
      if (flagastro) then
        call astro
      else
        call finalout
      endif
      if (flagintegral) call integral
      if (flagendf) call endf
  100 if (flagprod) call isoprod
      if (flagmain) call timer
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
