      subroutine incident
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 27, 2011
c | Task  : Main settings and basic cross sections for incident energy
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 yesno
c
c ************ Setting of some parameters for incident channel *********
c
c flagpreeq : flag for pre-equilibrium calculation
c flagmulpre: flag for multiple pre-equilibrium calculation
c flagmain  : flag for main output
c Einc      : incident energy in MeV
c yesno     : function to assign y or n to logical value
c flagwidth : flag for width fluctuation calculation
c flagurr   : flag for output of unresolved resonance parameters
c nbins     : number of continuum excitation energy bins
c
c Multiple pre-equilibrium emission is always set off if there is no
c primary pre-equilibrium emission.
c
      if (.not.flagpreeq) flagmulpre=.false.
c
c Write the energy dependent flags to the output file.
c
      if (flagmain) then
        if (Einc.ge.0.001) then
          write(*,'(/," ########## RESULTS FOR E=",f9.5,
     +      " ##########"/)') Einc
        else
          write(*,'(/," ########## RESULTS FOR E=",1p,e12.5,
     +      " ##########"/)') Einc
        endif
        write(*,'(" Energy dependent input flags"/)')
        write(*,'(" Width fluctuations (flagwidth)            : ",a1)')
     +    yesno(flagwidth)
        write(*,'(" Unresolved resonance parameters (flagurr) : ",a1)')
     +    yesno(flagurr)
        write(*,'(" Preequilibrium (flagpreeq)                : ",a1)')
     +    yesno(flagpreeq)
        write(*,'(" Multiple preequilibrium (flagmulpre)      : ",a1)')
     +    yesno(flagmulpre)
        write(*,'(" Number of continuum excitation energy bins:",i3)')
     +    nbins
      endif
c
c *** Calculation of total, reaction, elastic cross section,
c     transmission coefficients and elastic angular distribution for
c     incident energy ***
c
c k0           : index for incident particle
c coullimit    : energy limit for charged particle OMP calculation
c incidentecis : subroutine for ECIS calculation for incident energy
c incidentread : subroutine to read ECIS results for incident energy
c incidentnorm : subroutine for normalization of reaction cross sections
c                and transmission coefficients for incident channel
c flaginitpop  : flag for initial population distribution
c incidentgamma: subroutine for incident photons
c strength     : model for E1 gamma-ray strength function
c parinclude   : logical to include outgoing particle
c flagcomp     : flag for compound nucleus calculation
c radwidtheory : subroutine for theoretical calculation of total
c                radiative width
c gammanorm    : subroutine for normalization of gamma ray strength
c                functions
c spr          : subroutine for S, P and R' resonance parameters
c flaginverse  : flag for output of transmission coefficients and
c                inverse reaction cross sections
c incidentout  : subroutine for reaction output for incident channel
c
      if (k0.gt.1.and.Einc.lt.coullimit(k0)) return
      if (k0.gt.0) then
        call incidentecis
        call incidentread
        call incidentnorm
      else
        if (.not.flaginitpop) call incidentgamma
      endif
      if ((parinclude(0).or.flagcomp).and.Einc.le.100.) then
        call radwidtheory(0,0,Einc)
        if (strength.eq.1) call gammanorm(0,0)
      endif
      if (k0.eq.1.and.(parinclude(0).or.flagcomp)) call spr
      if (flaginverse) call incidentout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
