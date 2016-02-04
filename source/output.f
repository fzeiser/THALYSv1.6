      subroutine output
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 19, 2012
c | Task  : Output
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ******************************* Output *******************************
c
c flagmain     : flag for main output
c flaginitpop  : flag for initial population distribution
c totalout     : subroutine for output of total cross sections
c binaryout    : subroutine for output of binary cross sections
c productionout: subroutine for output of particle production cross
c                sections
c residualout  : subroutine for output of residual production cross
c                sections
c flagfission  : flag for fission
c fissionout   : subroutine for output of fission cross sections
c flagdisc     : flag for output of discrete state cross sections
c discreteout  : subroutine for output of cross sections for discrete
c                states
c flagchannels : flag for exclusive channels calculation
c channelsout  : subroutine for output of exclusive reaction channels
c flagspec     : flag for output of spectra
c spectraout   : subroutine for output of particle spectra
c flagrecoil   : flag for calculation of recoils
c recoilout    : subroutine for output of recoils
c flagang      : flag for output of angular distributions
c angleout     : subroutine for output of discrete angular distributions
c flagddx      : flag for output of double-differential cross sections
c ddxout       : subroutine for output of double-differential cross
c                sections
c flaggamdis   : flag for output of discrete gamma-ray intensities
c gamdisout    : subroutine for output of discrete gamma-ray intensities
c flagracap    : flag for radiative capture model
c racapout     : subroutine for output of radiative capture model
c
      if (flagmain) then
        if (.not.flaginitpop) then
          call totalout
          call binaryout
        endif
        call productionout
        call residualout
        if (flagfission) call fissionout
      endif
      if (flagdisc) call discreteout
      if (flagchannels) call channelsout
      if (flagspec) call spectraout
      if (flagrecoil) call recoilout
      if (flagang) call angleout
      if (flagddx) call ddxout
      if (flaggamdis) call gamdisout
      if (flagracap) call racapout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
