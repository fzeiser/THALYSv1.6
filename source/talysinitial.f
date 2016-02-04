      subroutine talysinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 23, 2012
c | Task  : Initialization of nuclear structure
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ********** Initialization of constants and nuclear structure *********
c
c particles   : subroutine to determine included light particles
c nuclides    : subroutine for properties of nuclides
c flagreaction: flag for calculation of nuclear reactions
c grid        : subroutine for energy and angle grid
c flagmain    : flag for main output
c mainout     : subroutine for main output
c
      call particles
      call nuclides
      if (flagreaction) call grid
      if (flagmain) call mainout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
