      subroutine msdcalc
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : MSD calculation
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ****************************** MSD model *****************************
c
c onestepB     : subroutine for continuum one-step direct cross sections
c flagonestep  : flag for continuum one-step direct only
c onecontinuumB: subroutine for one-step direct cross sections for MSD
c multistepA   : subroutine for multi-step direct cross sections
c multistepB   : subroutine for multi-step direct cross sections on
c                outgoing energy grid
c
      call onestepB
      if (.not.flagonestep) then
        call onecontinuumB
        call multistepA
        call multistepB
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
