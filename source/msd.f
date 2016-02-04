      subroutine msd
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 16, 2004
c | Task  : Multi-step direct model
c +---------------------------------------------------------------------
c
c ****************************** MSD model *****************************
c
c msdinit : subroutine for initialization of MSD model parameters
c dwbaecis: subroutine for ECIS calculations of DWBA for MSD
c msdcalc : subroutine for MSD calculation
c msdtotal: subroutine for total multi-step direct cross sections
c msdout  : subroutine for output of multi-step direct cross sections
c
c The MSD model is implemented for neutrons and protons only.
c
      call msdinit
      call dwbaecis
      call msdcalc
      call msdtotal
      call msdout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
