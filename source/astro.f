      subroutine astro
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely and Stephane Hilaire
c | Date  : February 8, 2006
c | Task  : Astrophysical reaction rates
c +---------------------------------------------------------------------
c
c ************** Astrophysical reaction rate calculation ***************
c
c stellarrate: subroutine for calculation of reaction rate for a
c              Maxwell-Boltzmann distribution
c astroout   : subroutine for output of astrophysical results
c
      call stellarrate
      call astroout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
