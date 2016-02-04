      subroutine basicxs(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 23, 2009
c | Task  : Basic cross sections and transmission coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp
c
c ******************* ECIS calculations and output *********************
c
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c flagomponly : flag to execute ONLY an optical model calculation
c flagcomp    : flag for compound nucleus calculation
c basicinitial: subroutine for initialization of arrays for basic cross
c               sections
c inverse     : subroutine for ECIS calculation of total, reaction and
c               elastic cross sections and transmission coefficients
c               for outgoing energy grid.
c
c The transmission coefficients and inverse reaction cross sections for
c the outgoing energy grid need to be calculated only once, for the
c maximal incident energy.
c
      if (flagomponly.and..not.flagcomp) return
      call basicinitial
      call inverse(Zcomp,Ncomp)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
