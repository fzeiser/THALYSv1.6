      subroutine direct
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 14, 2008
c | Task  : Calculation of direct inelastic cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ****************** Direct cross section calculation ******************
c
c k0         : index for incident particle
c Ltarget    : excited level of target
c flagomponly: flag to execute ONLY an optical model calculation
c directecis : subroutine for ECIS calculation of direct cross section
c directread : subroutine to read ECIS results for direct cross section
c flaggiant  : flag for collective contribution from giant resonances
c giant      : subroutine for giant resonance contribution
c flagdirect : flag for output of direct reaction results
c directout  : subroutine for output of direct reaction cross sections
c
      if (k0.eq.0) return
      if (Ltarget.ne.0) return
      if (flagomponly) return
      call directecis
      call directread
      if (flaggiant) call giant
      if (flagdirect) call directout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
