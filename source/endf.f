      subroutine endf
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 19, 2013
c | Task  : Cross sections and information for ENDF-6 file
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ************************** ECIS calculation **************************
c
c numinc      : number of incident energies
c endfinfo    : subroutine for info for ENDF-6 file
c endfenergies: subroutine for energy grid for ENDF-6 file
c k0          : index for incident particle
c endfecis    : subroutine for ECIS calculation for incident particle on
c               ENDF-6 energy grid
c flagecisinp : flag for existence of ecis input file
c endfread    : subroutine to read ECIS results for incident particle on
c               ENDF-6 energy grid
c
      if (numinc.eq.1) return
      call endfinfo
      call endfenergies
      if (k0.ne.0) then
        call endfecis
        if (flagecisinp) call endfread
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
