      function cmsd(Ein)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 4, 2012
c | Task  : Normalization factor for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      real         cmsd,Ein,factor,M2c
c
c ***************** Normalization factor for MSD ***********************
c
c cmsd       : normalization factor for MSD
c Ein        : incident energy
c preeqadjust: logical for energy-dependent pre-eq adjustment
c adjust     : subroutine for energy-dependent parameter adjustment
c factor     : multiplication factor
c M2constant : constant for matrix element in exciton model (here used
c              for MSD model)
c Atarget    : mass number of target nucleus
c
      if (preeqadjust) then
        key='m2constant'
        call adjust(Ein,key,0,0,0,0,factor)
        M2c=factor*M2constant
      else
        M2c=M2constant
      endif
      cmsd=M2c*3.3e6/(Atarget**3)/Ein
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
