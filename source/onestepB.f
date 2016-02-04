      subroutine onestepB
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : Continuum one-step direct cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,iang
      real    cmsd
c
c ************* Calculate continuum one-step cross sections ************
c
c parskip   : logical to skip outgoing particle
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c msdstep   : continuum n-step direct cross section
c cmsd      : normalization factor for MSD
c eninccm   : center-of-mass incident energy in MeV
c msdstep1  : continuum one-step direct cross section (unnormalized)
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c msdstepad : continuum n-step direct angular distribution
c msdstepad1: continuum one-step direct angular distribution
c             (unnormalized)
c
      do 10 type=1,2
        if (parskip(type)) goto 10
        do 20 nen=ebegin(type),eend(type)
          msdstep(type,1,nen)=cmsd(eninccm)*msdstep1(type,nen)
          if (flagddx) then
            do 30 iang=0,nanglecont
              msdstepad(type,1,nen,iang)=cmsd(eninccm)*
     +          msdstepad1(type,nen,iang)
   30       continue
          endif
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
