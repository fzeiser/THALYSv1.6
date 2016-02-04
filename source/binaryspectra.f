      subroutine binaryspectra
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 28, 2004
c | Task  : Creation of binary spectra
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,iang
c
c *************** Interpolate decay on emission spectrum ***************
c
c binemission: subroutine for compound emission cross sections for
c              binary reaction
c parskip    : logical to skip outgoing particle
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c xscomp     : compound emission spectrum
c xsemis     : cross section for emission from compound nucleus
c xsbinemis  : cross section for emission from first compound nucleus
c xspreeq    : preequilibrium cross section per particle type and
c              outgoing energy
c xsgr       : smoothed giant resonance cross section
c flagddx    : flag for output of double-differential cross sections
c flagrecoil : flag for calculation of recoils
c nanglecont : number of angles for continuum
c xscompad   : compound emission angular distribution
c fourpi     : 4.*pi
c xsbinemisad: angular distribution for emission from first compound
c              nucleus
c xspreeqad  : preequilibrium angular distribution per particle type
c xsgrad     : smoothed giant resonance angular distribution
c
      call binemission
      do 10 type=0,6
        if (parskip(type)) goto 10
        do 20 nen=ebegin(type),eend(type)
          xscomp(type,nen)=xsemis(type,nen)
          xsbinemis(type,nen)=xscomp(type,nen)+xspreeq(type,nen)+
     +      xsgr(type,nen)
          if (flagddx.or.flagrecoil) then
            do 30 iang=0,nanglecont
              xscompad(type,nen,iang)=xscomp(type,nen)/fourpi
              xsbinemisad(type,nen,iang)=xscompad(type,nen,iang)+
     +          xspreeqad(type,nen,iang)+xsgrad(type,nen,iang)
   30       continue
          endif
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
