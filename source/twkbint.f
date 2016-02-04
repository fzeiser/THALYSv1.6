      function twkbint(efis,ibar,Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : December 18, 2009
c | Task  : Interpolation of WKB penetrability
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ibar,Zix,Nix,nen
      real    twkbint,efis,Ea,Eb,Ta,Tb,Tf,Ewkb(0:nbinswkb)
c
c ************************* Interpolation ******************************
c
      do nen=0,nbinswkb
        Ewkb(nen)=Uwkb(Zix,Nix,nen)
      enddo
      if (efis.gt.Ewkb(nbinswkb)) then
        twkbint=1.
      else
        call locate(Ewkb,0,nbinswkb,efis,nen)
        Ea=Ewkb(nen)
        Eb=Ewkb(nen+1)
        Ta=Twkb(Zix,Nix,nen,ibar)
        Tb=Twkb(Zix,Nix,nen+1,ibar)
        call pol1(Ea,Eb,Ta,Tb,efis,Tf)
        twkbint=Tf
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
