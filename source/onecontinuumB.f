      subroutine onecontinuumB
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : One-step direct cross sections for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer itype,type,Zix,Nix,nen1,nen2,iang
      real    cmsd
c
c ************* Calculate continuum one-step cross sections ************
c
c parskip   : logical to skip outgoing particle
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c msdbins2  : number of energy points for MSD calculation
c Emsdin    : incident MSD energy
c specmass  : specific mass for target nucleus and all particles
c Emsd      : MSD energy grid
c xscont    : continuum one-step direct cross section for MSD
c cmsd      : normalization factor for MSD
c xscont1   : continuum one-step direct cross section for MSD
c             (unnormalized)
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c xscontad  : continuum one-step direct angular distribution for MSD
c xscontad1 : continuum one-step direct angular distribution for MSD
c             (unnormalized)
c
      do 10 itype=1,2
        if (parskip(itype)) goto 10
        do 20 type=1,2
          if (parskip(type)) goto 20
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          do 30 nen1=0,msdbins2
            Emsdin=real(specmass(Zix,Nix,itype)*Emsd(nen1))
            do 40 nen2=nen1,msdbins2
              xscont(itype,type,nen1,nen2)=cmsd(Emsdin)*
     +          xscont1(itype,type,nen1,nen2)
              if (flagddx) then
                do 50 iang=0,nanglecont
                  xscontad(itype,type,nen1,nen2,iang)=cmsd(Emsdin)*
     +              xscontad1(itype,type,nen1,nen2,iang)
   50           continue
              endif
   40       continue
   30     continue
   20   continue
   10 continue
c
c ******************* First-step cross sections for MSD ****************
c
c msdstep0  : n-step cross section for MSD
c k0        : index of incident particle
c msdstepad0: n-step angular distribution for MSD
c
      do 110 type=1,2
        if (parskip(type)) goto 110
        do 120 nen2=0,msdbins2
          msdstep0(type,1,nen2)=xscont(k0,type,0,nen2)
          if (flagddx) then
            do 140 iang=0,nanglecont
              msdstepad0(type,1,nen2,iang)=
     +          xscontad(k0,type,0,nen2,iang)
  140       continue
          endif
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
