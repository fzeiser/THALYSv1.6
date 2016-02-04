      subroutine astroinit
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Stephane Goriely
c | Date  : December 13, 2013
c | Task  : Initialization of astrophysics quantities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer in,iz,nen,i
c
c **** Initialization of arrays of astrophysics interest ***************
c
c numN        : maximal number of neutrons away from the initial
c               compound nucleus
c numZ        : maximal number of protons away from the initial
c               compound nucleus
c numT        : number of temperatures
c numinc      : number of incident energies
c xsastro     : cross section for astrophysical calculation
c rateastro   : thermonuclear reaction rate factor
c macsastro   : Maxwellian-averaged thermonuclear reaction cross section
c partf       : integrated partition function
c rateastrofis: thermonuclear reaction rate factor for fission
c macsastro   : thermonuclear reaction cross section
c macsastrofis: thermonuclear reaction cross section for fission
c xsastrofis  : astrophysical fission cross section
c
      do 10 in=0,numN
        do 20 iz=0,numZ
          do 30 nen=1,numinc
            xsastro(iz,in,nen)=0.
   30     continue
          do 40 i=1,numT
            rateastro(iz,in,i)=0.
            macsastro(iz,in,i)=0.
   40     continue
   20   continue
   10 continue
      do 50 i=1,numT
        partf(i)=0.
        rateastrofis(i)=0.
        rateastroracap(i)=0.
        macsastrofis(i)=0.
        macsastroracap(i)=0.
   50 continue
      do 60 nen=1,numenin
        xsastrofis(nen)=0.
   60 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
