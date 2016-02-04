      subroutine optical(Zix,Nix,kopt,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Determination of optical potential
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer  Zix,Nix,kopt
      real     eopt
c
c ********************** Call optical model module *********************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c kopt       : index for fast particle
c eopt       : incident energy
c opticalnp  : optical model for neutrons and protons
c opticalcomp: optical model for composite particles
c
      eopt=max(eopt,0.)
      if (kopt.le.2) then
        call opticalnp(Zix,Nix,kopt,eopt)
      else
        call opticalcomp(Zix,Nix,kopt,eopt)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
