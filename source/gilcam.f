      double precision function gilcam(Zix,Nix,ald,Eex,P,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 7, 2005
c | Task  : Gilbert-Cameron level density formula
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer          Zix,Nix,ibar
      real             ald,Eex,P
      double precision fermi,expo
c
c *********************** Level density formula ************************
c
c gilcam : level density
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c ald    : level density parameter
c Eex    : excitation energy
c P      : pairing energy
c ibar   : fission barrier number, zero for states on ground state
c Exmatch: matching point for Ex
c fermi  : function for Fermi gas level density formula
c pair   : pairing energy
c expo   : help variable
c E0     : constant of temperature formula
c T      : nuclear temperature
c
      if (Eex.gt.Exmatch(Zix,Nix,ibar)) then
c
c 1. Fermi Gas
c
        gilcam=fermi(Zix,Nix,ald,Eex,P,ibar)
      else
c
c 2. Constant temperature
c
        expo=min((Eex-E0(Zix,Nix,ibar))/T(Zix,Nix,ibar),300.)
        gilcam=exp(expo)/T(Zix,Nix,ibar)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
