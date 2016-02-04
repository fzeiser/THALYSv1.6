      double precision function fermi(Zix,Nix,ald,Eex,P,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 3, 2013
c | Task  : Fermi gas level density formula
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      integer          Zix,Nix,ibar
      real             ald,Eex,P,U,sigma,spincut,denom
      double precision factor
c
c *********************** Total level density **********************
c
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c ald         : level density parameter
c Eex         : excitation energy
c P           : pairing energy
c ibar        : fission barrier
c U           : excitation energy minus pairing energy
c factor,denom: help variable
c sigma       : square root of spin cutoff factor
c spincut     : spin cutoff factor
c AA          : mass number of nucleus
c
      U=Eex-P
      if (U.gt.0.) then
        factor=min(2.*sqrt(ald*U),700.)
        sigma=sqrt(spincut(Zix,Nix,ald,Eex,ibar))
        denom=12.*sqrt(2.)*sigma*(ald**0.25)*(U**1.25)
        fermi=exp(factor)/denom
      else
        fermi=1.
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
