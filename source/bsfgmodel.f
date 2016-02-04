      double precision function bsfgmodel(Zix,Nix,ald,Eex,P,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2006
c | Task  : Back-shifted Fermi gas level density formula
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      integer          Zix,Nix,ibar
      real             ald,Eex,P,sigma,spincut,U,an,ap,T2
      double precision invfermi,fermi,term,expo,deninv
c
c *********************** Level density formula ************************
c
c bsfgmodel: level density
c Zix      : charge number index for residual nucleus
c Nix      : neutron number index for residual nucleus
c ald      : level density parameter
c Eex      : excitation energy
c P        : pairing energy
c ibar     : fission barrier number, zero for states on ground state
c sigma    : square root of spin cutoff factor
c spincut  : spin cutoff factor
c an       : neutron level density parameter
c ap       : proton level density parameter
c U        : excitation energy minus pairing energy
c invfermi : help variable
c fermi    : function for Fermi gas level density formula
c T2       : square of temperature
c expo,term: help variables
c deninv   : help variable
c
c Back-shifted Fermi gas
c
c We apply the method of Grossjean and Feldmeier, as implemented by
c Demetriou and Goriely, to avoid the unphysical divergence near zero
c energy. The contribution given by the 1./term goes rapidly to zero
c with increasing excitation energy.
c
      sigma=sqrt(spincut(Zix,Nix,ald,Eex,ibar))
      an=0.5*ald
      ap=0.5*ald
      term=exp(1.)/24.*(an+ap)**2/sqrt(an*ap)/sigma
      U=Eex-P
      if (U.gt.0.) then
        invfermi=1./fermi(Zix,Nix,ald,Eex,P,ibar)
        T2=U/ald
        expo=4.*an*ap*T2
        if (expo.lt.80.) then
          deninv=invfermi+1./(term*exp(expo))
        else
          deninv=invfermi
        endif
        bsfgmodel=1./deninv
      else
        bsfgmodel=term
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
