      function spindis(Zix,Nix,Eex,ald,Rspin,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 17, 2004
c | Task  : Wigner spin distribution
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      integer Zix,Nix,ibar
      real    Eex,ald,Rspin,sigma22,spincut,spindis
c
c *********************** Wigner formula ******************************
c
c spindis: Wigner spin distribution
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c Eex    : excitation energy
c ald    : level density parameter
c Rspin  : spin
c ibar   : fission barrier
c sigma22: 2 * spin cutoff factor
c spincut: spin cutoff factor
c
      sigma22=2.*spincut(Zix,Nix,ald,Eex,ibar)
      spindis=(2.*Rspin+1.)/sigma22*exp(-(Rspin+0.5)**2/sigma22)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
