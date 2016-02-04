      function ignatyuk(Zix,Nix,Eex,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 22, 2011
c | Task  : Energy dependent level density parameter a
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer Zix,Nix,ibar,A
      real    ignatyuk,Eex,U,aldlim,damp,fU,aldlow,expo,qfermi
c
c *********************** Level density formula ************************
c
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c ignatyuk: function for energy dependent level density parameter a
c Eex     : excitation energy
c ibar    : fission barrier
c U       : excitation energy minus pairing energy
c delta   : energy shift
c alimit  : asymptotic level density parameter
c aldlim  : asymptotic level density parameter
c fU,damp : help variables
c gammald : gamma-constant for asymptotic level density parameter
c deltaW  : shell correction in nuclear mass
c
c Formalism from Ignatyuk et al. Sov. Jour. Nuc. Phys. 21 (1975), 255.
c
      U=Eex-delta(Zix,Nix,ibar)
      aldlim=alimit(Zix,Nix)
c
c 1. For very low Eex, i.e. U < 0, we use the first order Taylor
c    expansion
c
      if (U.le.0.) then
        damp=(1.+deltaW(Zix,Nix,ibar)*gammald(Zix,Nix))
      else
c
c 2. Higher energies
c
        expo=gammald(Zix,Nix)*U
        fU=1.
        if (abs(expo).le.80.) fU=1.-exp(-expo)
        damp=1.+fU*deltaW(Zix,Nix,ibar)/U
      endif
c
c Only for fission models with damping of collective effects in
c effective level density model. (Re-installed from TALYS-0.64).
c Fermi distribution for asymptotic level density parameter. This takes
c the damping of collective effects into account in a phenomenological
c way. The level density parameters are then equal to A/13 in the
c high-energy limit rather than A/8.
c
c flagcolldamp: flag for damping of collective effects in effective
c               level density (without explicit collective enhancement)
c               Only used for Bruyeres-le-Chatel (Pascal Romain) fission
c               model
c AA,A        : mass number of residual nucleus
c aldlow      : lower limit of a
c expo        : help variable
c Ufermi      : energy of Fermi distribution for damping of ground-state
c cfermi      : width of Fermi distribution for damping of ground-state
c qfermi      : Fermi distribution
c
      if (flagcolldamp) then
        A=AA(Zix,Nix,0)
        aldlow=A/13.
        expo=(U-Ufermi)/cfermi
        qfermi=0.
        if (expo.gt.-80.) qfermi=1./(1.+exp(-expo))
        aldlim=aldlow*qfermi+aldlim*(1.-qfermi)
      endif
      ignatyuk=max(aldlim*damp,1.)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
