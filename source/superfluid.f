      double precision function superfluid(Zix,Nix,ald,Eex,P,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : June 22, 2007
c | Task  : Superfluid model level density formula
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer          Zix,Nix,ibar
      real             ald,Eex,P,U,phi2,Df,phi1,Tf,Sf,sigma,spincut
      double precision fermi
c
c *********************** Level density formula ************************
c
c superfluid: level density
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c ald       : level density parameter
c Eex       : excitation energy
c P         : pairing energy
c ibar      : fission barrier number, zero for states on ground state
c U         : effective excitation energy
c pair      : total pairing correction
c Pshift    : adjustable pairing shift
c Ucrit     : critical U
c fermi     : function for Fermi gas level density formula
c phi2      : help variable
c Df        : determinant
c Dcrit     : critical determinant
c phi1      : phi function of superfluid model
c Tf        : temperature
c Tcrit     : critical temperature
c Sf        : entropy
c Scrit     : critical entropy
c sigma     : spin cutoff factor
c spincut   : spin cutoff factor
c
c Superfluid model
c
      U=Eex+pair(Zix,Nix)+Pshift(Zix,Nix,ibar)
      if (U.gt.0.) then
        if (U.gt.Ucrit(Zix,Nix)) then
          superfluid=fermi(Zix,Nix,ald,Eex,P,ibar)
        else
          phi2=1.-U/Ucrit(Zix,Nix)
          Df=Dcrit(Zix,Nix)*(1.-phi2)*(1.+phi2)*(1.+phi2)
          phi1=sqrt(phi2)
          Tf=2.*Tcrit(Zix,Nix)*phi1/log((phi1+1.)/(1.-phi1))
          Sf=Scrit(Zix,Nix)*Tcrit(Zix,Nix)/Tf*(1.-phi2)
          sigma=sqrt(spincut(Zix,Nix,ald,Eex,ibar))
          superfluid=exp(dble(Sf))/sqrt(Df)/sqrttwopi/sigma
        endif
      else
        superfluid=1.
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
