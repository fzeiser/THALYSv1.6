      function aldmatch(Zix,Nix,Eex,ibar)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire, Marieke Duijvestijn and Arjan Koning
c | Date  : April 7, 2006
c | Task  : Effective level density parameter
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,ibar,A,jmax,j
      real             aldmatch,Eex,rjbegin,rj,aldref,ignatyuk,spindis,
     +                 sigma,spincut,Krot,Kvib,Kcoll,aldacc,ald1,ald2,
     +                 dald,aldmid,sigmamid
      double precision rhosum,factor,fermi,rhoref,fdiff,fmid
      parameter (jmax=40)
c
c ************ Search for effective level density parameter ************
c
c The effective level density parameter is obtained in three steps:
c
c 1. Create total level density in Fermi gas region
c
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c AA,A         : mass number of residual nucleus
c rjbegin,rj   : help variables
c rhosum,factor: help variables
c aldref       : level density parameter
c ignatyuk     : function for energy dependent level density parameter a
c fermi        : function for Fermi gas level density formula
c pair         : pairing energy
c spindis      : Wigner spin distribution
c
      A=AA(Zix,Nix,0)
      rjbegin=0.5*mod(A,2)
      rj=rjbegin-1.
      rhosum=0.
      aldref=ignatyuk(Zix,Nix,Eex,ibar)
   10 rj=rj+1.
      factor=(2.*rj+1.)*fermi(Zix,Nix,aldref,Eex,pair(Zix,Nix),ibar)*
     +  spindis(Zix,Nix,Eex,aldref,rj,ibar)
      rhosum=rhosum+factor
      if (factor.gt.0.00001) goto 10
c
c 2. Apply a rotational enhancement to the total level density
c
c colenhance: subroutine for collective enhancement
c Krot      : rotational enhancement factor
c Kvib      : vibrational enhancement factor
c Kcoll     : total collective enhancement
c rhoref    : help variable
c
      call colenhance(Zix,Nix,Eex,aldref,ibar,Krot,Kvib,Kcoll)
      rhoref=Kcoll*rhosum
c
c 3. Determine effective level density parameter by equating the
c    rotational enhanced level density by a new effective
c    total level density.
c
c aldacc      : accuracy
c ald1,ald2   : boundary values for search
c ald         : level density parameter
c spincut     : spin cutoff factor
c sigma,rhonew: help variables
c rfit,isgn   : help variables
c sqrttwopi   : sqrt(2.*pi)
c aldmatch    : function to determine effective level density parameter
c
      aldacc=0.001
      ald1=0.5*aldref
      ald2=2.0*aldref
      sigma=sqrt(spincut(Zix,Nix,ald1,Eex,ibar))
      fdiff=fermi(Zix,Nix,ald1,Eex,pair(Zix,Nix),ibar)*sigma*
     +  sqrttwopi-rhoref
      if (fdiff.lt.0.) then
        aldmatch=ald1
        dald=ald2-ald1
      else
        aldmatch=ald2
        dald=ald1-ald2
      endif
      do 20 j=1,jmax
        dald=dald*0.5
        aldmid=aldmatch+dald
        sigmamid=sqrt(spincut(Zix,Nix,aldmid,Eex,ibar))
        fmid=fermi(Zix,Nix,aldmid,Eex,pair(Zix,Nix),ibar)*sigmamid*
     +    sqrttwopi-rhoref
        if (fmid.le.0.) aldmatch=aldmid
        if (abs(dald).lt.aldacc.or.fmid.eq.0.) return
   20 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
