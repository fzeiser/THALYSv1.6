      subroutine matching(Zix,Nix,Exm,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 16, 2007
c | Task  : Determine matching between temperature and Fermi-gas region
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      external match
      integer  Zix,Nix,ibar,A,nseg,Nb
      real     Exm,ald,ignatyuk,xacc,x1,x2,match,rtbis,xb1(2),xb2(2),
     +         Exm1,Exm2
c
c ************************ Search for zeroes ***************************
c
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c Exm     : matching energy
c ibar    : fission barrier
c AA,A    : mass number of residual nucleus
c Exmemp  : empirical estimate for matching point for Ex
c Exmemp  : empirical values for matching energy
c pair    : total pairing correction
c ald     : level density parameter
c ignatyuk: function for energy dependent level density parameter a
c xacc,...: help variables
c Nb      : number of solutions
c match   : function to search for zero crossings of the function
c zbrak   : function to bracket the function
c
      A=AA(Zix,Nix,0)
      ald=ignatyuk(Zix,Nix,Exmemp,ibar)
      xacc=0.0001
c
c Set possible region for solution
c
      x1=max(2.25/ald+pair(Zix,Nix),0.)+0.11
      x2=19.+300./A
      nseg=100
      Nb=2
      call zbrak(match,x1,x2,nseg,xb1,xb2,Nb)
c
c Look for 0,1 or 2 solutions
c
      if (Nb.eq.0) Exm=0.
      if (Nb.eq.1) Exm=rtbis(match,xb1(1),xb2(1),xacc)
c
c If there are 2 solutions, we choose the one closest to the empirical
c expression.
c
      if (Nb.eq.2) then
        Exm1=rtbis(match,xb1(1),xb2(1),xacc)
        Exm2=rtbis(match,xb1(2),xb2(2),xacc)
        if (abs(Exm1-Exmemp).gt.abs(Exm2-Exmemp)) then
          Exm=Exm2
        else
          Exm=Exm1
        endif
      endif
c
c If the solution is unphysical, we choose the empirical expression.
c
      if (Exm.lt.x1.or.Exm.gt.x2) Exm=0.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
