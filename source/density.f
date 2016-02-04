      function density(Zix,Nix,Eex,Rspin,parity,ibar,ldmod)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 21, 2013
c | Task  : Level density
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
c Note that the parity is not used explicitly for analytical level
c densities. For those, an equidistant parity distribution is assumed.
c
      include "talys.cmb"
      character*80     key
      integer          Zix,Nix,parity,ibar,ldmod,jj,nex2
      real             Eex,Rspin,ald,ignatyuk,spindis,factor,pt,ct,
     +                 Eshift,expo
      double precision density,densitytot,eb,ee,ldb,lde,ldtab,cctable
c
c ************************** Level density *****************************
c
c density   : level density
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c Eex       : excitation energy
c Rspin     : spin
c parity    : parity
c ibar      : fission barrier number, zero for states on ground state
c ldmod     : level density model
c ldexist   : flag for existence of level density table
c ald       : level density parameter
c ignatyuk  : function for energy dependent level density parameter a
c densitytot: total level density
c pardis    : parity distribution
c spindis   : Wigner spin distribution
c
c 1. Gilbert and Cameron
c 2. Back-shifted Fermi gas
c 3. Superfluid model
c
      density=0.
      if (Eex.lt.0.) goto 100
      if (ldmod.le.3.or..not.ldexist(Zix,Nix,ibar)) then
        ald=ignatyuk(Zix,Nix,Eex,ibar)
        density=densitytot(Zix,Nix,Eex,ibar,ldmod)*pardis*
     +    spindis(Zix,Nix,Eex,ald,Rspin,ibar)
      else
c
c 4. Tabulated level densities
c
c Eshift       : shifted excitation energy
c ctable,ptable: constant to adjust tabulated level densities
c adjust       : subroutine for energy-dependent parameter adjustment
c factor       : multiplication factor
c ldadjust     : logical for energy-dependent level density adjustment
c Edensmax     : maximum energy on level density table
c locate       : subroutine to find value in ordered table
c edens        : energy grid for tabulated level densities
c nendens      : number of energies for level density grid
c eb,ee,ldb,lde: help variables
c ldtable      : level density from table
c expo         : exponent
c cctable      : constant to adjust tabulated level densities
c
        jj=min(29,int(Rspin))
c
c Possible adjustment of final level densities
c
        if (ldadjust(Zix,Nix)) then
          key='ptable'
          call adjust(Eex,key,Zix,Nix,0,ibar,factor)
          pt=ptable(Zix,Nix,ibar)+factor-1.
          key='ctable'
          call adjust(Eex,key,Zix,Nix,0,ibar,factor)
          ct=ctable(Zix,Nix,ibar)+factor-1.
        else
          pt=ptable(Zix,Nix,ibar)
          ct=ctable(Zix,Nix,ibar)
        endif
        Eshift=Eex-pt
c
c Interpolate from tables
c
        if (Eshift.le.0.) goto 100
        if (Eshift.le.Edensmax(Zix,Nix)) then
          call locate(edens,0,nendens(Zix,Nix),Eshift,nex2)
          eb=edens(nex2)
          ee=edens(nex2+1)
          ldb=ldtable(Zix,Nix,nex2,jj,parity,ibar)
          lde=ldtable(Zix,Nix,nex2+1,jj,parity,ibar)
        else
          eb=edens(nendens(Zix,Nix)-1)
          ee=edens(nendens(Zix,Nix))
          ldb=ldtable(Zix,Nix,nendens(Zix,Nix)-1,jj,parity,ibar)
          lde=ldtable(Zix,Nix,nendens(Zix,Nix),jj,parity,ibar)
        endif
        if (ldb.le.1..or.lde.le.1.) then
          ldtab=ldb+(Eshift-eb)/(ee-eb)*(lde-ldb)
        else
          ldtab=exp(log(ldb)+(Eshift-eb)/(ee-eb)*
     +     (log(lde)-log(ldb)))
        endif
        expo=min(ct*sqrt(Eshift),80.)
        cctable=exp(dble(expo))
        density=cctable*ldtab
      endif
  100 density=max(density,1.d-30)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
