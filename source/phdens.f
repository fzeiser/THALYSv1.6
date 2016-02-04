      function phdens(Zix,Nix,p,h,gs,Eex,Ewell,surfwell)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Particle-hole state density
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell
      integer Zix,Nix,p,h,n1,nex2
      real    phdens,gs,Eex,Ewell,Ap,fac1,factor,finitewell,eb,ee,ldb,
     +        lde,lldb,llde,ldtab
c
c ***** State density of Betak and Dobes, Z. Phys. A279 (1976) 319. ****
c
c phdens    : particle-hole state density
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c p         : particle number
c h         : hole number
c gs        : single-particle level density parameter
c Eex       : excitation energy
c Ewell     : depth of potential well
c surfwell  : flag for surface effects in finite well
c phmodel   : particle-hole state density model
c phexist1  : flag for existence of particle-hole state density table
c Apauli,Ap : Pauli blocking correction energy
c n1        : exciton number-1
c fac1      : help variable
c factor    : help variable
c nfac      : n!
c finitewell: function for correction for finite well depth
c
c In general, the finite depth of the hole is included. If the
c uncorrected state density is required, it should be specified by
c Ewell=0. (which thus actually means Ewell=infinity).
c
      phdens=0.
      if (p.lt.0.or.h.lt.0) return
      if (p+h.eq.0) return
      if (phmodel.eq.1.or..not.phexist1(Zix,Nix,p,h)) then
        Ap=Apauli(p,h)
        factor=(p*p+h*h+p+h)/(4.*gs)
        if (Ap+factor.ge.Eex) return
        n1=p+h-1
        fac1=nfac(p)*nfac(h)*nfac(n1)
        factor=gs**(p+h)/fac1
        phdens=factor*(Eex-Ap)**n1
        phdens=phdens*finitewell(p,h,Eex,Ewell,surfwell)
      else

c 2. Tabulated particle-hole state densities
c
c Ephdensmax   : maximum energy on particle-hole state density table
c nenphdens    : number of energies for particle-hole state density grid
c locate       : subroutine to find value in ordered table
c edens        : energy grid for tabulated level densities
c phtable1     : particle-hole state density from table
c eb,ee,ldb,...: help variables
c
        if (Eex.le.0.) return
        if (Eex.le.Ephdensmax) then
          call locate(edens,0,nenphdens,Eex,nex2)
          eb=edens(nex2)
          ee=edens(nex2+1)
          ldb=phtable1(Zix,Nix,p,h,nex2)
          lde=phtable1(Zix,Nix,p,h,nex2+1)
        else
          eb=edens(nenphdens-1)
          ee=edens(nenphdens)
          ldb=phtable1(Zix,Nix,p,h,nenphdens-1)
          lde=phtable1(Zix,Nix,p,h,nenphdens)
        endif
        if (ldb.gt.1..and.lde.gt.1.) then
          lldb=log(ldb)
          llde=log(lde)
          ldtab=exp(lldb+(Eex-eb)/(ee-eb)*(llde-lldb))
        else
          ldtab=ldb+(Eex-eb)/(ee-eb)*(lde-ldb)
        endif
        phdens=ldtab
      endif
      if (phdens.lt.1.e-10) phdens=0.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
