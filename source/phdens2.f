      function phdens2(Zix,Nix,ppi,hpi,pnu,hnu,gsp,gsn,Eex,Ewell,
     +  surfwell)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Two-component particle-hole state density
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell
      integer Zix,Nix,ppi,hpi,pnu,hnu,h,p,npi,nnu,n1,nex2
      real    phdens2,gsp,gsn,Eex,Ewell,Ap,factorn,factorp,fac1,factor,
     +        finitewell,eb,ee,ldb,lde,lldb,llde,ldtab
c
c 1. State density of Kalbach, Phys. Rev. C33, 818 (1986).
c
c phdens2   : two-component particle-hole state density
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c ppi       : proton particle number
c hpi       : proton hole number
c pnu       : neutron particle number
c hnu       : neutron hole number
c gsp       : single-particle proton level density parameter
c gsn       : single-particle neutron level density parameter
c Eex       : excitation energy
c Ewell     : depth of potential well
c surfwell  : flag for surface effects in finite well
c phmodel   : particle-hole state density model
c phexist2  : flag for existence of particle-hole state density table
c Apauli2,Ap: two-component Pauli blocking correction factor
c factorn,..: help variables
c h         : hole number
c p         : particle number
c npi       : proton exciton number
c nnu       : neutron exciton number
c n1        : exciton number-1
c fac1      : help variable
c nfac      : n!
c factor    : help variable
c finitewell: function for correction for finite well depth
c
c The finite depth of the hole is included. If the uncorrected state
c density is required, it should be specified by Ewell=0. (which thus
c actually means Ewell=infinity).
c
      phdens2=0.
      if (ppi.lt.0.or.hpi.lt.0.or.pnu.lt.0.or.hnu.lt.0) return
      if (ppi+hpi+pnu+hnu.eq.0) return
      if (phmodel.eq.1.or..not.phexist2(Zix,Nix,ppi,hpi,pnu,hnu)) then
        Ap=Apauli2(ppi,hpi,pnu,hnu)
        factorn=(pnu*pnu+hnu*hnu+pnu+hnu)/(4.*gsn)
        factorp=(ppi*ppi+hpi*hpi+ppi+hpi)/(4.*gsp)
        if (Ap+factorn+factorp.ge.Eex) return
        h=hpi+hnu
        p=ppi+pnu
        npi=ppi+hpi
        nnu=pnu+hnu
        n1=npi+nnu-1
        fac1=nfac(ppi)*nfac(hpi)*nfac(pnu)*nfac(hnu)*nfac(n1)
        factor=gsp**npi*gsn**nnu/fac1
        phdens2=factor*(Eex-Ap)**n1
        phdens2=phdens2*finitewell(p,h,Eex,Ewell,surfwell)
      else
c
c 2. Tabulated particle-hole state densities
c
c Ephdensmax   : maximum energy on particle-hole state density table
c nenphdens    : number of energies for particle-hole state density grid
c locate       : subroutine to find value in ordered table
c edens        : energy grid for tabulated level densities
c phtable2     : particle-hole state density from table
c eb,ee,ldb,...: help variables
c
        if (Eex.le.0.) return
        if (Eex.le.Ephdensmax) then
          call locate(edens,0,nenphdens,Eex,nex2)
          eb=edens(nex2)
          ee=edens(nex2+1)
          ldb=phtable2(Zix,Nix,ppi,hpi,pnu,hnu,nex2)
          lde=phtable2(Zix,Nix,ppi,hpi,pnu,hnu,nex2+1)
        else
          eb=edens(nenphdens-1)
          ee=edens(nenphdens)
          ldb=phtable2(Zix,Nix,ppi,hpi,pnu,hnu,nenphdens-1)
          lde=phtable2(Zix,Nix,ppi,hpi,pnu,hnu,nenphdens)
        endif
        if (ldb.gt.1..and.lde.gt.1.) then
          lldb=log(ldb)
          llde=log(lde)
          ldtab=exp(lldb+(Eex-eb)/(ee-eb)*(llde-lldb))
        else
          ldtab=ldb+(Eex-eb)/(ee-eb)*(lde-ldb)
        endif
        phdens2=ldtab
      endif
      if (phdens2.lt.1.e-10) phdens2=0.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
