      subroutine prodrates
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 28, 2013
c | Task  : Calculate reaction rates
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer   numint
      parameter (numint=100)
      integer   Ninte,nE,Zix,Nix,N,is,nen
      real      dE,dxdEsum,Eint(numint),dEdx,dxdE(numint),ratesum,
     +          Erpgrid(0:numen),E,Ea,Eb,xsa,xsb,xs
c
c *********** Determine integration grid and stopping power ************
c
c Ninte        : number of integration points
c dE           : energy step
c Ebeam        : incident energy in MeV for isotope production
c Eback        : lower end of energy range in MeV for isotope production
c dxdEsum      : integration sum
c Eint         : energy on integration grid
c stoppingpower: subroutine to calculate stopping power
c dEdx         : stopping power
c dxdE         : 1/stopping power
c targetdx     : effective thickness of target
c Vtar         : active target volume
c Mtar         : active target mass
c rhotarget    : target density
c Area         : target area in cm^2
c heat         : produced heat
c Ibeam        : beam current in mA
c
      Ninte=100
      dE=(Ebeam-Eback)/Ninte
      dxdEsum=0.
      do 10 nE=1,Ninte
        Eint(nE)=Eback+(nE-0.5)*dE
        call stoppingpower(Eint(nE),dEdx)
        dxdE(nE)=1./dEdx
        dxdEsum=dxdEsum+dxdE(nE)
   10 continue
      targetdx=dxdEsum*dE
      Vtar=Area*targetdx
      Mtar=rhotarget*Vtar
      heat=Ibeam*(Ebeam-Eback)
c
c ********************* Calculate reaction rates ***********************
c
c projnum  : number of incident particles [s^-1]
c parZ     : charge number of particle
c k0       : index of incident particle
c qelem    : elementary charge in C
c maxZ     : maximal number of protons away from the initial
c            compound nucleus
c maxN     : maximal number of neutrons away from the initial
c            compound nucleus
c Nisomer  : number of isomers for this nuclide
c prate    : production rate per isotope
c prodexist: flag for existence of residual product
c ratesum  : integration sum
c Nenrp,N  : number of incident energies for residual production cross
c            sections
c Erpgrid  : incident energy for cross section in MeV
c locate   : subroutine to find value in ordered table
c xsrp     : residual production cross section in mb
c pol1     : subroutine for interpolation of first order
c
      projnum=Ibeam/(1000.*parZ(k0)*qelem)
      do 110 Zix=-1,maxZ
        do 120 Nix=-1,maxN
          do 130 is=-1,Nisomer(Zix,Nix)
            prate(Zix,Nix,is)=0.
            ratesum=0.
            if (.not.prodexist(Zix,Nix,is)) goto 130
            N=Nenrp(Zix,Nix,is)
            if (N.eq.0) goto 130
            do 140 nen=1,N
              Erpgrid(nen)=Erp(Zix,Nix,is,nen)
  140       continue
            do 150 nE=1,Ninte
              E=Eint(nE)
              if (E.lt.Erpgrid(1)) goto 150
              if (E.gt.Erpgrid(N)) goto 160
              call locate(Erpgrid,1,N,E,nen)
              if (nen.eq.0) goto 150
              Ea=Erpgrid(nen)
              Eb=Erpgrid(nen+1)
              xsa=xsrp(Zix,Nix,is,nen)
              xsb=xsrp(Zix,Nix,is,nen+1)
              call pol1(Ea,Eb,xsa,xsb,E,xs)
              ratesum=ratesum+dxdE(nE)*xs
  150       continue
  160       prate(Zix,Nix,is)=projnum/Vtar*ratesum*dE*1.e-27
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C) 2010  A.J. Koning
