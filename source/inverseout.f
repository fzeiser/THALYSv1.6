      subroutine inverseout(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 5, 2006
c | Task  : Reaction output for outgoing channels
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,type,Zix,Nix,l,nen
      real    e
c
c **************** Transmission coefficients per energy ****************
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c flagtransen: flag for output of transmission coefficients per energy
c parskip    : logical to skip outgoing particle
c Zindex     : charge number index for residual nucleus
c Nindex     : neutron number index for residual nucleus
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c e          : outgoing energy in CM
c egrid      : outgoing energy grid
c specmass   : specific mass
c parname    : name of particle
c lmax       : maximal l-value for transmission coefficients
c Tjl        : transmission coefficients as a function of particle
c              type, energy, spin and l-value
c Tl         : transmission coefficients as a function of particle
c              type, energy and l-value (averaged over spin)
c
      write(*,'(/" ########## TRANSMISSION COEFFICIENTS AND",
     +  " INVERSE REACTION CROSS SECTIONS ##########")')
c
c For each energy, the whole set of transmission coefficients is given
c as a function of the l-value and spin value.
c
      if (flagtransen) then
        do 10 type=1,6
          if (parskip(type)) goto 10
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          do 20 nen=ebegin(type),eend(type)
            e=real(egrid(nen)/specmass(Zix,Nix,type))
            write(*,'(/" Transmission coefficients for incident ",a8,
     +        " at ",f9.5," MeV"/)') parname(type),e
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
            if (type.ne.3.and.type.ne.6.and.lmax(type,nen).ge.0) then
              write(*,'("   L   T(L-1/2,L)   T(L+1/2,L)    Tav(L)"/)')
              do 30 l=0,lmax(type,nen)
                write(*,'(1x,i3,1p,3e13.5)') l,Tjl(type,nen,-1,l),
     +            Tjl(type,nen,1,l),Tl(type,nen,l)
   30         continue
            endif
c
c 2. Spin 1 particles: Deuterons
c
            if (type.eq.3.and.lmax(type,nen).ge.0) then
              write(*,'("   L    T(L-1,L)     T(L,L)       ",
     +          "T(L+1,L)     Tav(L)"/)')
              do 40 l=0,lmax(type,nen)
                write(*,'(1x,i3,1p,4e13.5)') l,Tjl(type,nen,-1,l),
     +            Tjl(type,nen,0,l),Tjl(type,nen,1,l),Tl(type,nen,l)
   40         continue
            endif
c
c 3. Spin 0 particles: Alpha-particles
c
            if (type.eq.6.and.lmax(type,nen).ge.0) then
              write(*,'("   L     T(L)"/)')
              do 50 l=0,lmax(type,nen)
                write(*,'(1x,i3,1p,e13.5)') l,Tjl(type,nen,0,l)
   50         continue
            endif
   20     continue
   10   continue
c
c ************ Transmission coefficients per angular momentum **********
c
      else
c
c For each l-value, the whole set of transmission coefficients is given
c as a function of the energy and spin value.
c
        do 110 type=1,6
          if (parskip(type)) goto 110
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          if (ebegin(type).ge.eend(type)) goto 110
          do 120 l=0,lmax(type,eend(type))
            write(*,'(/" Transmission coefficients for incident ",a8,
     +        " and l= ",i2/)') parname(type),l
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
            if (type.ne.3.and.type.ne.6) then
              write(*,'("   Energy    T(L-1/2,L)   T(L+1/2,L)     ",
     +          "Tav(L)"/)')
              do 130 nen=ebegin(type),eend(type)
                e=real(egrid(nen)/specmass(Zix,Nix,type))
                write(*,'(1x,f9.5,1p,3e13.5)') e,
     +            Tjl(type,nen,-1,l),Tjl(type,nen,1,l),Tl(type,nen,l)
  130         continue
            endif
c
c 2. Spin 1 particles: Deuterons
c
            if (type.eq.3) then
              write(*,'("   Energy     T(L-1,L)     T(L,L)       ",
     +          "T(L+1,L)     Tav(L)"/)')
              do 140 nen=ebegin(type),eend(type)
                e=real(egrid(nen)/specmass(Zix,Nix,type))
                write(*,'(1x,f9.5,1p,4e13.5)') e,
     +            Tjl(type,nen,-1,l),Tjl(type,nen,0,l),
     +            Tjl(type,nen,1,l),Tl(type,nen,l)
  140         continue
            endif
c
c 3. Spin 0 particles: Alpha-particles
c
            if (type.eq.6) then
              write(*,'("   Energy      T(L)"/)')
              do 150 nen=ebegin(type),eend(type)
                e=real(egrid(nen)/specmass(Zix,Nix,type))
                write(*,'(1x,f9.5,1p,e13.5)') e,Tjl(type,nen,0,l)
  150         continue
            endif
  120     continue
  110   continue
      endif
c
c **************** Cross sections for inverse channels *****************
c
c xstot : total cross section (neutrons only)
c xsreac: reaction cross section
c xselas: total elastic cross section (neutrons only)
c xsopt : optical model reaction cross section
c
      do 210 type=1,6
        if (parskip(type)) goto 210
        if (ebegin(type).ge.eend(type)) goto 210
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        if (type.eq.1) then
          write(*,'(/" Total cross sections for ",a8/)') parname(1)
          write(*,'("     E        total      reaction    elastic",
     +      "   OMP reaction"/)')
          do 220 nen=ebegin(1),eend(type)
            e=real(egrid(nen)/specmass(Zix,Nix,type))
            write(*,'(1x,f9.5,1p,4e12.4)') e,xstot(1,nen),
     +        xsreac(1,nen),xselas(1,nen),xsopt(1,nen)
  220     continue
        else
          write(*,'(/" Total cross sections for ",a8/)') parname(type)
          write(*,'("     E       reaction  OMP reaction"/)')
          do 230 nen=ebegin(type),eend(type)
            e=real(egrid(nen)/specmass(Zix,Nix,type))
            write(*,'(1x,f9.5,1p,2e12.4)') e,xsreac(type,nen),
     +        xsopt(type,nen)
  230     continue
        endif
  210 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
