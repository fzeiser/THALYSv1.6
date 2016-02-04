      subroutine separation
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 5, 2009
c | Task  : Separation energies
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,type,Zr,Nr
c
c Zix,Zr : charge number index for residual nucleus
c maxZ   : maximal number of protons away from initial compound nucleus
c Nix,Nr : neutron number index for residual nucleus
c maxN   : maximal number of neutrons away from initial compound nucleus
c parZ   : charge number of particle
c parN   : neutron number of particle
c expmexc: experimental mass excess
c S      : separation energy per particle
c excmass: mass excess of particle in a.m.u.
c amu    : atomic mass unit in MeV
c thmexc : theoretical mass excess
c dumexc : theoretical mass excess from Duflo-Zuker formula
c
c For consistency, separation energies are always calculated using two
c nuclear masses of the same type, i.e. both experimental or both
c theoretical. Hence if (Zix,Nix) is in the Audi-Wapstra table but
c (Zr,Nr) is not, for both nuclides the theoretical masses are used.
c For the calculation of separation energies, Zix and Nix act as
c compound nucleus indices.
c
      do 10 Zix=0,maxZ+2
        do 10 Nix=0,maxN+2
          do 20 type=1,6
            Zr=Zix+parZ(type)
            Nr=Nix+parN(type)
            if (expmexc(Zix,Nix).ne.0..and.expmexc(Zr,Nr).ne.0.) then
              S(Zix,Nix,type)=real(expmexc(Zr,Nr)-expmexc(Zix,Nix)+
     +          excmass(type)*amu)
              goto 20
            endif
            if (thmexc(Zix,Nix).ne.0..and.thmexc(Zr,Nr).ne.0.) then
              S(Zix,Nix,type)=real(thmexc(Zr,Nr)-thmexc(Zix,Nix)+
     +          excmass(type)*amu)
            else
              S(Zix,Nix,type)=real(dumexc(Zr,Nr)-dumexc(Zix,Nix)+
     +          excmass(type)*amu)
            endif
   20     continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
