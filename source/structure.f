      subroutine structure(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 2, 2011
c | Task  : Nuclear structure parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix
c
c ********** Calculate and read various nuclear parameters *************
c
c Zix           : charge number index for residual nucleus
c Nix           : neutron number index for residual nucleus
c levels        : subroutine for discrete levels
c flagendf      : flag for information for ENDF-6 file
c primary       : flag to designate primary (binary) reaction
c gammadecay    : subroutine for scheme for discrete gamma decay
c deformpar     : subroutine for deformation parameters
c parinclude    : logical to include outgoing particle
c flagcomp      : flag for compound nucleus calculation
c resonancepar  : subroutine for s-wave resonance parameters
c gammapar      : subroutine for gamma ray parameters
c flagompall    : flag for new optical model calculation for all
c                 residual nuclides
c omppar        : subroutine for optical model parameters
c flagjlm       : flag for using semi-microscopic JLM OMP
c alphaomp      : alpha optical model (1=normal, 2= McFadden-Satchler,
c                 3-5= folding potential)
c radialtable   : subroutine for tabulated radial matter densities
c flagomponly   : flag to execute ONLY an optical model calculation
c flagfission   : flag for fission
c fissionpar    : subroutine for fission parameters
c densitypar    : subroutine for level density parameters
c ldmodel       : level density model
c densitytable  : subroutine for tabulated level densities
c densitymatch  : subroutine for level density matching solution
c phmodel       : particle-hole state density model
c phdensitytable: subroutine for tabulated particle-hole state densities
c k0            : index of incident particle
c parZ          : charge number of particle
c parN          : neutron number of particle
c thermalxs     : subroutine for cross sections at thermal energies
c flagpartable  : flag for output of model parameters on separate file
c partable      : subroutine to write model parameters per nucleus to
c                 separate file
c
c All the nuclear structure info is read and/or calculated for the
c nucleus under consideration.
c
      call levels(Zix,Nix)
      if (flagendf.and.primary) call gammadecay(Zix,Nix)
      call deformpar(Zix,Nix)
      if (parinclude(0).or.flagcomp) then
        call resonancepar(Zix,Nix)
        call gammapar(Zix,Nix)
      endif
      if ((Zix.le.2.and.Nix.le.2).or.flagompall) call omppar(Zix,Nix)
      if (flagjlm.or.alphaomp.ge.3) call radialtable(Zix,Nix)
      if (flagomponly.and..not.flagcomp) return
      if (flagfission) call fissionpar(Zix,Nix)
      call densitypar(Zix,Nix)
      if (ldmodel(Zix,Nix).ge.4) call densitytable(Zix,Nix)
      call densitymatch(Zix,Nix)
      if (phmodel.eq.2.and.Zix.le.numZph.and.Nix.le.numNph)
     +  call phdensitytable(Zix,Nix)
      if (k0.eq.1.and.Zix.eq.parZ(k0).and.Nix.eq.parN(k0))
     + call thermalxs
      if (flagpartable) call partable(Zix,Nix)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
