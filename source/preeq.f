      subroutine preeq
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 26, 2009
c | Task  : Pre-equilibrium reactions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real surface
c
c ********************** General initializations ***********************
c
c flagomponly: flag to execute ONLY an optical model calculation
c p0         : initial particle number
c parA       : mass number of particle
c k0         : index of incident particle
c h0         : initial hole number
c ppi0       : initial proton number
c parZ       : charge number of particle
c hpi0       : initial proton hole number
c pnu0       : initial neutron number
c parN       : neutron number of particle
c hnu0       : initial neutron hole number
c
      if (flagomponly) return
      p0=parA(k0)
      h0=0
      ppi0=parZ(k0)
      hpi0=0
      pnu0=parN(k0)
      hnu0=0
c
c Special case for photonuclear reactions
c
      if (k0.eq.0) then
        p0=1
        h0=0
        pnu0=1
        hnu0=0
        ppi0=0
        hpi0=0
      endif
c
c Initialization for surface interaction
c
c Efermi      : depth of Fermi well
c flagsurface : flag for surface effects in exciton model
c Esurf,Esurf0: well depth for surface interaction
c surface     : well depth for first hole
c Einc        : incident energy in MeV
c
      Esurf=Efermi
      if (flagsurface) then
        if(Esurf0.eq.-1.) then
          if (k0.eq.1) Esurf=surface(1,Einc)
          if (k0.eq.2) Esurf=surface(2,Einc)
        else
          Esurf=Esurf0
        endif
      endif
c
c *********************** Pre-equilibrium model ************************
c
c flagpeout   : flag for output of pre-equilibrium results
c xsflux      : cross section flux
c xsreacinc   : reaction cross section for incident channel
c xsdirdiscsum: total direct cross section
c xsgrsum     : total giant resonance cross section
c flag2comp   : flag for two-component pre-equilibrium model
c exciton     : subroutine for exciton model
c excitonout  : subroutine for output of exciton model parameters
c exciton2    : subroutine for two-component exciton model
c exciton2out : subroutine for output of two-component exciton model
c               parameters
c
c Pre-equilibrium models:
c preeqmode= 1: Exciton model: Analytical solution - matrix element
c preeqmode= 2: Exciton model: Numerical solution - matrix element
c preeqmode= 3: Exciton model: Numerical solution - optical model in
c               transition rates
c preeqmode= 4: Multi-step direct/Multi-step compound
c
      if (flagpeout)
     +  write(*,'(/," ########## PRE-EQUILIBRIUM ##########")')
c
c Correct reaction cross section for direct and giant resonance effects
c
      xsflux=xsreacinc-xsdirdiscsum-xsgrsum
      xsflux=max(xsflux,0.)
c
c Choice between one-component and two-component exciton model.
c
      if (.not.flag2comp) then
        call exciton
        if (flagpeout) call excitonout
      else
        call exciton2
        if (flagpeout) call exciton2out
      endif
c
c Quantum-mechanical pre-equilibrium models.
c
c preeqmode   : designator for pre-equilibrium model
c msd         : subroutine for multi-step direct model
c msdplusmsc  : subroutine for total quantum-mechanical pre-equilibrium
c               cross sections
c flagpecomp  : flag for Kalbach complex particle emission model
c preeqcomplex: subroutine for pre-equilibrium complex particle emission
c preeqcorrect: subroutine to correct pre-equilibrium cross sections for
c               direct effects
c
      if (preeqmode.eq.4) then
        call msd
        call msdplusmsc
      endif
c
c For pre-equilibrium reactions involving complex particles, pickup,
c stripping and knockout contributions are added.
c
      if (flagpecomp) call preeqcomplex
c
c Direct contributions are subtracted from the pre-equilibrium cross
c sections and pre-equilibrium cross sections are collapsed onto
c discrete states.
c
      if (k0.gt.0) call preeqcorrect
c
c The total pre-equilibrium cross sections are assembled.
c
c preeqtotal: subroutine for total pre-equilibrium cross sections
c flagang   : flag for output of angular distributions
c flagddx   : flag for output of double-differential cross sections
c flagrecoil: flag for calculation of recoils
c preeqang  : subroutine for pre-equilibrium angular distributions
c preeqout  : subroutine for output of pre-equilibrium cross sections
c
      call preeqtotal
      if (flagang.or.flagddx.or.flagrecoil) call preeqang
      if (flagpeout) call preeqout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
