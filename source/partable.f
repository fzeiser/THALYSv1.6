      subroutine partable(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 16, 2013
c | Task  : Write model parameters per nucleus to separate file
c +---------------------------------------------------------------------
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,Z,A,ibar,l
c
c ****************************** Z and A of nucleus ********************
c
c Zix : charge number index for residual nucleus
c Nix : neutron number index for residual nucleus
c ZZ,Z: charge number of residual nucleus
c AA,A: mass number of residual nucleus
c nuc : symbol of nucleus
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      write(11,'("##")')
      write(11,'("## Parameters for ",i3,a2)')  A,nuc(Z)
c
c ********************** Level density parameters **********************
c
c alev         : level density parameter
c aadjust....  : adjustable factors for level density parameters
c                (default 1.)
c gammald      : gamma-constant for asymptotic level density parameter
c pair         : total pairing correction
c D0           : experimental s-wave resonance spacing in eV
c ibar         : fission barrier
c nfisbar      : number of fission barrier parameters
c Pshift       : adjustable pairing shift
c deltaW       : shell correction in nuclear mass
c ldmodel      : level density model
c T            : nuclear temperature
c E0           : constant of temperature formula
c Exmatch      : matching point for Ex
c Ntop         : highest discrete level for temperature matching
c Nlow         : lowest discrete level for temperature matching
c s2adjust     : adjustable constant (Z,A,barrier-dependent) for spin
c                cutoff parameter
c flagcolall   : flag for collective enhancement of level density
c                for all nuclides
c Krotconstant : normalization constant for rotational enhancement
c ctable,ptable: constant to adjust tabulated level densities
c phmodel      : particle-hole state density model
c g            : single-particle level density parameter
c gp           : single-particle proton level density parameter
c gn           : single-particle neutron level density parameter
c
      write(11,'("##")')
      write(11,'("## Level density")')
      write(11,'("##")')
      write(11,'("a              ",2i4,f10.5)') Z,A,alev(Zix,Nix)
      write(11,'("aadjust        ",2i4,f10.5)') Z,A,aadjust(Zix,Nix)
      write(11,'("gammald        ",2i4,f10.5)') Z,A,gammald(Zix,Nix)
      write(11,'("pair           ",2i4,f10.5)') Z,A,pair(Zix,Nix)
      do 10 ibar=0,nfisbar(Zix,Nix)
        write(11,'("Pshift         ",2i4,f10.5,i4)') Z,A,
     +    Pshift(Zix,Nix,ibar),ibar
        write(11,'("deltaW         ",2i4,f10.5,i4)') Z,A,
     +    deltaW(Zix,Nix,ibar),ibar
        if (ldmodel(Zix,Nix).eq.1) then
          write(11,'("T              ",2i4,f10.5,i4)') Z,A,
     +      T(Zix,Nix,ibar),ibar
          write(11,'("E0             ",2i4,f10.5,i4)') Z,A,
     +      E0(Zix,Nix,ibar),ibar
          write(11,'("Exmatch        ",2i4,f10.5,i4)') Z,A,
     +      Exmatch(Zix,Nix,ibar),ibar
          write(11,'("Tadjust        ",2i4,f10.5,i4)') Z,A,
     +      Tadjust(Zix,Nix,ibar),ibar
          write(11,'("E0adjust       ",2i4,f10.5,i4)') Z,A,
     +      E0adjust(Zix,Nix,ibar),ibar
          write(11,'("Exmatchadjust  ",2i4,f10.5,i4)') Z,A,
     +      Exmatchadjust(Zix,Nix,ibar),ibar
        endif
        write(11,'("Ntop           ",2i4,2i4)') Z,A,Ntop(Zix,Nix,ibar),
     +   ibar
        write(11,'("Nlow           ",2i4,2i4)') Z,A,Nlow(Zix,Nix,ibar),
     +   ibar
        write(11,'("s2adjust       ",2i4,f10.5,i4)') Z,A,
     +    s2adjust(Zix,Nix,ibar),ibar
        write(11,'("ctable         ",2i4,f10.5)') Z,A,
     +    ctable(Zix,Nix,ibar)
        write(11,'("ptable         ",2i4,f10.5)') Z,A,
     +    ptable(Zix,Nix,ibar)
        if (flagcolall) write(11,'("Krotconstant   ",2i4,f10.5,i4)')
     +    Z,A,Krotconstant(Zix,Nix,ibar),ibar
   10 continue
      if (D0(Zix,Nix).ne.0.)
     +   write(11,'("D0             ",2i4,1p,e12.5)') Z,A,
     +   D0(Zix,Nix)*0.001
      if (phmodel.eq.1) then
        write(11,'("g              ",2i4,f10.5)') Z,A,g(Zix,Nix)
        write(11,'("gp             ",2i4,f10.5)') Z,A,gp(Zix,Nix)
        write(11,'("gn             ",2i4,f10.5)') Z,A,gn(Zix,Nix)
        write(11,'("gnadjust       ",2i4,f10.5)') Z,A,gnadjust(Zix,Nix)
        write(11,'("gpadjust       ",2i4,f10.5)') Z,A,gpadjust(Zix,Nix)
      endif
c
c ************************ Gamma-ray parameters ************************
c
c gamgam       : total radiative width in eV
c gammax       : number of l-values for gamma multipolarity
c sgr          : strength of GR
c egr          : energy of GR
c ggr          : width of GR
c egradjust....: adjustable factors for giant resonance parameters
c strength     : model for E1 gamma-ray strength function
c etable,ftable: constant to adjust tabulated strength functions
c ngr          : number of GR
c
      write(11,'("##")')
      write(11,'("## Gamma-ray")')
      write(11,'("##")')
      write(11,'("gamgam         ",2i4,f10.5)') Z,A,gamgam(Zix,Nix)
      write(11,'("gamgamadjust   ",2i4,f10.5)') Z,A,
     +  gamgamadjust(Zix,Nix)
      do 110 l=1,gammax
        if (strength.le.2.or.strength.eq.5) then
          write(11,'("sgr            ",2i4,f8.3," E",i1)') Z,A,
     +      sgr(Zix,Nix,1,l,1),l
          write(11,'("egr            ",2i4,f8.3," E",i1)') Z,A,
     +      egr(Zix,Nix,1,l,1),l
          write(11,'("ggr            ",2i4,f8.3," E",i1)') Z,A,
     +      ggr(Zix,Nix,1,l,1),l
          write(11,'("sgradjust      ",2i4,f8.3," E",i1)') Z,A,
     +      sgradjust(Zix,Nix,1,l,1),l
          write(11,'("egradjust      ",2i4,f8.3," E",i1)') Z,A,
     +      egradjust(Zix,Nix,1,l,1),l
          write(11,'("ggradjust      ",2i4,f8.3," E",i1)') Z,A,
     +      ggradjust(Zix,Nix,1,l,1),l
        else
          write(11,'("etable         ",2i4,f10.5)') Z,A,etable(Zix,Nix)
          write(11,'("ftable         ",2i4,f10.5)') Z,A,ftable(Zix,Nix)
        endif
        if (ngr(Zix,Nix,1,l).eq.2) then
          write(11,'("sgr            ",2i4,f8.3," E",i1," 2")') Z,A,
     +      sgr(Zix,Nix,1,l,2),l
          write(11,'("egr            ",2i4,f8.3," E",i1," 2")') Z,A,
     +      egr(Zix,Nix,1,l,2),l
          write(11,'("ggr            ",2i4,f8.3," E",i1," 2")') Z,A,
     +      ggr(Zix,Nix,1,l,2),l
          write(11,'("sgradjust      ",2i4,f8.3," E",i1," 2")') Z,A,
     +      sgradjust(Zix,Nix,1,l,2),l
          write(11,'("egradjust      ",2i4,f8.3," E",i1," 2")') Z,A,
     +      egradjust(Zix,Nix,1,l,2),l
          write(11,'("ggradjust      ",2i4,f8.3," E",i1," 2")') Z,A,
     +      ggradjust(Zix,Nix,1,l,2),l
        endif
        write(11,'("sgr            ",2i4,f8.3," M",i1)') Z,A,
     +    sgr(Zix,Nix,0,l,1),l
        write(11,'("egr            ",2i4,f8.3," M",i1)') Z,A,
     +    egr(Zix,Nix,0,l,1),l
        write(11,'("ggr            ",2i4,f8.3," M",i1)') Z,A,
     +    ggr(Zix,Nix,0,l,1),l
        write(11,'("sgradjust      ",2i4,f8.3," M",i1)') Z,A,
     +    sgradjust(Zix,Nix,0,l,1),l
        write(11,'("egradjust      ",2i4,f8.3," M",i1)') Z,A,
     +    egradjust(Zix,Nix,0,l,1),l
        write(11,'("ggradjust      ",2i4,f8.3," M",i1)') Z,A,
     +    ggradjust(Zix,Nix,0,l,1),l
  110 continue
c
c ************************** Fission parameters ************************
c
c flagfission: flag for fission
c nfisbar    : number of fission barrier parameters
c fbarrier   : height of fission barrier
c fbaradjust : adjustable factors for fission parameters
c fwidth     : width of fission barrier
c fismodelx  : fission model
c widthc2    : width of class2 states
c Rtransmom  : normalization constant for moment of inertia for
c              transition states
c Rclass2mom : normalization constant for moment of inertia for
c              class 2 states
c
      if (flagfission) then
        write(11,'("##")')
        write(11,'("## Fission parameters")')
        write(11,'("##")')
        do 210 ibar=1,nfisbar(Zix,Nix)
          if (fismodelx(Zix,Nix).ne.5) then
            write(11,'("fisbar         ",2i4,f10.5,i3)') Z,A,
     +        fbarrier(Zix,Nix,ibar),ibar
            write(11,'("fishw          ",2i4,f10.5,i3)') Z,A,
     +        fwidth(Zix,Nix,ibar),ibar
            write(11,'("fisbaradjust   ",2i4,f10.5,i3)') Z,A,
     +        fbaradjust(Zix,Nix,ibar),ibar
            write(11,'("fishwadjust    ",2i4,f10.5,i3)') Z,A,
     +        fwidthadjust(Zix,Nix,ibar),ibar
          endif
          if (fismodelx(Zix,Nix).eq.5.and.ibar.eq.1) then
            write(11,'("betafiscor     ",2i4,f10.5)') Z,A,
     +        betafiscor(Zix,Nix)
            write(11,'("vfiscor        ",2i4,f10.5)') Z,A,
     +        vfiscor(Zix,Nix)
          endif
          if (ibar.lt.nfisbar(Zix,Nix))
     +      write(11,'("class2width    ",2i4,f10.5,i3)') Z,A,
     +      widthc2(Zix,Nix,ibar),ibar
          write(11,'("Rtransmom      ",2i4,f10.5,i3)') Z,A,
     +      Rtransmom(Zix,Nix,ibar),ibar
          write(11,'("Rclass2mom     ",2i4,f10.5,i3)') Z,A,
     +      Rclass2mom(Zix,Nix,ibar),ibar
  210   continue
      endif
      write(11,'("##--------------------------------------------")')
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
