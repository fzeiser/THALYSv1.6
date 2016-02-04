      subroutine spectraout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Output of particle spectra
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*20 specfile
      integer      type,nen
c
c ***************************** Spectra ********************************
c
c parskip    : logical to skip outgoing particle
c xsparticle : total particle production cross section
c parname    : name of particle
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c espec      : outgoing energy grid
c xssumout   : cross section summed over mechanisms
c xsdiscout  : total smoothed cross section for discrete state
c xspreeqout : preequilibrium cross section per particle type and
c              outgoing energy
c xsmpreeqout: multiple pre-equilibrium emission spectrum
c xscompout  : compound emission cross section
c flagrecoil : flag for calculation of recoils
c flaglabddx : flag for calculation of DDX in LAB system
c iejlab     : number of ejectile lab bins
c Eejlab     : center of ejectile lab bin
c xsejlab    : LAB ejectile spectrum
c xsejlabint : LAB energy-integrated spectrum
c
      write(*,'(/" 7. Composite particle spectra")')
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (xsparticle(type).eq.0.) goto 10
        write(*,'(/" Spectra for outgoing ",a8/)') parname(type)
        write(*,'("  Energy   Total       Direct    Pre-equil.",
     +    "  Mult. preeq  Compound"/)')
        do 20 nen=ebegin(type),eendout(type)
          write(*,'(1x,f7.3,1p,5e12.5)') espec(type,nen),
     +      xssumout(type,nen),xsdiscout(type,nen),xspreeqout(type,nen),
     +      xsmpreeqout(type,nen),xscompout(type,nen)
   20   continue
        if (flagrecoil.and.flaglabddx) then
          write(*,'(/" LAB spectra for outgoing ",a8/)') parname(type)
          write(*,'("  Energy   Cross section"/)')
          do 30 nen=1,iejlab(type)
            write(*,'(1x,f7.3,1p,e12.5)') Eejlab(type,nen),
     +        xsejlab(type,nen)
   30     continue
          write(*,'(/" Energy-integrated cross section:",1p,e12.5/)')
     +      xsejlabint(type)
        endif
c
c Write results to separate file
c
c filespectrum: designator for spectrum on separate file
c natstring   : string extension for file names
c iso         : counter for isotope
c Einc        : incident energy in MeV
c specfile    : file with spectrum
c parsym      : symbol of particle
c preeqratio  : pre-equilibrium ratio
c
        if (filespectrum(type)) then
          specfile=' spec000.000.tot'//natstring(iso)
          write(specfile(1:1),'(a1)') parsym(type)
          write(specfile(6:12),'(f7.3)') Einc
          write(specfile(6:8),'(i3.3)') int(Einc)
          open (unit=1,status='unknown',file=specfile)
          write(1,'("# ",a1," + ",i3,a2,": ",a8," spectrum")')
     +      parsym(k0),Atarget,Starget,parname(type)
          write(1,'("# E-incident = ",f7.3)') Einc
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') eendout(type)-ebegin(type)+1
          write(1,'("# E-out    Total       Direct    Pre-equil.",
     +      "  Mult. preeq  Compound  Pre-eq ratio")')
          do 40 nen=ebegin(type),eendout(type)
            write(1,'(f7.3,1p,6e12.5)')
     +        espec(type,nen),xssumout(type,nen),xsdiscout(type,nen),
     +        xspreeqout(type,nen),xsmpreeqout(type,nen),
     +        xscompout(type,nen),preeqratio(type,nen)
   40     continue
          close (unit=1)
          if (flagrecoil.and.flaglabddx) then
            specfile(14:16)='lab'
            open (unit=1,status='unknown',file=specfile)
            write(1,'("# ",a1," + ",i3,a2,": ",a8,
     +        " spectrum in LAB frame")') parsym(k0),Atarget,
     +        Starget,parname(type)
            write(1,'("# E-incident = ",f7.3)') Einc
            write(1,'("# ")')
            write(1,'("# # energies =",i3)') iejlab(type)
            write(1,'("# E-out    Total")')
            do 50 nen=1,iejlab(type)
              write(1,'(f7.3,1p,e12.5)') Eejlab(type,nen),
     +          xsejlab(type,nen)
   50       continue
            close (unit=1)
          endif
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
