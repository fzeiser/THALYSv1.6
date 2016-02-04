      subroutine binaryout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 19, 2012
c | Task  : Output of binary cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*10 binfile
      integer      type,nen
      real         xsrac,xsc
c
c ******************* Binary non-elastic channels **********************
c
c flagcompo   : flag for output of cross section components
c parskip     : logical to skip outgoing particle
c parname     : name of particle
c xsracape    : direct radiative capture cross section
c xsbinary    : cross section from initial compound to residual nucleus
c xsdirdisctot: direct cross section summed over discrete states
c xspreeqtot  : preequilibrium cross section per particle type
c xsgrtot     : total smoothed giant resonance cross section
c               incident channel
c xscomp      : compound cross section per particle type
c
      write(*,'(/" 2. Binary non-elastic cross sections ",
     +  "(non-exclusive)")')
      if (flagcompo) then
        write(*,'(36x," Direct  Preequilibrium Compound  Dir. Capt.")')
      else
        write(*,'()')
      endif
      do 10 type=-1,6
        if (parskip(type)) goto 10
        if (type.eq.0) then
          xsrac=xsracape
        else
          xsrac=0.
        endif
        if (flagcompo.and.type.ge.0) then
          xsc=max(xsbinary(type)-xsdirdisctot(type)-xspreeqtot(type)-
     +      xsgrtot(type)-xsrac,0.)
          write(*,'(1x,a8,"=",1p,e12.5,12x,4e12.5)') parname(type),
     +      xsbinary(type),xsdirdisctot(type),
     +      xspreeqtot(type)+xsgrtot(type),xsc,xsrac
        else
          write(*,'(1x,a8,"=",1p,e12.5)') parname(type),xsbinary(type)
        endif
   10 continue
c
c Write results to separate file
c
c filetotal : flag for total cross sections on separate file
c numinclow : number of incident energies below Elow
c eninc,Einc: incident energy in MeV
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c Ztarget   : charge number of target nucleus
c numinc    : number of incident energies
c
      if (filetotal) then
        binfile='binary.tot'
        if (nin.eq.numinclow+1) then
          open (unit=1,status='unknown',file=binfile)
          write(1,'("# ",a1," + ",i3,a2," Binary cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') numinc
          write(1,'("#    E       ",7(2x,a8,1x))')
     +      (parname(type),type=0,6)
          do 20 nen=1,numinclow
            write(1,'(1p,e10.3,2x,7e11.4)') eninc(nen),
     +        (fxsbinary(nen,type),type=0,6)
  20     continue
        else
          open (unit=1,status='old',file=binfile)
          do 30 nen=1,nin+4
            read(1,*,end=40,err=40)
  30     continue
        endif
        write(1,'(1p,e10.3,2x,7e11.4)') Einc,(xsbinary(type),type=0,6)
  40    close (unit=1)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
