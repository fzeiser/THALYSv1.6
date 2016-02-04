      subroutine totalout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 19, 2012
c | Task  : Output of total cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*18 totfile
      integer      nen
c
c ********************* Total cross sections ***************************
c
c Einc        : incident energy in MeV
c nin         : counter for incident energy
c eninccm     : center-of-mass incident energy in MeV
c k0          : index of incident particle
c xstotinc    : total cross section (neutrons only) for incident channel
c xselasinc   : total elastic cross section (neutrons only) for incident
c               channel
c xsreacinc   : reaction cross section for incident channel
c xscompel    : compound elastic cross section
c xsnonel     : non-elastic cross section
c xsdirdiscsum: total direct cross section
c xspreeqsum  : total preequilibrium cross section summed over particles
c flaggiant   : flag for collective contribution from giant resonances
c xsgrsum     : sum over giant resonance cross sections
c xscompnonel : total compound non-elastic cross section
c xselastot   : total elastic cross section (shape + compound)
c
      if (Einc.ge.0.001) then
        write(*,'(/" ########### REACTION SUMMARY FOR E=",f9.5,
     +    " ###########"/)') Einc
      else
        write(*,'(/" ########### REACTION SUMMARY FOR E=",1p,e12.5,
     +    " ###########"/)') Einc
      endif
      write(*,'(" Center-of-mass energy: ",f7.3/)') eninccm
      write(*,'(" 1. Total (binary) cross sections"/)')
      if (k0.eq.1) write(*,'(" Total           =",1p,e12.5)') xstotinc
      if (k0.eq.1) write(*,'("   Shape elastic   =",1p,e12.5)')
     +  xselasinc
      write(*,'("   Reaction        =",1p,e12.5)') xsreacinc
      write(*,'("     Compound elastic=",1p,e12.5)') xscompel
      write(*,'("     Non-elastic     =",1p,e12.5)') xsnonel
      write(*,'("       Direct          =",1p,e12.5)') xsdirdiscsum
      write(*,'("       Pre-equilibrium =",1p,e12.5)') xspreeqsum
      if (flagracap) write(*,'("       Direct Capture  =",1p,e12.5)')
     +  xsracape
      if (flaggiant) write(*,'("       Giant resonance =",1p,e12.5)')
     +  xsgrsum
      write(*,'("       Compound non-el =",1p,e12.5)') xscompnonel
      if (k0.eq.1) write(*,'("     Total elastic   =",1p,e12.5)')
     +  xselastot
c
c Write results to separate file
c
c filetotal: flag for total cross sections on separate file
c natstring: string extension for file names
c iso      : counter for isotope
c numinclow: number of incident energies below Elow
c parsym   : symbol of particle
c k0       : index of incident particle
c Atarget  : mass number of target nucleus
c Starget  : symbol of target nucleus
c Ztarget  : charge number of target nucleus
c numinc   : number of incident energies
c
      if (filetotal) then
        totfile='total.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,status='unknown',file=totfile)
          write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') numinc
          write(1,'("#    E      Non-elastic  Elastic     Total",
     +      "     Comp. el.  Shape el.  Reaction",
     +      " Comp. nonel   Direct   Pre-equil.  Dir. Capt.")')
          do 10 nen=1,numinclow
            write(1,'(1p,e10.3,2x,9e11.4)') eninc(nen),fxsnonel(nen),
     +        fxselastot(nen),fxstotinc(nen),fxscompel(nen),
     +        fxselasinc(nen),fxsreacinc(nen),fxscompnonel(nen),
     +        fxsdirdiscsum(nen),fxspreeqsum(nen)
   10     continue
        else
          open (unit=1,status='old',file=totfile)
          do 20 nen=1,nin+4
            read(1,*,end=30,err=30)
   20     continue
        endif
        write(1,'(1p,e10.3,2x,10e11.4)') Einc,xsnonel,xselastot,
     +    xstotinc,xscompel,xselasinc,xsreacinc,xscompnonel,
     +    xsdirdiscsum,xspreeqsum,xsracape
   30   close (unit=1)
c
c Total cross sections (i.e. from OMP) only
c
        totfile='totalxs.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,status='unknown',file=totfile)
          write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') numinc
          write(1,'("#    E      Cross section")')
          do 40 nen=1,numinclow
            write(1,'(1p,e10.3,2x,e11.4)') eninc(nen),fxstotinc(nen)
   40     continue
        else
          open (unit=1,status='old',file=totfile)
          do 50 nen=1,nin+4
            read(1,*,end=60,err=60)
   50     continue
        endif
        write(1,'(1p,e10.3,2x,e11.4)') Einc,xstotinc
   60   close (unit=1)
c
c Elastic cross sections only
c
        totfile='elastic.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,status='unknown',file=totfile)
          write(1,'("# ",a1," + ",i3,a2," Elastic cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') numinc
          write(1,'("#    E      Cross section")')
          do 70 nen=1,numinclow
            write(1,'(1p,e10.3,2x,e11.4)') eninc(nen),fxselastot(nen)
   70     continue
        else
          open (unit=1,status='old',file=totfile)
          do 80 nen=1,nin+4
            read(1,*,end=90,err=90)
   80     continue
        endif
        write(1,'(1p,e10.3,2x,e11.4)') Einc,xselastot
   90   close (unit=1)
c
c Nonelastic cross sections only
c
        totfile='nonelastic.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,status='unknown',file=totfile)
          write(1,'("# ",a1," + ",i3,a2," Nonelastic cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') numinc
          write(1,'("#    E      Cross section")')
          do 100 nen=1,numinclow
            write(1,'(1p,e10.3,2x,e11.4)') eninc(nen),fxsnonel(nen)
  100     continue
        else
          open (unit=1,status='old',file=totfile)
          do 110 nen=1,nin+4
            read(1,*,end=120,err=120)
  110     continue
        endif
        write(1,'(1p,e10.3,2x,e11.4)') Einc,xsnonel
  120   close (unit=1)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
