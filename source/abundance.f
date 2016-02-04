      subroutine abundance
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Natural abundances
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*1  ch
      character*4  abchar
      character*90 abfile
      integer      i,i2,ia
      real         ab,abtot
c
c ****************************** Abundances ****************************
c
c nlines : number of input lines
c inline : input line
c ch     : character
c abfile : isotopic abundance file
c abchar : help variable
c Ztarget: charge number of target nucleus
c Starget: symbol of target nucleus
c path   : directory containing structure files to be read
c lenpath: length of pathname
c ia     : mass number from abundance table
c abun,ab: isotopic abundance
c isotope: isotope number of residual nucleus
c isonum : number of isotopes
c abtot  : summed abundances for normalization
c nuc    : symbol of nucleus
c
c Note that for non-natural elements we take the longest-lived isotope
c as default.
c
c 1. Isotopic abundances from user file
c
c
      do 10 i=1,nlines
        if (inline(i)(1:10).eq.'abundance ') then
          do 20 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=200) abfile
              goto 100
            endif
   20     continue
        endif
   10 continue
c
c 2. Isotopic abundances from abundance directory
c
      abchar='z   '
      write(abchar(2:4),'(i3.3)') Ztarget
      abfile=path(1:lenpath)//'abundance/'//abchar
      inquire (file=abfile,exist=lexist)
      if (.not.lexist) goto 200
c
c Read abundances from file
c
  100 open (unit=2,status='old',file=abfile)
      i=1
  110 read(2,'(4x,i4,f11.6)',end=120,err=210) ia,ab
      isotope(i)=ia
      abun(i)=0.01*ab
      i=i+1
      goto 110
  120 close (unit=2)
      isonum=i-1
c
c Normalize abundances to 1.
c
      abtot=0.
      do 130 i=1,isonum
        abtot=abtot+abun(i)
  130 continue
      do 140 i=1,isonum
        abun(i)=abun(i)/abtot
  140 continue
      write(*,'(/" Calculation for multi-isotope case"/)')
      write(*,'("  Isotope Abundance"/)')
      do 150 i=1,isonum
        write(*,'(2x,i3,a2,f11.6)') isotope(i),Starget,abun(i)
  150 continue
c
c Create file name extensions
c
c numiso   : maximum number of isotopes per element
c natstring: string extension for file names
c
      do 160 i=1,isonum
        write(natstring(i)(1:4),'(".",i3.3)') isotope(i)
  160 continue
      return
  200 write(*,'(" TALYS-error: No natural isotopes for this",
     +  " element, the mass keyword must be different from 0")')
      stop
  210 write(*,'(" TALYS-error: Format error in abundance file ",a72)')
     +  abfile
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
