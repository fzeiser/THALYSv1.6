      subroutine endfinfo
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 25, 2012
c | Task  : Info for ENDF-6 file
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      flagel
      character*1  yesno
      character*40 endfenergyfile
      integer      nen
c
c ****************** Reaction information for ENDF-6 file **************
c
c k0         : index of incident particle
c Ztarget    : charge number of target nucleus
c Atarget    : mass number of target nucleus
c tarmass    : mass of target nucleus
c targetspin : spin of target
c targetE    : energy of target
c Ltarget    : excited level of target
c Liso       : isomeric number of target
c numinc     : number of incident energies
c eninc      : incident energy in MeV
c nlevmax    : maximum number of included discrete levels for target
c energyfile : file with incident energies
c yesno      : y or n function
c flagendfdet: flag for detailed ENDF-6 information per channel
c Nrescue    : number of energies for adjustment factors
c flagrecoil : flag for calculation of recoils
c flagurrendf: flag for URR info to ENDF
c
      open (unit=1,status='unknown',file='tefal.inf')
      write(1,'(i3,"       : projectile type")') k0
      write(1,'(i3,"       : Z of target")') Ztarget
      write(1,'(i3,"       : A of target")') Atarget
      write(1,'(f10.6,": mass of target in a.m.u.")') tarmass
      write(1,'(f4.1,"      : spin of target")') targetspin
      write(1,'(f10.6,": energy of target")') targetE
      write(1,'(2i4,"  : level and isomeric number of target")')
     +  Ltarget,Liso
      endfenergyfile='energies.endf'
      open (unit=2,status='unknown',file=endfenergyfile)
      do 10 nen=1,numinc
        write(2,'(1pe12.5)') eninc(nen)
   10 continue
      close (unit=2)
      write(1,'(i3,"       : number of incident energies")') numinc
      write(1,'(i3,"       : number of discrete levels")') nlevmax
      write(1,'(a40,": file with incident energies")') endfenergyfile
      write(1,'(a1,"         : detailed ENDF-6 information",
     + " per channel")') yesno(flagendfdet)
      if (Nrescue(1,-1).ne.0) then
        flagel=.false.
      else
        flagel=.true.
      endif
      write(1,'(a1,"         : keep elastic cross section",
     + " in normalization")') yesno(flagel)
      write(1,'(a1,"         : recoils")') yesno(flagrecoil)
      write(1,'(a1,"         : urr")') yesno(flagurrendf)
      close (unit=1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
