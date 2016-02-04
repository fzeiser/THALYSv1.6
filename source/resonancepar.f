      subroutine resonancepar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 25, 2012
c | Task  : S-wave resonance parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*4  reschar
      character*90 resfile
      integer      Zix,Nix,Z,A,ia
      real         D0f,dD0f,gamgamf,dgamgamf
c
c ********** Resonance spacings and total radiative widths *************
c
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c ZZ,Z   : charge number of residual nucleus
c AA,A   : mass number of residual nucleus
c reschar: help variable
c resfile: resonance parameter file
c path   : directory containing structure files to be read
c lenpath: length of pathname
c
c Read experimental values from resonance parameter file
c Resonance parameters from the table can always be overruled
c by a value given in the input file.
c
c 1. Inquire whether file is present
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      reschar='z   '
      write(reschar(2:4),'(i3.3)') Z
      resfile=path(1:lenpath)//'resonances/'//reschar
      inquire (file=resfile,exist=lexist)
      if (.not.lexist) goto 30
      open (unit=2,status='old',file=resfile)
c
c 2. Search for the isotope under consideration and read information
c
c ia     : mass number from resonance table
c D0     : experimental s-wave resonance spacing in eV
c dD0    : uncertainty in D0
c gamgam : experimental total radiative width in eV
c dgamgam: uncertainty in gamgam
c D0f,...: help variables
c
   10 read(2,'(4x,i4,2e9.2,10x,2f9.5)',end=20) ia,D0f,dD0f,gamgamf,
     +  dgamgamf
      if (A.ne.ia) goto 10
      if (dD0f.ne.0..and.D0(Zix,Nix).eq.0.) dD0(Zix,Nix)=dD0f*1000.
      if (D0f.ne.0..and.D0(Zix,Nix).eq.0.) D0(Zix,Nix)=D0f*1000.
      if (dgamgamf.ne.0..and.gamgam(Zix,Nix).eq.0.) dgamgam(Zix,Nix)=
     +  dgamgamf
      if (gamgamf.ne.0..and.gamgam(Zix,Nix).eq.0.)
     +  gamgam(Zix,Nix)=gamgamf
   20 close (unit=2)
c
c 3. Tabulated value for total radiative width or systematics
c    (Kopecky, 2002)
c
c gamkopecky  : radiative width in eV by spline fit of Kopecky
c gamgamadjust: adjustable factor for radiative parameters
c             (default 1.)
c
   30 if (gamgam(Zix,Nix).eq.0.) then
        if (A.ge.40.and.A.le.250) then
          gamgam(Zix,Nix)=gamkopecky(A)
        else
          gamgam(Zix,Nix)=min(1593./(A*A),10.)
        endif
      endif
      gamgam(Zix,Nix)=gamgamadjust(Zix,Nix)*gamgam(Zix,Nix)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
