      subroutine thermalxs
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 17, 2008
c | Task  : Cross sections at thermal energies
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*4  thchar
      character*90 thfile
      integer      Z,A,ia
      real         xs,xsp,xsalpha,ald,Spair
c
c ********** Resonance spacings and total radiative widths *************
c
c ZZ,Z   : charge number of residual nucleus
c AA,A   : mass number of residual nucleus
c thchar : help variable
c thfile : thermal cross section file
c path   : directory containing structure files to be read
c lenpath: length of pathname
c
c Read experimental values from thermal cross section file
c Values from the table can always be overruled by values given in the
c input file.
c
c 1. Inquire whether file is present
c
      Z=ZZ(0,0,1)
      A=AA(0,0,1)
      thchar='z   '
      write(thchar(2:4),'(i3.3)') Z
      thfile=path(1:lenpath)//'thermal/'//thchar
      inquire (file=thfile,exist=lexist)
      if (.not.lexist) goto 30
      open (unit=2,status='old',file=thfile)
c
c 2. Search for the isotope under consideration and read information
c
c ia          : mass number from resonance table
c xscaptherm  : thermal capture cross section
c xsptherm    : thermal (n,p) cross section
c xsalphatherm: thermal (n,a) cross section
c xs,....     : help variables
c
      xs=0.
      xsp=0.
      xsalpha=0.
   10 read(2,'(4x,i4,8x,3(e9.2,9x))',end=20) ia,xs,xsp,xsalpha
      if (A.ne.ia) goto 10
      if (xs.ne.0..and.xscaptherm.eq.0.) xscaptherm=xs
      if (xsp.ne.0..and.xsptherm.eq.0.) xsptherm=xsp
      if (xsalpha.ne.0..and.xsalphatherm.eq.0.) xsalphatherm=xsalpha
   20 close (unit=2)
c
c 2. Systematics
c
c Kopecky's value for (n,gamma) cross section at thermal energy.
c (J. Kopecky, M.G. Delfini, H.A.J. van der Kamp and D. Nierop:
c Revisions and extensions of neutron capture cross-sections in
c the European Activation File EAF-3, ECN-C--92-051, July 1992.)
c
c alev,ald: level density parameter
c Spair   : help variable
c S       : separation energy per particle
c pair    : pairing energy
c
   30 if (xscaptherm.eq.0.) then
        ald=alev(0,0)
        Spair=S(0,0,1)-pair(0,0)
        Spair=max(Spair,1.)
        xscaptherm=1.5e-3*(ald*Spair)**3.5
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
