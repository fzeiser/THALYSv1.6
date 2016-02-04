      subroutine phdensitytablejp(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely & Arjan Koning
c | Date  : December 23, 2013
c | Task  : Tabulated spin- and parity-dependent particle-hole level
c |         densities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist
      character*4      denchar
      character*6      phdir
      character*90     denfile
      integer          Zix,Nix,Z,A,ia,parity,nex,J
      real             ephjpgrid
      double precision ldtot,ld2j1(0:numJ)
c
c *********** Tabulated level densities from Goriely *******************
c
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c ZZ,Z         : charge number of residual nucleus
c AA,A,ia      : mass number of residual nucleus
c ldmodelracap : level density model for direct radiative captures
c nenphdens    : number of energies for level density grid
c denfile      : level density parameter file
c parity       : parity
c edensphjp    : energy grid of ph spin- and parity-dependent level density from table
c phdensjp     : ph spin- and parity-dependent level density from table
c phdenstot    : total ph level density from table
c
      if (ldmodelracap.ne.1) return
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      if (k0.eq.1) then
        phdir='ph0011'
      else
        ldmodelracap=2
        return
      endif
ctest if (k0.eq.2) phdir='ph1100'
ctest if (k0.eq.3) phdir='ph1111'
ctest if (k0.eq.4) phdir='ph1122'
ctest if (k0.eq.5) phdir='ph2211'
ctest if (k0.eq.6) phdir='ph2222'
      denchar='z   '
      write(denchar(2:4),'(i3.3)') Z
      denfile=path(1:lenpath)//'density/phjp/'//phdir//'/'
     +  //denchar
c
c Check existence of file and read data from the tables.
c
      inquire (file=denfile,exist=lexist)
      if (lexist) then
        open (unit=2,status='old',file=denfile)
        do 10 parity=1,-1,-2
   20     read(2,'(/31x,i3//)',end=10,err=10) ia
          if (A.ne.ia) then
            do 30 nex=1,nenphdens+1
              read(2,'()')
   30       continue
            goto 20
          else
            do 40 nex=1,nenphdens
              read(2,'(1x,f6.2,17x,e9.2,9x,30e9.2)',err=100)
     +          ephjpgrid,ldtot,(ld2j1(J),J=0,29)
c
c Determination of the mass-asymmetric enhancement factor for fission
c barrier
c
              phdenstot(Zix,Nix,nex)=phdenstot(Zix,Nix,nex)+ldtot
              edensphjp(Zix,Nix,nex)=ephjpgrid
              do 50 J=0,29
                phdensjp(Zix,Nix,nex,J,parity)=ld2j1(J)
   50         continue
   40       continue
            read(2,'()')
          endif
   10   continue
        close (unit=2)
      else
        write(*,'("Input ph file not available:",a90)') denfile
        write(*,'("  Change of ldmodelracap=1 to ldmodelracap=2")')
        ldmodelracap=2
      endif
c
      return
  100 write(*,'(" TALYS-error: Wrong ph level density table for",
     +  " Z=",i3," A=",i3)') Z,A
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
