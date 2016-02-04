      subroutine phdensitytable(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 17, 2009
c | Task  : Tabulated particle-hole state densities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*4  denchar
      character*90 denfile
      integer      Zix,Nix,Z,A,ia,nex,i,nex2
      real         phden2(0:numdens,numconf),phden(0:numdens,numconf),
     +             sum
c
c ******** Tabulated particle-hole level densities from Hilaire ********
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c AA,A,ia    : mass number of residual nucleus
c nenphdens  : number of energies for particle-hole state density grid
c phexist1   : flag for existence of particle-hole state density table
c phexist2   : flag for existence of particle-hole state density table
c denfile    : level density parameter file
c flag2comp  : flag for two-component pre-equilibrium model
c phtable2   : particle-hole state density from table
c edens      : energy grid for tabulated level densities
c phtable1   : particle-hole state density from table
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      denchar='z   '
      write(denchar(2:4),'(i3.3)') Z
      denfile=path(1:lenpath)//'density/ph/'//denchar
c
c Check existence of file and read data from the tables.
c
      inquire (file=denfile,exist=lexist)
      if (.not.lexist) return
      open (unit=2,status='old',file=denfile)
   20 read(2,'(/31x,i3)',end=10) ia
      if (A.ne.ia) then
        do 30 nex=1,nenphdens+5
          read(2,'()')
   30   continue
        goto 20
      else
        if (flag2comp) then
          do 40 i=1,Nphconf2
            phexist2(Zix,Nix,ppitable(i),hpitable(i),pnutable(i),
     +        hnutable(i))=.true.
   40     continue
        else
          do 50 i=1,Nphconf1
            phexist1(Zix,Nix,pptable(i),hhtable(i))=.true.
   50     continue
        endif
        read(2,'(///)')
        do 60 nex=1,nenphdens
          read(2,'(7x,86(e9.2))',err=200)
     +      (phden2(nex,i),i=1,Nphconf2),(phden(nex,i),i=1,Nphconf1)
   60   continue
        read(2,'()')
        if (flag2comp) then
          do 70 nex=1,nenphdens
            do 80 i=1,Nphconf2
              sum=0.
              do 90 nex2=max(nex-2,1),min(nex+2,nenphdens)
                sum=sum+phden2(nex2,i)*0.5*
     +            (edens(min(nex2+1,nenphdens))-edens(max(nex2-1,1)))
   90         continue
              phtable2(Zix,Nix,ppitable(i),hpitable(i),pnutable(i),
     +          hnutable(i),nex)=sum/(edens(min(nex+2,nenphdens))-
     +          edens(max(nex-2,1)))
   80       continue
   70     continue
        else
          do 100 nex=1,nenphdens
            do 110 i=1,Nphconf1
              sum=0.
              do 120 nex2=max(nex-2,1),min(nex+2,nenphdens)
                sum=sum+phden(nex2,i)*0.5*
     +            (edens(min(nex2+1,nenphdens))-edens(max(nex2-1,1)))
  120         continue
              phtable1(Zix,Nix,pptable(i),hhtable(i),nex)=
     +          sum/(edens(min(nex+2,nenphdens))-edens(max(nex-2,1)))
  110       continue
  100     continue
        endif
      endif
   10 return
  200 write(*,'(" TALYS-error: Wrong particle-hole density table for",
     +  " Z=",i3," A=",i3)') Z,A
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
