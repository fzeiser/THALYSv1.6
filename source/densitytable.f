      subroutine densitytable(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Marieke Duijvestijn
c | Date  : June 2, 2012
c | Task  : Tabulated level densities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist
      character*4      denchar
      character*90     denfile
      integer          Zix,Nix,Z,A,ibar,nloop,ldmod,ploop,ia,parity,
     +                 nex,J
      real             Ktriax,Eex,ald,ignatyuk,spincut,term
      double precision pardisloc,ldtot,ld2j1(0:numJ)
c
c *********** Tabulated level densities from Goriely *******************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c AA,A,ia    : mass number of residual nucleus
c flagfission: flag for fission
c nfisbar    : number of fission barrier parameters
c ldmodel    : level density model
c nloop,ploop: help variables
c pardisloc  : variable to account for parity distribution
c nendens    : number of energies for level density grid
c ldexist    : flag for existence of level density table
c denfile    : level density parameter file
c parity     : parity
c ldtable    : level density from table
c ldtottableP: total level density per parity from table
c ldtottable : total level density from table
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      denchar='z   '
      write(denchar(2:4),'(i3.3)') Z
      nloop=0
      if (flagfission) nloop=nfisbar(Zix,Nix)
      ldmod=ldmodel(Zix,Nix)
      if (ldmod.eq.5) then
        ploop=-1
        pardisloc=1.
      else
        ploop=1
        pardisloc=0.5
      endif
      do 10 ibar=0,nloop
        ldexist(Zix,Nix,ibar)=.false.
        denfile='                                                      '
c
c Ground state
c
        if (ibar.eq.0) then
          if (ldmod.eq.4) then
            denfile=
     +       path(1:lenpath)//'density/ground/goriely/'//denchar//'.tab'
          endif
          if (ldmod.eq.5) then
            denfile=
     +       path(1:lenpath)//'density/ground/hilaire/'//denchar//'.tab'
          endif
          if (ldmod.eq.6) then
            denfile=path(1:lenpath)//'density/ground/hilaireD1M/'//
     +        denchar//'.tab'
          endif
        endif
c
c First barrier
c
        if (ibar.eq.1) then
          if (ldmod.eq.4) then
            denfile=path(1:lenpath)//'density/fission/goriely/inner/'
     +        //denchar
          else
            denfile=path(1:lenpath)//'density/fission/hilaire/Max1/'
     +        //denchar
          endif
        endif
c
c Second barrier
c
        if (ibar.eq.2) then
          if (ldmod.eq.4) then
            denfile=path(1:lenpath)//'density/fission/goriely/outer/'
     +        //denchar
          else
            denfile=path(1:lenpath)//'density/fission/hilaire/Max2/'
     +        //denchar
          endif
        endif
c
c Third barrier
c
        if (ibar.eq.3) then
          if (ldmod.eq.5) then
            denfile=path(1:lenpath)//'density/fission/hilaire/Max3/'
     +        //denchar
          else
            goto 10
          endif
        endif
c
c Check existence of file and read data from the tables.
c
        inquire (file=denfile,exist=lexist)
        if (lexist) then
          open (unit=2,status='old',file=denfile)
          do 20 parity=1,ploop,-2
   30       read(2,'(/31x,i3//)',end=80) ia
            if (A.ne.ia) then
              do 40 nex=1,nendens(Zix,Nix)+1
                read(2,'()')
   40         continue
              goto 30
            else
              ldexist(Zix,Nix,ibar)=.true.
              do 50 nex=1,nendens(Zix,Nix)
                read(2,'(24x,e9.2,9x,30e9.2)',err=100)
     +            ldtot,(ld2j1(J),J=0,29)
c
c Determination of the mass-asymmetric enhancement factor for fission
c barrier
c
c Ktriax  : level density enhancement factor for triaxial shapes
c ald     : level density parameter
c ignatyuk: function for energy dependent level density parameter a
c twopi   : 2.*pi
c spincut : spin cutoff factor
c
                Ktriax=1.
                if (ibar.gt.0) then
                  if (axtype(Zix,Nix,ibar).eq.2) Ktriax=2.
                  if (axtype(Zix,Nix,ibar).ge.3) then
                    Eex=edens(nex)
                    ald=ignatyuk(Zix,Nix,Eex,ibar)
                    term=sqrt(spincut(Zix,Nix,ald,Eex,ibar))
                    if (axtype(Zix,Nix,ibar).eq.3)
     +                Ktriax=0.5*sqrt(twopi)*term
                    if (axtype(Zix,Nix,ibar).eq.4)
     +                Ktriax=sqrt(twopi)*term
                    if (axtype(Zix,Nix,ibar).eq.5)
     +                Ktriax=2.*sqrt(twopi)*term
                  endif
                endif
                ldtottableP(Zix,Nix,nex,parity,ibar)=
     +            pardisloc*ldtot*Ktriax
                if (ploop.eq.1) ldtottableP(Zix,Nix,nex,-1,ibar)=
     +            pardisloc*ldtot*Ktriax
                ldtottable(Zix,Nix,nex,ibar)=
     +            ldtottable(Zix,Nix,nex,ibar)+ldtot
                do 60 J=0,29
                  ldtable(Zix,Nix,nex,J,parity,ibar)=
     +              pardisloc*ld2j1(J)*Ktriax
   60           continue
                if (ploop.eq.1) then
                  do 70 J=0,29
                    ldtable(Zix,Nix,nex,J,-1,ibar)=
     +                pardisloc*ld2j1(J)*Ktriax
   70             continue
                endif
   50         continue
              read(2,'()')
            endif
   20     continue
   80     close (unit=2)
        endif
c
c Special case: make parity-independent level densities from
c parity-dependent tables (e.g. for testing the impact of
c parity-dependence).
c
        if (ldmod.eq.5.and.ldexist(Zix,Nix,ibar).and..not.flagparity)
     +    then
          do 110 nex=1,nendens(Zix,Nix)
            ldtottableP(Zix,Nix,nex,1,ibar)=0.5*
     +        (ldtottableP(Zix,Nix,nex,-1,ibar)+
     +        ldtottableP(Zix,Nix,nex,1,ibar))
            do 120 J=0,29
              ldtable(Zix,Nix,nex,J,1,ibar)=0.5*
     +          (ldtable(Zix,Nix,nex,J,-1,ibar)+
     +          ldtable(Zix,Nix,nex,J,1,ibar))
  120       continue
  110     continue
        endif
   10 continue
      return
  100 write(*,'(" TALYS-error: Wrong level density table for",
     +  " Z=",i3," A=",i3)') Z,A
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
