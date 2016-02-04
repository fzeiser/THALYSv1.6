      subroutine fissionpar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire, Marieke Duijvestijn and Arjan Koning
c | Date  : February 16, 2012
c | Task  : Fission parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*4  fischar
      character*90 fisfile,hbsfile,c2file
      integer      Zix,Nix,fislocal,Z,N,A,ia,i,j,il,modz,modn,nbar
      real         bar1,bar2,hw1,hw2,egs,lbar0,esp,bb,vv
c
c ****************** Read fission barrier parameters *******************
c
c Note that next to the chosen fission model (fismodel), there is always
c an alternative fission model (fismodelalt) which comes into play if
c fission parameters for the first choice model are not available.

c Determine whether barrier parameters have been provided in input
c
c Zix      : charge number index for residual nucleus
c Nix      : neutron number index for residual nucleus
c nfisbar  : number of fission barrier parameters
c fbarrier : height of fission barrier
c fwidth   : width of fission barrier
c fislocal : fission model
c fismodelx: fission model
c
      nfisbar(Zix,Nix)=0
      do 10 i=1,numbar
        if (fbarrier(Zix,Nix,i).ne.0..or.fwidth(Zix,Nix,i).ne.0.)
     +    nfisbar(Zix,Nix)=nfisbar(Zix,Nix)+1
   10 continue
   20 fislocal=fismodelx(Zix,Nix)
c
c Fission parameters from database
c
c ZZ,Z   : charge number of residual nucleus
c NN,N   : neutron number of residual nucleus
c AA,A   : mass number of residual nucleus
c fischar: help variable
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      fischar='z   '
      write(fischar(2:4),'(i3.3)') Z
c
c Fismodel 1: Experimental parameters
c
c fisfile     : fission file
c path        : directory containing structure files to be read
c lenpath     : length of pathname
c axtype      : type of axiality of barrier
c                  1: axial symmetry
c                  2: left-right asymmetry
c                  3: triaxial and left-right symmetry
c                  4: triaxial no left-right symmetry
c                  5: no symmetry
c bar1,bar2   : inner and outer barrier heights
c hw1,hw2     : inner and outer barrier curvatures
c
c Fission barriers from database may have been overruled by user input.
c As starting point, we take the RIPL values as compiled by Maslov.
c
      if (fislocal.eq.1) then
        fisfile=path(1:lenpath)//'fission/barrier/'//fischar
        inquire (file=fisfile,exist=lexist)
        if (.not.lexist) goto 120
        open (unit=2,status='old',file=fisfile)
  110   read(2,'(4x,i4,4x,1x,2(f8.2),5x,2(f8.2))',end=120)
     +    ia,bar1,hw1,bar2,hw2
        if (A.ne.ia) goto 110
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=bar1
        if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=hw1
        if (fbarrier(Zix,Nix,2).eq.0.) fbarrier(Zix,Nix,2)=bar2
        if (fwidth(Zix,Nix,2).eq.0.) fwidth(Zix,Nix,2)=hw2
        if (nfisbar(Zix,Nix).ne.3) nfisbar(Zix,Nix)=2
        close (unit=2)
  120   if (nfisbar(Zix,Nix).eq.0) fislocal=2
      endif
c
c Fismodel 2: Mamdouh parameters
c
      if (fislocal.eq.2) then
        fisfile=path(1:lenpath)//'fission/mamdouh/'//fischar
        inquire (file=fisfile,exist=lexist)
        if (.not.lexist) goto 100
        open (unit=2,status='old',file=fisfile)
  210   read(2,'(4x,i4,2(24x,f8.2))',end=100) ia,bar1,bar2
        if (A.ne.ia) goto 210
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=bar1
        if (fbarrier(Zix,Nix,2).eq.0.) fbarrier(Zix,Nix,2)=bar2
        if (fbarrier(Zix,Nix,1).eq.0..or.fbarrier(Zix,Nix,2).eq.0.) then
          nfisbar(Zix,Nix)=1
        else
          nfisbar(Zix,Nix)=2
        endif
        close (unit=2)
      endif
  100 if (fismodelx(Zix,Nix).le.2.and.nfisbar(Zix,Nix).eq.0)
     +  fislocal=fismodelalt
c
c Fismodel 3: Sierk
c
c barsierk: subroutine for fission barrier heights, rotating gs energy
c           and lbar0
c egs     : rotating ground state energy
c lbar0   : l-value for which barrier height becomes zero
c il      : angular momentum
c
      if (fislocal.eq.3) then
        il=0
        nfisbar(Zix,Nix)=1
        call barsierk(Z,A,il,bar1,egs,lbar0)
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=bar1
        if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=0.24
      endif
c
c Fismodel 4: Rotating Liquid Drop Model
c
c rldm: subroutine for saddle point energies, rotating gs energy
c esp : saddle point energy
c
      if (fislocal.eq.4) then
        il=0
        nfisbar(Zix,Nix)=1
        call rldm(Z,A,il,egs,esp)
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=esp-egs
        if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=0.24
      endif
c
c Fismodel 5: WKB approximation
c
c Read the potential energy curve and call the WKB subroutine
c
c nbeta     : number of beta values
c betafiscor: adjustable factor for fission path width
c vfiscor   : adjustable factor for fission path height
c betafis   : fission path width
c vfis      : fission path height
c nbinswkb  : integration step for WKB calculation
c nbins0    : number of continuum excitation energy bins
c wkb       : subroutine for WKB approximation for fission
c
      if (fislocal.eq.5) then
        nfisbar(Zix,Nix)=0
        fisfile=path(1:lenpath)//'fission/hfbpath/'//fischar
        inquire (file=fisfile,exist=lexist)
        if (lexist) then
          open (unit=2,status='old',file=fisfile)
  300     read(2,'(/11x,i4,12x,i4//)',end=330) ia,nbeta
          if (A.ne.ia) then
            do 310 i=1,nbeta
              read(2,'()')
  310       continue
            goto 300
          else
            do 320 i=1,nbeta
              read(2,'(f10.3,20x,f10.3)') bb,vv
              betafis(i)=betafiscor(Zix,Nix)*bb
              vfis(i)=vfiscor(Zix,Nix)*vv
  320       continue
            if (nbins0.eq.0) then
              nbinswkb=30
            else
              nbinswkb=nbins0
            endif
            call wkb(Z,A,Zix,Nix,nbar)
            nfisbar(Zix,Nix)=nbar
          endif
          goto 340
  330     fismodelx(Zix,Nix)=fismodelalt
          axtype(Zix,Nix,2)=2
          close (unit=2)
          goto 20
  340     close (unit=2)
        else
          fismodelx(Zix,Nix)=fismodelalt
          axtype(Zix,Nix,2)=2
          goto 20
        endif
      endif
c
c Read fission states
c
c modz,modn: help variables
c hbsfile  : file with head band transition states
c c2file   : file with class 2 states
c
      modz=mod(Z,2)
      modn=mod(N,2)
      if (modz.eq.0) then
        if (modn.eq.0) then
          hbsfile=path(1:lenpath)//'fission/states/hbstates.ee'
          c2file=path(1:lenpath)//'fission/states/class2states.ee'
        else
          hbsfile=path(1:lenpath)//'fission/states/hbstates.eo'
          c2file=path(1:lenpath)//'fission/states/class2states.eo'
        endif
      else
        if (modn.eq.0) then
          hbsfile=path(1:lenpath)//'fission/states/hbstates.oe'
          c2file=path(1:lenpath)//'fission/states/class2states.oe'
        else
          hbsfile=path(1:lenpath)//'fission/states/hbstates.oo'
          c2file=path(1:lenpath)//'fission/states/class2states.oo'
        endif
      endif
c
c Use user-defined files for head band and class 2 transition states
c
      if (hbtransfile(Zix,Nix)(1:1).ne.' ') hbsfile=hbtransfile(Zix,Nix)
      if (clas2file(Zix,Nix)(1:1).ne.' ')  c2file=clas2file(Zix,Nix)
c
c Read head band transition states
c
c flaghbstate: flag for head band states in fission
c nfistrhb   : number of head band transition states for barrier
c fecont     : start of continuum energy
c numlev     : maximum number of included discrete levels
c efistrhb   : energy of head band transition states
c jfistrhb   : spin of head band transition states
c pfistrhb   : parity of head band transition states
c
      if (flaghbstate) then
        open (unit=2,status='old',file=hbsfile)
        do 410 i=1,nfisbar(Zix,Nix)
          read(2,'(4x,i4,f8.3)',end=410,err=410)
     +      nfistrhb(Zix,Nix,i),fecont(Zix,Nix,i)
          if (nfistrhb(Zix,Nix,i).gt.numlev) then
            write(*,'(" TALYS-error: there are more than",i3,
     +        " head band states in file ",a73)') numlev,hbsfile
            write(*,'(" numlev in talys.cmb should be increased")')
            stop
          endif
          do 420 j=1,nfistrhb(Zix,Nix,i)
            read(2,'(4x,f11.6,f6.1,i5)',end=410,err=410)
     +        efistrhb(Zix,Nix,i,j),jfistrhb(Zix,Nix,i,j),
     +        pfistrhb(Zix,Nix,i,j)
 420      continue
 410    continue
        close (unit=2)
      endif
c
c Class2 states
c
c flagclass2: flag for class2 states in fission
c nclass2   : number of sets of class2 states
c nfisc2hb  : number of class2 states for barrier
c efisc2hb  : energy of class2 states
c jfisc2hb  : spin of class2 states
c pfisc2hb  : parity of class2 states
c
      if (flagclass2) then
        open (unit=2,status='old',file=c2file)
        nclass2(Zix,Nix)=nfisbar(Zix,Nix)-1
        do 430 i=1,nclass2(Zix,Nix)
          read(2,'(4x,i4)',end=430,err=430) nfisc2hb(Zix,Nix,i)
          if (nfisc2hb(Zix,Nix,i).gt.numlev) then
            write(*,'(" TALYS-error: there are more than",i3,
     +        " class 2 states in file ",a73)') numlev,c2file
            write(*,'(" numlev in talys.cmb should be increased")')
            stop
          endif
          do 440 j=1,nfisc2hb(Zix,Nix,i)
            read(2,'(4x,f11.6,f6.1,i5)',end=430,err=430)
     +        efisc2hb(Zix,Nix,i,j),jfisc2hb(Zix,Nix,i,j),
     +        pfisc2hb(Zix,Nix,i,j)
  440     continue
  430   continue
        close (unit=2)
      endif
c
c ************************* Default parameters *************************
c
c minertia    : moment of inertia of fission barrier deformation
c Rtransmom   : normalization constant for moment of inertia for
c               transition states
c Irigid      : rigid body value of moment of inertia
c minertc2    : moment of inertia for class2 states
c Rclass2mom  : normalization constant for moment of inertia for
c               class 2 states
c
      if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=1.
      if (fwidth(Zix,Nix,2).eq.0.) fwidth(Zix,Nix,2)=0.6
      if (nfisbar(Zix,Nix).eq.1.and.fbarrier(Zix,Nix,1).eq.0.) then
        fbarrier(Zix,Nix,1)=fbarrier(Zix,Nix,2)
        fwidth(Zix,Nix,1)=fwidth(Zix,Nix,2)
      endif
      do 510 i=1,numbar
        minertia(Zix,Nix,i)=Rtransmom(Zix,Nix,i)*Irigid(Zix,Nix,i)
        minertc2(Zix,Nix,i)=Rclass2mom(Zix,Nix,i)*Irigid(Zix,Nix,i)
  510 continue
c
c ********** Rotational bands on transition and class2 states **********
c
c rotband   : subroutine to build rotational bands on transition states
c rotclass2 : subroutine to build rotational bands on class2 states
c
      call rotband(Zix,Nix)
      if (flagclass2) call rotclass2(Zix,Nix)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
