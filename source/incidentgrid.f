      subroutine incidentgrid(fname,fexist)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : February 28, 2013
c | Task  : Predefined incident energy grids
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      fexist
      character*14 fname
      integer      nen,i,type,ktype,lenE,pos,iE1,iE2,Nnen
      real         E(numenin),Ein,degrid,Eeps,Emin,Emax,Erest
c
c ************************ Basic incident energy grid ******************
c
c eninc,Ein,Eeps: energies of basic energy grid in MeV
c degrid        : energy increment
c Nnen          : total number of energies
c
c The general incident energy grid we use is:
c
c        1.e-11,2.53e-8,  1.e-6,  1.e-5,  1.e-4 MeV
c         0.001,  0.002,  0.004,  0.007 MeV
c          0.01,   0.02,   0.04,   0.07 MeV
c         0.1-   1 MeV : dE= 0.1 MeV
c         1  -   8 MeV : dE= 0.2 MeV
c         8  -  15 MeV : dE= 0.5 MeV
c        15  -  30 MeV : dE= 1.0 MeV
c        30  -  60 MeV : dE= 2.0 MeV
c        60  -  80 MeV : dE= 5.0 MeV
c        80  - 160 MeV : dE=10.0 MeV
c        160 - 600 MeV : dE=20.0 MeV
c        600 -1000 MeV : dE=50.0 MeV
c
c This grid ensures that the excitation functions are sufficiently
c smooth at low energies, while at higher energies a somewhat coarser
c energy grid can be used.
c For various pre-defined energy grids (energy filenames), subsets of
c this grid can be taken.
c
      data (E(nen),nen=1,14) /
     +   1.e-11,2.53e-8,  1.e-6,  1.e-5,  1.e-4,  0.001,  0.002,  0.004,
     +    0.007,   0.01,   0.02,   0.04,   0.07,    0.1/
c
c Test for existence of pre-defined energy grid
c
      fexist=.false.
      do 10 type=0,6
        if (fname(1:1).eq.parsym(type)) then
          ktype=type
          goto 20
        endif
   10 continue
      return
   20 do 30 i=1,14
        if (fname(i:i).eq.'-') then
          pos=i
          goto 40
        endif
   30 continue
      return
   40 do 50 i=pos+2,10
        if (fname(i:i+4).eq.'.grid') then
          lenE=i-1
          goto 60
        endif
   50 continue
      return
   60 read(fname(2:pos-1),*,err=300,end=300) iE1
      read(fname(pos+1:lenE),*,err=300,end=300) iE2
      Emin=real(iE1)
      Emax=real(iE2)
c
c Set basic energy grid
c
      numinc=0
      nen=14
      Ein=0.1
      degrid=0.1
  110 Ein=Ein+degrid
      Eeps=Ein+1.e-4
      if (Ein.gt.Emaxtalys) goto 120
      if (nen.eq.numenin) goto 120
      nen=nen+1
      E(nen)=Ein
      if (Eeps.gt.1.) degrid=0.2
      if (Eeps.gt.8.) degrid=0.5
      if (Eeps.gt.15.) degrid=1.
      if (Eeps.gt.30.) degrid=2.
      if (Eeps.gt.60.) degrid=5.
      if (Eeps.gt.80.) degrid=10.
      if (Eeps.gt.160.) degrid=20.
      if (Eeps.gt.600.) degrid=50.
      goto 110
  120 Nnen=nen
c
c Fill energy grid
c
c enincmin: minimum incident energy
c enincmax: maximum incident energy
c numinc  : number of incident energies
c
      i=0
      do 210 nen=1,Nnen
        Ein=E(nen)
        Eeps=Ein+1.e-4
        if (Eeps.lt.Emin) goto 210
        if (ktype.ne.1) then
          if (Eeps.lt.1.) goto 210
          Erest=Eeps-real(int(Eeps))
          if (Erest.gt.0.1) goto 210
        endif
        if (Ein.le.Emax+1.e-4) then
          i=i+1
          eninc(i)=Ein
        else
          goto 220
        endif
  210 continue
  220 numinc=i
      enincmin=eninc(1)
      enincmax=eninc(numinc)
      if (numinc.gt.0) fexist=.true.
  300 return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
