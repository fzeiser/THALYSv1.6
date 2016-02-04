      subroutine gammapar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 8, 2013
c | Task  : Gamma ray parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*4  gamchar
      character*80 key
      character*90 gamfile
      integer      Zix,Nix,Z,A,N,ia,irad,l,nen,it
      real         eg1,sg1,gg1,eg2,sg2,gg2,egamref,enum,denom,ee,et,ft,
     +             factor,fe1(numTqrpa),fstrength,temp,dtemp
c
c ***************** Default giant resonance parameters *****************
c
c We use indices to indicate the type of radiation, the multipolarity
c and the number of the peak (which sometimes is two). They are
c (0(M) or 1(E), l, number), i.e. egr(0,1,1) means a constant for
c M1-radiation and egr(1,2,1) a constant for E2-radiation.
c
c 1. Read experimental E1 values from GR parameter file
c
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c ZZ,Z   : charge number of residual nucleus
c AA,A   : mass number of residual nucleus
c N      : neutron number of residual nucleus
c gamchar: help variable
c gamfile: giant resonance parameter file
c path   : directory containing structure files to be read
c lenpath: length of pathname
c ia     : mass number from GR table
c eg1,...: help variables
c egr    : energy of GR
c ggr    : width of GR
c sgr    : strength of GR
c
c GDR parameters from the table can always be overruled by a value
c given in the input file.
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      N=A-Z
      gamchar='z   '
      write(gamchar(2:4),'(i3.3)') Z
      gamfile=path(1:lenpath)//'gamma/gdr/'//gamchar
      inquire (file=gamfile,exist=lexist)
      if (.not.lexist) goto 20
      open (unit=2,status='old',file=gamfile)
   10 read(2,'(4x,i4,6f8.2)',end=20) ia,eg1,sg1,gg1,eg2,sg2,gg2
      if (A.ne.ia) goto 10
      if (egr(Zix,Nix,1,1,1).eq.0.) egr(Zix,Nix,1,1,1)=eg1
      if (sgr(Zix,Nix,1,1,1).eq.0.) sgr(Zix,Nix,1,1,1)=sg1
      if (ggr(Zix,Nix,1,1,1).eq.0.) ggr(Zix,Nix,1,1,1)=gg1
      if (egr(Zix,Nix,1,1,2).eq.0.) egr(Zix,Nix,1,1,2)=eg2
      if (sgr(Zix,Nix,1,1,2).eq.0.) sgr(Zix,Nix,1,1,2)=sg2
      if (ggr(Zix,Nix,1,1,2).eq.0.) ggr(Zix,Nix,1,1,2)=gg2
   20 close (unit=2)
c
c 1. Default GR parameterization compiled by Kopecky in RIPL
c    IAEA-TECDOC-1034, August 1998.
c
c onethird: 1/3
c pi      : pi
c kgr     : constant for gamma-ray strength function
c pi2h2c2 : 1/(pi*pi*clight*clight*hbar**2) in mb**-1.MeV**-2
c
c E1-radiation
c
      if (egr(Zix,Nix,1,1,1).eq.0.) egr(Zix,Nix,1,1,1)=
     +  31.2*A**(-onethird)+20.6*A**(-(1./6.))
      if (ggr(Zix,Nix,1,1,1).eq.0.) ggr(Zix,Nix,1,1,1)=
     +  0.026*(egr(Zix,Nix,1,1,1)**1.91)
      if (sgr(Zix,Nix,1,1,1).eq.0.) sgr(Zix,Nix,1,1,1)=
     +  1.2*120.*Z*N/(A*pi*ggr(Zix,Nix,1,1,1))
      kgr(Zix,Nix,1,1)=onethird*pi2h2c2
c
c E2-radiation
c
      if (egr(Zix,Nix,1,2,1).eq.0.) egr(Zix,Nix,1,2,1)=
     +  63.*A**(-onethird)
      if (ggr(Zix,Nix,1,2,1).eq.0.) ggr(Zix,Nix,1,2,1)=
     +  6.11-0.012*A
      if (sgr(Zix,Nix,1,2,1).eq.0.) sgr(Zix,Nix,1,2,1)=
     +  1.4e-4*(Z**2)*egr(Zix,Nix,1,2,1)/
     +  (A**onethird*ggr(Zix,Nix,1,2,1))
      kgr(Zix,Nix,1,2)=0.2*pi2h2c2
c
c E3-6 radiation
c
c gammax: number of l-values for gamma multipolarity
c
      do 110 l=3,gammax
        if (egr(Zix,Nix,1,l,1).eq.0.) egr(Zix,Nix,1,l,1)=
     +    egr(Zix,Nix,1,l-1,1)
        if (ggr(Zix,Nix,1,l,1).eq.0.) ggr(Zix,Nix,1,l,1)=
     +    ggr(Zix,Nix,1,l-1,1)
        if (sgr(Zix,Nix,1,l,1).eq.0.) sgr(Zix,Nix,1,l,1)=
     +    sgr(Zix,Nix,1,l-1,1)*8.e-4
        kgr(Zix,Nix,1,l)=pi2h2c2/(2*l+1.)
  110 continue
c
c Check whether number of giant resonances is two (in which case
c all parameters must be specified)
c
c irad    : variable to indicate M(=0) or E(=1) radiation
c ngr     : number of GR
c strength: model for E1 gamma-ray strength function
c
      do 130 irad=0,1
        do 130 l=1,gammax
          if (egr(Zix,Nix,irad,l,2).ne.0..and.
     +      sgr(Zix,Nix,irad,l,2).ne.0..and.
     +      egr(Zix,Nix,irad,l,2).ne.0.)  ngr(Zix,Nix,irad,l)=2
  130 continue
      if (strength.le.2.or.strength.eq.5) goto 300
c
c
c ***************** HFbcs or HFB QRPA strength functions ***************
c
c For Goriely's HFbcs or HFB QRPA strength function we overwrite the
c E1 strength function with tabulated results, if available.
c
c nTqrpa       : number of temperatures for QRPA
c numgamqrpa   : number of energies for QRPA strength function
c eqrpa        : energy grid for QRPA strength function
c fe1qrpa,fe1  : tabulated QRPA strength function
c adjust       : subroutine for energy-dependent parameter adjustment
c factor       : multiplication factor
c gamadjust    : logical for energy-dependent gamma adjustment
c etable,ftable: constant to adjust tabulated strength functions
c Tqrpa        : temperature for QRPA
c qrpaexist    : flag for existence of tabulated QRPA strength functions
c
      nTqrpa=1
      if (strength.eq.3)
     +  gamfile=path(1:lenpath)//'gamma/hfbcs/'//gamchar
      if (strength.eq.4)
     +  gamfile=path(1:lenpath)//'gamma/hfb/'//gamchar
      inquire (file=gamfile,exist=lexist)
      if (.not.lexist) goto 210
      open (unit=2,status='old',file=gamfile)
  220 read(2,'(10x,i4)',end=210) ia
      read(2,*)
      if (ia.ne.A) then
        do 230 nen=1,numgamqrpa+1
          read(2,*)
  230   continue
        goto 220
      endif
      do 240 nen=1,numgamqrpa
        read(2,'(f9.3,1p,20e12.3)') ee,(fe1(it),it=1,nTqrpa)
        if (gamadjust(Zix,Nix)) then
          key='etable'
          call adjust(ee,key,Zix,Nix,0,0,factor)
          et=etable(Zix,Nix)+factor-1.
          key='ftable'
          call adjust(ee,key,Zix,Nix,0,0,factor)
          ft=ftable(Zix,Nix)+factor-1.
        else
          et=etable(Zix,Nix)
          ft=ftable(Zix,Nix)
        endif
        eqrpa(Zix,Nix,nen)=ee+et
        do 250 it=1,nTqrpa
          fe1qrpa(Zix,Nix,nen,it)=onethird*pi2h2c2*fe1(it)*ft
  250   continue
  240 continue
      if (nTqrpa.gt.1) then
        dtemp=0.2
        temp=-dtemp
        do 260 it=1,nTqrpa
          temp=temp+dtemp
          Tqrpa(it)=temp
  260   continue
      endif
      qrpaexist(Zix,Nix)=.true.
  210 close (unit=2)
c
c ############################# M radiation ############################
c
c M1 radiation: strength function is related to that of E1
c
c egamref   : help variable
c enum,denom: help variables
c strengthM1: model for M1 gamma-ray strength function
c
  300 if (egr(Zix,Nix,0,1,1).eq.0.) egr(Zix,Nix,0,1,1)=
     +  41.*A**(-onethird)
      if (ggr(Zix,Nix,0,1,1).eq.0.) ggr(Zix,Nix,0,1,1)=4.
      kgr(Zix,Nix,0,1)=onethird*pi2h2c2
      if (sgr(Zix,Nix,0,1,1).eq.0.) then
        egamref=7.
        enum=kgr(Zix,Nix,0,1)*egamref*ggr(Zix,Nix,0,1,1)**2
        denom=(egamref**2-egr(Zix,Nix,0,1,1)**2)**2+
     +    (ggr(Zix,Nix,0,1,1)*egamref)**2
        if (strengthM1.eq.1) then
          factor=1.58e-9*A**0.47
        else
          factor=fstrength(Zix,Nix,0.,egamref,1,1)/(0.0588*A**0.878)
        endif
        sgr(Zix,Nix,0,1,1)=factor*denom/enum
      endif
c
c M2-6 radiation
c
      do 310 l=2,gammax
        if (egr(Zix,Nix,0,l,1).eq.0.) egr(Zix,Nix,0,l,1)=
     +    egr(Zix,Nix,0,l-1,1)
        if (ggr(Zix,Nix,0,l,1).eq.0.) ggr(Zix,Nix,0,l,1)=
     +    ggr(Zix,Nix,0,l-1,1)
        if (sgr(Zix,Nix,0,l,1).eq.0.) sgr(Zix,Nix,0,l,1)=
     +    sgr(Zix,Nix,0,l-1,1)*8.e-4
        kgr(Zix,Nix,0,l)=pi2h2c2/(2*l+1.)
  310 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
