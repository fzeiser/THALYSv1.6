      subroutine omppar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 15, 2012
c | Task  : Optical model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*1  omptype
      character*4  ompchar
      character*72 optmodfile
      character*90 ompfile
      integer      Zix,Nix,Z,N,A,k,ia,i,nomp,ii,nen
      real         e
c
c ************** Read optical model parameters from database ***********
c
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c ZZ,Z        : charge number of residual nucleus
c NN,N        : neutron number of residual nucleus
c AA,A        : mass number of residual nucleus
c omptype     : type of optical model (spherical or coupled)
c flaglocalomp: flag for local (y) or global (n) optical model
c ompchar     : help variable
c ompfile     : optical model parameter file
c optmodfileN : optical model parameter file for neutrons
c optmodfileP : optical model parameter file for protons
c path        : directory containing structure files to be read
c lenpath     : length of pathname
c ia          : mass number from level file
c nomp        : number of particles for optical model parameters
c ef          : Fermi energy
c rc0         : Coulomb radius
c rv0,av0     : real volume radius, diffuseness
c v1,v2,v3    : components for V
c w1,w2       : components for W
c rvd0,avd0   : real surface radius, diffuseness
c d1,d2,d3    : components for Wd
c rvso0,avso0 : real spin-orbit radius, diffuseness
c vso1,vso2   : components for Vso
c wso1,wso2   : components for Wso
c flagdisp    : flag for dispersive optical model
c flagjlm     : flag for using semi-microscopic JLM OMP
c disp        : flag for dispersive optical model
c colltype    : type of collectivity (D, V or R)
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      omptype=' '
      if (.not.flaglocalomp) goto 100
      ompchar='z   '
      write(ompchar(2:4),'(i3.3)') Z
      do 10 k=1,2
        if (k.eq.1) then
          if (optmodfileN(Zix)(1:1).ne.' ') then
            ompfile=optmodfileN(Zix)
          else
            ompfile=path(1:lenpath)//'optical/neutron/'//ompchar
          endif
        else
          if (optmodfileP(Zix)(1:1).ne.' ') then
            ompfile=optmodfileP(Zix)
          else
            ompfile=path(1:lenpath)//'optical/proton/'//ompchar
          endif
        endif
        omptype=' '
        inquire (file=ompfile,exist=lexist)
        if (.not.lexist) goto 10
        open (unit=2,status='old',file=ompfile)
   20   read(2,'(4x,2i4,3x,a1)',end=60) ia,nomp,omptype
        if (A.ne.ia) then
          do 30 i=1,nomp
            do 30 ii=1,4
              read(2,'()')
   30     continue
          goto 20
        endif
        omptype=' '
        do 40 i=1,nomp
          read(2,'(4x,f7.2,f8.3)') ef(Zix,Nix,k),rc0(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f9.6,f6.1,f7.1)') rv0(Zix,Nix,k),
     +      av0(Zix,Nix,k),v1(Zix,Nix,k),v2(Zix,Nix,k),v3(Zix,Nix,k),
     +      w1(Zix,Nix,k),w2(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f7.2)') rvd0(Zix,Nix,k),
     +      avd0(Zix,Nix,k),d1(Zix,Nix,k),d2(Zix,Nix,k),d3(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f6.1,f7.1)') rvso0(Zix,Nix,k),
     +      avso0(Zix,Nix,k),vso1(Zix,Nix,k),vso2(Zix,Nix,k),
     +      wso1(Zix,Nix,k),wso2(Zix,Nix,k)
          if (flagdisp) then
            if (nomp.eq.1.or.flagjlm) then
              disp(Zix,Nix,k)=.false.
            else
              disp(Zix,Nix,k)=.true.
            endif
          else
            goto 60
          endif
   40   continue
   60   close (unit=2)
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,k)=0.85*d1(Zix,Nix,k)
   10 continue
c
c ************************* Global optical model ***********************
c
c 1. Neutrons
c
c ompglobal: flag for use of global optical model
c onethird : 1/3
c twothird : 2/3
c
c Test if local OMP has been assigned.
c
  100 if (rv0(Zix,Nix,1).eq.0.) then
        if (flagdisp) disp(Zix,Nix,1)=.false.
        ompglobal(Zix,Nix,1)=.true.
        ef(Zix,Nix,1)=-11.2814+0.02646*A
        rv0(Zix,Nix,1)=1.3039-0.4054*A**(-onethird)
        av0(Zix,Nix,1)=0.6778-1.487e-4*A
        v1(Zix,Nix,1)=59.30-21.0*real(N-Z)/A-0.024*A
        v2(Zix,Nix,1)=7.228e-3-1.48e-6*A
        v3(Zix,Nix,1)=1.994e-5-2.0e-8*A
        w1(Zix,Nix,1)=12.195+0.0167*A
        w2(Zix,Nix,1)=73.55+0.0795*A
        rvd0(Zix,Nix,1)=1.3424-0.01585*A**onethird
        avd0(Zix,Nix,1)=0.5446-1.656e-4*A
        d1(Zix,Nix,1)=16.0-16.0*real(N-Z)/A
        d2(Zix,Nix,1)=0.0180+3.802e-3/(1.+exp((A-156.)/8.0))
        d3(Zix,Nix,1)=11.5
        vso1(Zix,Nix,1)=5.922+0.0030*A
        vso2(Zix,Nix,1)=0.0040
        rvso0(Zix,Nix,1)=1.1854-0.647*A**(-onethird)
        avso0(Zix,Nix,1)=0.59
        wso1(Zix,Nix,1)=-3.1
        wso2(Zix,Nix,1)=160.
        rc0(Zix,Nix,1)=0.
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,1)=0.85*d1(Zix,Nix,1)
      endif
c
c 2. Protons
c
      if (rv0(Zix,Nix,2).eq.0.) then
        if (flagdisp) disp(Zix,Nix,2)=.false.
        ompglobal(Zix,Nix,2)=.true.
        ef(Zix,Nix,2)=-8.4075+0.01378*A
        rv0(Zix,Nix,2)=1.3039-0.4054*A**(-onethird)
        av0(Zix,Nix,2)=0.6778-1.487e-4*A
        v1(Zix,Nix,2)=59.30+21.0*real(N-Z)/A-0.024*A
        v2(Zix,Nix,2)=7.067e-3+4.23e-6*A
        v3(Zix,Nix,2)=1.729e-5+1.136e-8*A
        w1(Zix,Nix,2)=14.667+0.009629*A
        w2(Zix,Nix,2)=73.55+0.0795*A
        rvd0(Zix,Nix,2)=1.3424-0.01585*A**onethird
        avd0(Zix,Nix,2)=0.5187+5.205e-4*A
        d1(Zix,Nix,2)=16.0+16.0*real(N-Z)/A
        d2(Zix,Nix,2)=0.0180+3.802e-3/(1.+exp((A-156.)/8.0))
        d3(Zix,Nix,2)=11.5
        vso1(Zix,Nix,2)=5.922+0.0030*A
        vso2(Zix,Nix,2)=0.0040
        rvso0(Zix,Nix,2)=1.1854-0.647*A**(-onethird)
        avso0(Zix,Nix,2)=0.59
        wso1(Zix,Nix,2)=-3.1
        wso2(Zix,Nix,2)=160.
        rc0(Zix,Nix,2)=1.198+0.697*A**(-twothird)+12.994*A**(-5./3.)
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,2)=0.85*d1(Zix,Nix,2)
      endif
c
c ************** Optical model file from user input file ***************
c
c numNph    : maximal number of neutrons away from the initial
c             compound nucleus for multiple pre-equilibrium emission
c numZph    : maximal number of protons away from the initial
c             compound nucleus for multiple pre-equilibrium emission
c optmod    : file with optical model parameters
c optmodfile: file with optical model parameters
c omplines  : number of lines in optical model file
c numomp    : maximum number of lines in optical model file
c eomp      : energies on optical model file
c vomp      : optical model parameters from file
c
      if (Zix.le.numZph.and.Nix.le.numNph) then
        do 210 k=1,6
          optmodfile=optmod(Zix,Nix,k)
          if (optmodfile(1:1).ne.' ') then
            open (unit=2,status='old',file=optmodfile)
            read(2,*) omplines(Zix,Nix,k)
            if (omplines(Zix,Nix,k).gt.numomp) then
              write(*,'(" TALYS-error: number of lines in OMP file > ",
     +          i4)') numomp
              stop
            endif
            eomp(Zix,Nix,k,0)=0.
            do 220 nen=1,omplines(Zix,Nix,k)
              read(2,*,err=300) eomp(Zix,Nix,k,nen),
     +          (vomp(Zix,Nix,k,nen,ii),ii=1,19)
  220       continue
            close (unit=2)
            if (omplines(Zix,Nix,k).eq.1) then
              omplines(Zix,Nix,k)=2
              eomp(Zix,Nix,k,2)=eomp(Zix,Nix,k,1)
              do 230 ii=1,19
                vomp(Zix,Nix,k,2,ii)=vomp(Zix,Nix,k,1,ii)
  230         continue
            endif
          endif
  210   continue
      endif
c
c Set OMP parameters for extension up to 1 GeV
c
c Zindex  : charge number index for residual nucleus
c Nindex  : neutron number index for residual nucleus
c Ejoin   : joining energy for high energy OMP
c enincmax: maximum incident energy
c optical : subroutine for determination of optical potential
c V0      : V at zero MeV
c Vjoin   : V at joining energy
c Wjoin   : W at joining energy
c
      do 310 k=1,2
        w3(Zix,Nix,k)=25.-0.0417*A
        w4(Zix,Nix,k)=250.
        if (Zix.eq.Zindex(0,0,k).and.Nix.eq.Nindex(0,0,k).and.
     +    enincmax.gt.Ejoin(k)) then
          e=0.
          call optical(Zix,Nix,k,e)
          V0(k)=v
          e=Ejoin(k)
          call optical(Zix,Nix,k,e)
          Vjoin(k)=v
          Wjoin(k)=w
        endif
  310 continue
      return
  300 write(*,'(" TALYS-error: Format error in ",a72)') optmodfile
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
