      subroutine endfecis
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 19, 2013
c | Task  : ECIS calculation for incident particle on ENDF-6 energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      jlmloc,rotational,vibrational
      character*13 outfile
      integer      Zix,Nix,i1,i,ii,nen
      real         e
c
c ********************** Set ECIS input parameters *********************
c
c legendre        : logical for output of Legendre coefficients
c hint            : integration step size h
c rmatch          : matching radius
c projmass,parmass: mass of projectile
c spin,parspin    : spin of incident particle
c resmass,tarmass : mass of nucleus
c prodZ           : product of charges of projectile and target nucleus
c Ztarget         : charge number of target nucleus
c parZ            : charge number of particle
c k0              : index of incident particle
c Nband           : number of vibrational bands
c angbeg          : first angle
c anginc          : angle increment
c angend          : last angle
c flagendfecis    : flag for new ECIS calculation for ENDF-6 files
c Zindex          : charge number index for residual nucleus
c Nindex          : neutron number index for residual nucleus
c
c Specific ECIS flags:
c ecis2(9)=T  : output of total, reaction, elastic and inelastic c.s.
c ecis2(13)=F : output of transmission coefficients
c ecis2(14)=F : no output of elastic angular distribution
c
      legendre=.false.
      hint=0.
      rmatch=0.
      projmass=parmass(k0)
      spin=parspin(k0)
      resmass=tarmass
      prodZ=real(Ztarget*parZ(k0))
      Nband=0
      angbeg=0.
      anginc=180.
      angend=180.
c
c Loop over energies on ENDF-6 energy grid.
c
      if (flagendfecis)
     +  open (unit=9,status='unknown',file='ecisendf.inp')
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
c
c Standard ECIS inputs for phenomenological optical potentials
c
c ecis1,ecis2: 100 input flags ('T' or 'F') for ECIS
c jlmloc     : flag for JLM OMP
c colltype   : type of collectivity (D, V or R)
c flagrot    : flag for use of rotational optical model per
c              outgoing particle, if available
c rotational : flag for rotational input
c vibrational: flag for vibrational input
c title      : title of ECIS input file
c ncoll      : number of nuclear states
c npp        : number of optical potentials
c iterm      : number of iterations
c idvib      : identifier for existence of vibrational state inside
c              rotational model
c Elevel     : energy of level
c tarspin    : spin of target nucleus
c tarparity  : parity of target nucleus
c jlmexist   : flag for existence of tabulated radial matter density
c nrad       : number of radial points
c
c Some input flags for ECIS are energy dependent for the rotational
c model so ecis1 will be defined inside the energy loop.
c
      ecis1='FFFFFTFFFFFFFFFFFFFFFFFFTFFTFFFFFFFFFFFFFFFFFFFFFF'
      ecis2='FFFFFFFFTFFFFFFFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
c
c 1. Spherical nucleus
c
      jlmloc=.false.
      if (colltype(Zix,Nix).eq.'S'.or..not.flagrot(k0)) then
        rotational=.false.
        vibrational=.false.
        title='Spherical optical model                           '
        ncoll=1
        npp=1
        iterm=1
        idvib(1)=1
        Elevel(1)=0.
        tarspin=0.
        tarparity='+'
        if (jlmexist(Zix,Nix,k0)) then
          ecis1(7:7)='T'
          ecis1(15:15)='T'
          ecis1(29:29)='T'
          ecis1(41:41)='T'
          hint=0.1
          rmatch=18.
          nrad=182
          jlmloc=.true.
        endif
      else
c
c 2. Deformed nucleus
c
c jdis       : spin of level
c cparity    : parity (character)
c parlev     : parity of level
c ndef       : number of collective levels
c indexlevel : level index
c leveltype  : type of level (rotational (R) or vibrational (V))
c vibband    : band number of level
c maxband    : highest vibrational level added to rotational model
c edis       : energy of level
c Jlevel     : spin of level
c jdis       : spin of level
c Plevel     : parity of level
c iph,iphonon: phonon (1 or 2)
c Jband,lband: angular momentum
c iband      : band number of level
c Kmag,Kband : magnetic quantum number
c defpar     : deformation parameter
c vibbeta    : vibrational deformation parameter
c flagstate  : flag for optical model potential for each excited state
c
        iterm=0
        tarspin=jdis(Zix,Nix,0)
        tarparity=cparity(parlev(Zix,Nix,0))
        i1=0
        do 10 i=1,ndef(Zix,Nix)
          ii=indexlevel(Zix,Nix,i)
          if (leveltype(Zix,Nix,ii).ne.'V'.and.
     +      leveltype(Zix,Nix,ii).ne.'R') goto 10
          if (colltype(Zix,Nix).eq.'R'.and.
     +      vibband(Zix,Nix,i).gt.maxband) goto 10
          i1=i1+1
          idvib(i1)=vibband(Zix,Nix,i)
          Elevel(i1)=edis(Zix,Nix,ii)
          Jlevel(i1)=jdis(Zix,Nix,ii)
          Plevel(i1)=cparity(parlev(Zix,Nix,ii))
          iph(i1)=iphonon(Zix,Nix,i)
          if (iph(i1).eq.2) ecis1(2:2)='T'
          Jband(i1)=lband(Zix,Nix,i)
          iband(i1)=vibband(Zix,Nix,i)
          Kmag(i1)=Kband(Zix,Nix,i)
          vibbeta(i1)=defpar(Zix,Nix,i)
          Nband=max(Nband,iband(i1))
   10   continue
        ncoll=i1
        if (flagstate) then
          npp=ncoll
        else
          npp=1
        endif
        ecis1(12:12)='T'
c
c 2a. Vibrational model
c
        if (colltype(Zix,Nix).eq.'V') then
          rotational=.false.
          vibrational=.true.
          title='Vibrational optical model                         '
          do 20 i=1,ncoll
            idvib(i)=0
   20     continue
        else
c
c 2b. Rotational model
c
c Nrotbeta: number of deformation parameters for rotational nucleus
c nrot    : number of deformation parameters
c rotbeta : deformation parameters for rotational nucleus
c rotpar  : deformation parameters for rotational nucleus
c iqm     : largest order of deformation
c iqmax   : maximum l-value of multipole expansion
c
          rotational=.true.
          vibrational=.false.
          ecis1(1:1)='T'
          Nrotbeta=nrot(Zix,Nix)
          do 30 i=1,Nrotbeta
            rotbeta(i)=rotpar(Zix,Nix,i)
   30     continue
          if (colltype(Zix,Nix).eq.'R') then
            title='Symmetric rotational optical model                '
            iqm=2*Nrotbeta
          else
            title='Asymmetric rotational optical model               '
            ecis1(3:3)='T'
            iqm=2*(Nrotbeta-1)
          endif
          iqmax=8
        endif
      endif
c
c **************** ECIS input files for several energies ***************
c
c deftype  : deformation length (D) or parameter (B)
c flagrel  : flag for relativistic kinematics
c disp     : flag for dispersive optical model
c efer     : Fermi energy
c w2disp,..: constants for imaginary potentials
c flagecisinp: flag for existence of ecis input file
c nen6     : total number of energies
c e        : energy in MeV
c e6       : energies of ENDF-6 energy grid in MeV
c coullimit: energy limit for charged particle OMP calculation
c njmax    : maximal number of j-values in ECIS
c Atarget  : mass number of target nucleus
c onethird : 1/3
c numl     : maximum l-value (set in talys.cmb)
c
      if (deftype(Zix,Nix).eq.'B') ecis1(6:6)='F'
      if (flagrel) ecis1(8:8)='T'
      if (disp(Zix,Nix,k0)) then
        ecis1(10:10)='T'
        efer=ef(Zix,Nix,k0)
        w2disp=w2(Zix,Nix,k0)
        d3disp=d3(Zix,Nix,k0)
        d2disp=d2(Zix,Nix,k0)
      endif
      flagecisinp=.false.
      do 110 nen=1,nen6
        e=real(e6(nen))
        if (k0.gt.1.and.e.lt.coullimit(k0)) goto 110
c
c We use a simple formula to estimate the required number of j-values:
c    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
c and we always take a minimum of njmax=20.
c
        njmax=int(2.4*1.25*(real(Atarget)**onethird)*0.22*
     +    sqrt(projmass*e))
        njmax=max(njmax,20)
        njmax=min(njmax,numl-2)
        if (jlmloc) njmax=1600
c
c *************** Calculate optical potential parameters ***************
c
c optical: subroutine for determination of optical potential
c
        call optical(Zix,Nix,k0,e)
c
c ******************* Write ECIS input file ****************************
c
c flagcoulomb: flag for Coulomb excitation calculation with ECIS
c soswitch   : switch for deformed spin-orbit calculation and sequential
c              iterations in ECIS
c coulbar    : Coulomb barrier
c ecisinput  : subroutine to create ECIS input file
c
c For rotational nuclei, the switch at soswitch MeV needs to be made
c according to Pascal Romain.
c
        if (colltype(Zix,Nix).eq.'R'.and.flagrot(k0)) then
          if (k0.gt.1.and.flagcoulomb) ecis1(11:11)='T'
          if ((k0.eq.1.and.e.le.soswitch).or.
     +        (k0.gt.1.and.e.le.coulbar(k0))) then
            ecis1(13:13)='F'
            ecis1(21:21)='T'
            ecis1(42:42)='T'
          else
            ecis1(13:13)='T'
            ecis1(21:21)='F'
            ecis1(42:42)='F'
          endif
          if (k0.gt.1.and.e.le.0.05*coulbar(k0).and.
     +      e.le.2.*Elevel(ncoll)) e=0.1*Elevel(ncoll)
          if (flagrel) ecis1(8:8)='T'
        endif
        flagecisinp=.true.
        call ecisinput(Zix,Nix,k0,e,rotational,vibrational,jlmloc)
  110 continue
      if (.not.flagendfecis) return
      if (.not.flagecisinp) return
      write(9,'("fin")')
      close (unit=9)
c
c ************ ECIS calculation for outgoing energies ******************
c
c flagoutecis: flag for output of ECIS results
c outfile    : output file
c nulldev    : null device
c ecist      : subroutine ecis, adapted for TALYS
c ecisstatus : status of ECIS file
c
      if (flagoutecis) then
        outfile='ecisendf.out '
      else
        outfile=nulldev
      endif
      call ecist('ecisendf.inp ',outfile,'ecis.endfcs  ',
     +  'ecis.endfin  ','null         ','null         ','null         ')
      open (unit=9,status='unknown',file='ecisendf.inp')
      close (unit=9,status=ecisstatus)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
