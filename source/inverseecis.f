      subroutine inverseecis(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 19, 2013
c | Task  : ECIS calculation for outgoing particles and energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      jlmloc,rotational,vibrational
      character*13 outfile
      integer      Zcomp,Ncomp,type,Zix,Nix,Z,A,i1,i,ii,nen
      real         e
c
c ********************** Set ECIS input parameters *********************
c
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c flaginvecis : logical for calculating inverse channel OMP
c legendre    : logical for output of Legendre coefficients
c hint        : integration step size h
c rmatch      : matching radius
c anginc      : angle increment
c angend      : last angle
c jlmloc      : flag for JLM OMP
c jlmexist    : flag for existence of tabulated radial matter density
c flagoutomp  : flag for output of optical model parameters
c numjlm      : maximum number of radial points
c rhojlmp     : density for protons
c rhojlmn     : density for neutrons
c flageciscalc: flag for new ECIS calculation for outgoing particles
c               and energy grid
c flagecisinp : flag for existence of ecis input file
c parskip     : logical to skip outgoing particle
c Zindex      : charge number index for residual nucleus
c Nindex      : neutron number index for residual nucleus
c ZZ,Z        : charge number of residual nucleus
c AA,A        : mass number of residual nucleus
c angbeg      : first angle
c nuc         : symbol of nucleus
c parname     : name of particle
c
c Specific ECIS flags:
c ecis2(9)=T  : output of total, reaction, elastic and inelastic c.s.
c ecis2(13)=T : output of transmission coefficients
c ecis2(14)=F : no output of elastic angular distribution
c
      flaginvecis=.true.
      legendre=.false.
      hint=0.
      rmatch=0.
      anginc=180.
      angend=180.
c
c Loop over all particle types and energies on standard energy grid.
c
      if (flagoutomp) then
        write(*,'(/" ######### OPTICAL MODEL PARAMETERS ##########")')
        if (jlmexist(Zcomp,Ncomp,1).or.jlmexist(Zcomp,Ncomp,6)) then
          write(*,'(/" Radial densities"/)')
          write(*,'(" Radius   Protons     Neutrons"/)')
          do 10 i=1,numjlm
            write(*,'(f7.3,1p,2e12.5)') 0.1*real(i),
     +        rhojlmp(Zcomp,Ncomp,i,1),rhojlmn(Zcomp,Ncomp,i,1)
   10     continue
        endif
      endif
      flagecisinp=.false.
      if (flageciscalc)
     +  open (unit=9,status='unknown',file='ecisinv.inp')
      do 110 type=1,6
        if (parskip(type)) goto 110
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        Z=ZZ(Zcomp,Ncomp,type)
        A=AA(Zcomp,Ncomp,type)
        if (type.eq.1) then
          angbeg=0.
        else
          angbeg=0.00001
        endif
        if (jlmexist(Zix,Nix,type).and.(colltype(Zix,Nix).eq.'S'.
     +    or..not.flagrot(type))) then
          jlmloc=.true.
        else
          jlmloc=.false.
        endif
c
c Output of optical model parameters, if requested.
c
        if (flagoutomp.and..not.jlmloc) then
          write(*,'(/11x,a8," on ",i3,a2/)') parname(type),A,nuc(Z)
          write(*,'("  Energy",4x,"V",5x,"rv",4x,"av",4x,"W",5x,
     +      "rw",4x,"aw",4x,"Vd",3x,"rvd",3x,"avd",4x,"Wd",
     +      3x,"rwd",3x,"awd",3x,"Vso",3x,"rvso",2x,"avso",
     +      2x,"Wso",3x,"rwso",2x,"awso",2x,"rc",/)')
        endif
c
c Standard ECIS inputs for phenomenological optical potentials
c
c ecis1,ecis2: 100 input flags ('T' or 'F') for ECIS
c Nband      : number of vibrational bands
c
c Some input flags for ECIS are energy dependent for the rotational
c model so ecis1 will be defined inside the energy loop.
c
        ecis1='FFFFFTFFFFFFFFFFFFFFFFFFTFFTFFFFFFFFFFFFFFFFFFFFFF'
        ecis2='FFFFFFFFTFFFTFFFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
        Nband=0
c
c 1. Spherical nucleus
c
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
c nrad       : number of radial points
c
        jlmloc=.false.
        if (colltype(Zix,Nix).eq.'S'.or..not.flagrot(type)) then
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
          if (jlmexist(Zix,Nix,type)) then
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
c maxband    : highest vibrational band added to rotational model
c edis       : energy of level
c Jlevel     : spin of level
c Plevel     : parity of level
c iph,iphonon: phonon (1 or 2)
c Kmag,Kband : magnetic quantum number
c iband      : band number of level
c Jband,lband: angular momentum
c vibbeta    : vibrational deformation parameter
c defpar     : deformation parameter
c flagstate  : flag for optical model potential for each excited state
c
          iterm=0
          tarspin=jdis(Zix,Nix,0)
          tarparity=cparity(parlev(Zix,Nix,0))
          i1=0
          do 120 i=1,ndef(Zix,Nix)
            ii=indexlevel(Zix,Nix,i)
            if (leveltype(Zix,Nix,ii).ne.'V'.and.
     +        leveltype(Zix,Nix,ii).ne.'R') goto 120
            if (colltype(Zix,Nix).eq.'R'.and.
     +        vibband(Zix,Nix,i).gt.maxband) goto 120
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
  120     continue
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
            do 130 i=1,ncoll
              idvib(i)=0
  130       continue
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
            do 140 i=1,Nrotbeta
              rotbeta(i)=rotpar(Zix,Nix,i)
  140       continue
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
c deftype         : deformation length (D) or parameter (B)
c flagrel         : flag for relativistic kinematics
c disp            : flag for dispersive optical model
c efer            : Fermi energy
c w2disp,.........: constants for imaginary potentials
c projmass,parmass: mass of projectile
c spin,parspin    : spin of incident particle
c nucmass,resmass : mass of nucleus
c prodZ           : product of charges of projectile and target nucleus
c parZ            : charge number of particle
c ebegin          : first energy point of energy grid
c eendmax         : last energy point of energy grid for maximum
c                   incident energy
c e               : energy in MeV
c egrid           : outgoing energy grid
c specmass        : specific mass
c
        if (deftype(Zix,Nix).eq.'B') ecis1(6:6)='F'
        if (flagrel) ecis1(8:8)='T'
        if (disp(Zix,Nix,type)) then
          ecis1(10:10)='T'
          efer=ef(Zix,Nix,type)
          w2disp=w2(Zix,Nix,type)
          d3disp=d3(Zix,Nix,type)
          d2disp=d2(Zix,Nix,type)
        endif
        projmass=parmass(type)
        spin=parspin(type)
        resmass=nucmass(Zix,Nix)
        prodZ=real(Z*parZ(type))
        do 210 nen=ebegin(type),eendmax(type)
          e=real(egrid(nen)/specmass(Zix,Nix,type))
c
c We use a simple formula to estimate the required number of j-values:
c    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
c and we always take a minimum of njmax=20.
c
c njmax   : maximal number of j-values in ECIS
c onethird: 1/3
c numl    : maximum l-value (set in talys.cmb)
c
          njmax=int(2.4*1.25*(real(A)**onethird)*0.22*
     +      sqrt(projmass*e))
          njmax=max(njmax,20)
          njmax=min(njmax,numl-2)
          if (jlmloc) njmax=1600
c
c *************** Calculate optical potential parameters ***************
c
c optical: subroutine for determination of optical potential
c v,rv,..: optical model parameters
c
          call optical(Zix,Nix,type,e)
          if (flagoutomp.and..not.jlmloc) then
            write(*,'(1x,f7.3,1x,6(f6.2,f6.3,f6.3),f6.3)')
     +        e,v,rv,av,w,rw,aw,vd,rvd,avd,wd,rwd,awd,vso,rvso,
     +        avso,wso,rwso,awso,rc
          endif
          if (.not.flageciscalc) goto 210
c
c ******************* Write ECIS input file ****************************
c
c soswitch : switch for deformed spin-orbit calculation and sequential
c            iterations in ECIS
c coulbar  : Coulomb barrier
c ecisinput: subroutine to create ECIS input file
c
c For rotational nuclei, the switch at soswitch MeV needs to be made
c according to Pascal Romain.
c
          if (colltype(Zix,Nix).eq.'R'.and.flagrot(type)) then
            if ((type.eq.1.and.e.le.soswitch).or.
     +        (type.gt.1.and.e.le.coulbar(type))) then
              ecis1(13:13)='F'
              ecis1(21:21)='T'
              ecis1(42:42)='T'
            else
              ecis1(13:13)='T'
              ecis1(21:21)='F'
              ecis1(42:42)='F'
            endif
            if (type.gt.1.and.e.le.0.05*coulbar(type).and.
     +        e.le.2.*Elevel(ncoll)) e=0.1*Elevel(ncoll)
            if (flagrel) ecis1(8:8)='T'
            if (disp(Zix,Nix,type)) ecis1(10:10)='T'
          endif
          flagecisinp=.true.
          call ecisinput(Zix,Nix,type,e,rotational,vibrational,jlmloc)
  210   continue
  110 continue
      flaginvecis=.false.
      if (.not.flageciscalc) return
      if (.not.flagecisinp) return
      write(9,'("fin")')
      close (unit=9)
c
c ************ ECIS calculation for outgoing energies ******************
c
c flagoutecis: flag for output of ECIS results
c outfile    : output file
c nulldev    : null device
c csfile     : file with inverse reaction cross sections
c ecist      : subroutine ecis, adapted for TALYS
c transfile  : file with transmission coefficients
c ecisstatus : status of ECIS file
c
      if (flagoutecis) then
        outfile='ecisinv.out  '
      else
        outfile=nulldev
      endif
      call ecist('ecisinv.inp  ',outfile,csfile,
     +    'ecis.invin   ',transfile,'null         ','null         ')
      open (unit=9,status='unknown',file='ecisinv.inp')
      close (unit=9,status=ecisstatus)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
