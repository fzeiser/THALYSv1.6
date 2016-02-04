      subroutine raynalcomp
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 27, 2013
c | Task  : ECIS calculation of compound cross sections (reference only)
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ilevel,type,Zix,Nix,Z,NLmax,nex,ildens
      real    ethrcm,econt
c
c In addition to a calculation by TALYS, a compound nucleus run by ECIS
c is performed. The results will however not be used for TALYS but are
c just for comparison.
c
c ********************** Set ECIS input parameters *********************
c
c Specific ECIS flags:
c ecis2(9)=T  : output of total, reaction, elastic and inelastic c.s.
c ecis2(14)=T : output of angular distributions
c ecis2(15)=T : output of Legendre coefficients
c ecis2(31)=T : compound
c ecis2(33)=T : No Engelbrecht-Weidenmuller transformation
c ecis2(34)=T : compound
c
c title      : title of ECIS input file
c ecis1,ecis2: 100 input flags ('T' or 'F') for ECIS
c flagrel    : flag for relativistic kinematics
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c disp       : flag for dispersive optical model
c k0         : index of incident particle
c
      open (unit=1,status='unknown',file='eciscomp.inp')
      title='Compound cross sections by ECIS                   '
      ecis1='FFTFFTFFFFFFFFTFFFFFFFFFFFFTFFFFFFFFFFFFFFFFFFFFFF'
      ecis2='FFFFFFFFTFFFFTTFTTTFTTTFTFFFFFTFTTFFFFFFFFFFTFFFFF'
      if (flagrel) ecis1(8:8)='T'
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      if (disp(Zix,Nix,k0)) ecis1(10:10)='T'
c
c Gamma emission
c
c parinclude: logical to include outgoing particle
c
      if (parinclude(0)) ecis2(36:36)='T'
c
c Width fluctuations
c
c flagwidth: flag for width fluctuation calculation
c bz1      : elastic enhancement factor
c
      if (flagwidth) then
        bz1=0.
      else
        bz1=1.
        ecis2(37:37)='F'
      endif
c
c ******************** Enumerate discrete levels ***********************
c
c ilevel          : level counter
c parskip         : logical to skip outgoing particle
c ZZ,Z            : charge number of residual nucleus
c ethrcm          : threshold energy
c eninccm         : center-of-mass incident energy in MeV
c Q               : Q-value for target nucleus
c NLmax,Nlast     : last discrete level
c typecomp        : particle type
c edis,elevelcomp : energy of level
c jdis,jcomp      : spin of level
c parlev,pcomp    : parity of level
c spincomp,parspin: spin of incident particle
c ejeccomp,parmass: mass of projectile
c masscomp,nucmass: mass of nucleus with indices (Z,N)
c prodZcomp       : product of charges
c parZ            : charge number of particle
c
      ilevel=0
      do 10 type=1,6
        if (parskip(type)) goto 10
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        Z=ZZ(0,0,type)
        ethrcm=eninccm+Q(type)
        NLmax=min(14,Nlast(Zix,Nix,0))
        do 20 nex=0,NLmax
          if (type.eq.k0.and.nex.eq.0) goto 20
          if (edis(Zix,Nix,nex).gt.ethrcm) goto 10
          ilevel=ilevel+1
          typecomp(ilevel)=type
          elevelcomp(ilevel)=edis(Zix,Nix,nex)
          jcomp(ilevel)=jdis(Zix,Nix,nex)
          if (parlev(Zix,Nix,nex).eq.1) then
            pcomp(ilevel)='+'
          else
            pcomp(ilevel)='-'
          endif
          spincomp(ilevel)=parspin(type)
          ejeccomp(ilevel)=parmass(type)
          masscomp(ilevel)=nucmass(Zix,Nix)
          prodZcomp(ilevel)=parZ(type)*Z
   20   continue
   10 continue
c
c ************** Add continua for photons and particles ****************
c
c 1. Photons
c
c aldcomp,alev: level density parameter with indices (Z,N)
c Umcomp      : matching point for U
c Exmatch     : matching point for Ex
c pair        : pairing energy
c tempcomp,T  : nuclear temperature
c E0comp,E0   : constant of temperature formula
c tgo         : slow s-wave neutron gamma width/spacing
c swaveth     : theoretical strength function for s-wave
c
      if (parinclude(0)) then
        Zix=0
        Nix=0
        aldcomp(0)=alev(Zix,Nix)
        Umcomp(0)=max(Exmatch(Zix,Nix,0)-pair(Zix,Nix),0.)
        tempcomp(0)=T(Zix,Nix,0)
        E0comp(0)=E0(Zix,Nix,0)
        Excomp(0)=Umcomp(0)+pair(Zix,Nix)
        tgo=swaveth(Zix,Nix)
      endif
c
c 2. Particles
c
c ildens : counter for continua
c econt  : first continuum energy
c deltaEx: excitation energy bin for population arrays
c Excomp : Matching Ex
c pair   : pairing energy
c
      ildens=0
      do 30 type=1,6
        if (parskip(type)) goto 30
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        ethrcm=eninccm+Q(type)
        econt=edis(Zix,Nix,NLmax)+0.5*deltaEx(Zix,Nix,NLmax+1)
        if (econt.lt.ethrcm) then
          ilevel=ilevel+1
          typecomp(ilevel)=type
          elevelcomp(ilevel)=econt
          jcomp(ilevel)=0.
          pcomp(ilevel)='+'
          spincomp(ilevel)=parspin(type)
          ejeccomp(ilevel)=parmass(type)
          masscomp(ilevel)=nucmass(Zix,Nix)
          prodZcomp(ilevel)=parZ(type)*Z
          ildens=ildens+1
          aldcomp(ildens)=alev(Zix,Nix)
          Umcomp(ildens)=max(Exmatch(Zix,Nix,0)-pair(Zix,Nix),0.)
          tempcomp(ildens)=T(Zix,Nix,0)
          E0comp(ildens)=E0(Zix,Nix,0)
          Excomp(ildens)=Umcomp(ildens)+pair(Zix,Nix)
        endif
   30 continue
c
c ************************ ECIS parameters *****************************
c
c ncoll : number of nuclear states
c iterm : number of iterations
c npp   : number of optical potentials
c hint  : integration step size h
c rmatch: matching radius
c nsp1  : number of uncoupled states and continua
c nsp2  : number of uncoupled states with angular distribution
c ncont : number of continua
c
      elevelcomp(0)=0.
      typecomp(0)=k0
      ncoll=1
      iterm=1
      npp=ilevel+1
      hint=0.
      rmatch=0.
      nsp1=ilevel
      nsp2=ilevel-ildens
      ncont=ildens
c
c We use a simple formula to estimate the required number of j-values:
c    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
c and we always take a minimum of njmax=20.
c
c projmass       : mass of projectile
c njmax          : maximal number of j-values in ECIS
c Atarget        : mass number of target nucleus
c onethird       : 1/3
c Einc           : incident energy in MeV
c numl           : maximum l-value (set in talys.cmb)
c parN           : neutron number of particle
c cparity        : parity of target (character)
c tarparity      : parity of target
c resmass,tarmass: mass of target nucleus
c prodZ          : product of charges of projectile and target nucleus
c Ztarget        : charge number of target nucleus
c angbeg         : first angle
c anginc         : angle increment
c nangle         : number of angles
c angend         : last angle
c
      projmass=parmass(k0)
      njmax=int(2.4*1.25*(real(Atarget)**onethird)*0.22*
     +  sqrt(projmass*Einc))
      njmax=max(njmax,20)
      njmax=min(njmax,numl)
      tarparity=cparity(parlev(parZ(k0),parN(k0),0))
      spin=parspin(k0)
      resmass=tarmass
      prodZ=real(Ztarget*parZ(k0))
      if (k0.eq.1) then
        angbeg=0.
      else
        angbeg=0.00001
      endif
      anginc=180./nangle
      angend=180.
c
c ******************* Write ECIS input file ****************************
c
c eciscompound: subroutine to create ECIS input file for compound cross
c               section
c
      call eciscompound
      write(1,'("fin")')
      close (unit=1)
c
c ************************** ECIS calculation **************************
c
c flagoutecis: flag for output of ECIS results
c ecist      : subroutine ecis, adapted for TALYS
c nulldev    : null device
c
      if (flagoutecis) then
        call ecist('eciscomp.inp ','eciscomp.out ',
     +    'ecis.comcs   ','ecis.comin   ','null         ',
     +    'ecis.comang  ','ecis.comleg  ')
      else
        call ecist('eciscomp.inp ',nulldev,
     +    'ecis.comcs   ','ecis.comin   ','null         ',
     +    'ecis.comang  ','ecis.comleg  ')
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
