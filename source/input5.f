      subroutine input5
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 11, 2013
c | Task  : Read input for fifth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      flagassign
      character*1  ch
      character*80 word(40),key,value,cval
      integer      Zix,Nix,ibar,irad,Z,N,oddZ,oddN,lval,type,mt,is,igr,
     +             i,class,iz,ia,ival,type2,omptype,nr,k,l,iword,ilev0,
     +             nbr,ilev1
      real         alphaldall,betaldall,gammashell1all,br,sum,val,
     +             Pshiftconstantall
c
c ************** Defaults for fifth set of input variables *************
c
c k0            : index for incident particle
c flagendf      : flag for information for ENDF-6 file
c eninclow      : minimal incident energy for nuclear model calculations
c Atarget       : mass number of target nucleus
c strength      : model for E1 gamma-ray strength function
c gnorm         : gamma normalization factor
c Rspincut      : adjustable constant (global) for spin cutoff factor
c spincutmodel  : model for spin cutoff factor for ground state
c shellmodel    : model for shell correction energies
c kvibmodel     : model for vibrational enhancement
c gammashell2   : gamma-constant for asymptotic level density parameter
c flagcol       : flag for collective enhancement of level density
c flagcolldamp  : flag for damping of collective effects in effective
c                 level density (without explicit collective
c                 enhancement). Only used for Bruyeres-le-Chatel
c                 (Pascal Romain) fission model
c pairconstant  : constant for pairing energy systematics
c Ufermi        : energy of Fermi distribution for damping of
c                 ground-state rotational effects
c cfermi        : width of Fermi distribution for damping of
c                 ground-state rotational effects
c Ufermibf      : energy of Fermi distribution for damping of barrier
c                 rotational effects
c cfermibf      : width of Fermi distribution for damping of barrier
c                 rotational effects
c Kph           : constant for single-particle level density parameter
c                 (g=A/Kph)
c M2constant    : overall constant for matrix element in exciton model
c M2limit       : constant for asymptotic value for matrix element
c M2shift       : constant for energy shift for matrix element
c Rpinu,Rnupi.  : ratio for two-component matrix element
c Esurf0        : well depth for surface interaction
c Rgamma        : adjustable parameter for pre-equilibrium gamma decay
c elwidth       : width of elastic peak in MeV
c xscaptherm    : thermal capture cross section
c xsptherm      : thermal (n,p) cross section
c xsalphatherm  : thermal (n,a) cross section
c Zix           : charge number index for residual nucleus
c numZ          : maximal number of protons away from the initial
c                 compound nucleus
c Nix           : neutron number index for residual nucleus
c numN          : maximal number of neutrons away from the initial
c                 compound nucleus
c ldmodel       : level density model
c alphad        : alpha-constant for asymptotic level density parameter
c betald        : beta-constant for asymptotic level density parameter
c gammashell1   : gamma-constant for asymptotic level density parameter
c Pshiftconstant: global constant for pairing shift
c alev          : level density parameter
c alimit        : asymptotic level density parameter
c gammald       : gamma-constant for asymptotic level density parameter
c pair          : total pairing correction
c ibar          : fission barrier
c numbar        : number of fission barriers
c Pshift        : adjustable pairing shift
c deltaW        : shell correction in nuclear mass
c Exmatch       : matching point for Ex
c T             : nuclear temperature
c E0            : constant of temperature formula
c Nlow          : lowest discrete level for temperature matching
c Ntop          : highest discrete level for temperature matching
c Krotconstant  : normalization constant for rotational enhancement
c s2adjust      : adjustable constant (Z,A,barrier-dependent) for spin
c                 cutoff parameter
c ctable,ptable : constant to adjust tabulated level densities
c cglobal       : global constant to adjust tabulated level densities
c pglobal       : global constant to adjust tabulated level densities
c g             : single-particle level density parameter
c gp            : single-particle proton level density parameter
c gn            : single-particle neutron level density parameter
c gamgam        : total radiative width in eV
c D0            : experimental s-wave resonance spacing in eV
c etable,ftable : constant to adjust tabulated strength functions
c ldadjust      : logical for energy-dependent level density adjustment
c gamadjust     : logical for energy-dependent gamma adjustment
c fisadjust     : logical for energy-dependent fission adjustment
c irad          : variable to indicate M(=0) or E(=1) radiation
c lval          : multipolarity
c numgam        : maximum number of l-values for gamma multipolarity
c igr           : giant resonance
c egr           : energy of GR
c ggr           : width of GR
c sgr           : strength of GR
c epr           : energy of PR
c gpr           : width of PR
c tpr           : strength of PR
c egradjust.....: adjustable factors for giant resonance parameters
c                 (default 1.)
c fiso          : correction factor for isospin forbidden transitions
c aadjust,...   : adjustable factors for level density parameters
c                 (default 1.)
c axtype        : type of axiality of barrier
c                 1: axial symmetry
c                 2: left-right asymmetry
c                 3: triaxial and left-right symmetry
c                 4: triaxial no left-right symmetry
c                 5: no symmetry
c fbarrier      : height of fission barrier
c fwidth        : width of fission barrier
c fbaradjust,...: adjustable factors for fission parameters
c                 (default 1.)
c betafiscor    : adjustable factor for fission path width
c Zinit         : charge number of initial compound nucleus
c Ninit         : neutron number of initial compound nucleus
c vfiscor       : adjustable factor for fission path height
c fismodel      : fission model
c Rtransmom     : normalization constant for moment of inertia for
c                 transition states
c Rclass2mom    : normalization constant for moment of inertia for
c                 class 2 states
c widthc2       : width of class2 states
c Ninit         : neutron number of initial compound nucleus
c hbtransfile   : file with head band transition states
c class2file    : file with class 2 transition states
c levelfile     : discrete level file
c deformfile    : deformation parameter file
c optmodfileN   : optical model parameter file for neutrons
c optmodfileP   : optical model parameter file for protons
c radialfile    : radial matter density file
c ompenergyfile : file with energies for OMP calculation (ENDF files
c                 only)
c radialmodel   : model for radial matter densities (JLM OMP only)
c massnucleus   : mass of nucleus in amu as read from user input file
c massexcess    : mass excess in MeV as read from user input file
c beta2         : deformation parameter
c msdbins       : number of energy points for DWBA calculation for MSD
c Emsdmin       : minimal outgoing energy for MSD calculation
c Cstrip        : adjustable parameter for stripping/pick-up reactions
c Cknock        : adjustable parameter for knockout reactions
c Cbreak        : adjustable parameter for breakup reactions
c nbranch       : number of branching levels
c branchlevel   : level to which branching takes place
c branchratio   : gamma-ray branching ratio to level
c v1adjust....  : adjustable factors for OMP (default 1.)
c preeqadjust   : logical for energy-dependent pre-eq. adjustment
c Nadjust       : number of adjustable parameters
c adjustkey     : keyword for local adjustment
c adjustfile    : file for local adjustment
c adjustpar     : local adjustment parameters
c adjustix      : local adjustment index
c nenadjust     : number of tabulated energies of local adjustment
c Eadjust       : tabulated energy of local adjustment
c Dadjust       : tabulated depth of local adjustment
c ompadjustF    : logical for local OMP adjustment
c ompadjustN    : number of energy ranges for local OMP adjustment
c ompadjustE1   : start energy of local OMP adjustment
c ompadjustE2   : end energy of local OMP adjustment
c ompadjustD    : depth of local OMP adjustment
c ompadjusts    : variance of local OMP adjustment
c ompadjustp    : logical for energy-dependent OMP adjustment
c Ejoin         : joining energy for high energy OMP
c Vinfadjust    : adjustable factor for high energy limit of
c                 real central potential
c tljadjust     : logical for energy-dependent Tlj adjustment
c tladjust      : adjustable factor for Tlj (default 1.)
c jlmmode       : option for JLM imaginary potential normalization
c flagrescue    : flag for final rescue: normalization to data
c rescuefile    : file with incident energy dependent adjustment factors
c grescue       : global multiplication factor for incident energy
c                 dependent adjustment factors
c flaglabddx    : flag for calculation of DDX in LAB system
c aradialcor    : adjustable parameter for shape of DF alpha potential
c adepthcor     : adjustable parameter for depth of DF alpha potential
c nanglerec     : number of recoil angles
c RprimeU       : potential scattering radius
c nTmax         : effective number of temperatures for Maxwellian
c                 average
c astroT9       : temperature, in 10^9 K, for Maxwellian average
c astroE        : energy, in MeV, for Maxwellian average
c Ebeam         : incident energy in MeV for isotope production
c Eback         : lower end of energy range in MeV for isotope
c                 production
c radiounit     : unit for radioactivity: Bq, kBq, MBq, Gbq,
c                 mCi, Ci or kCi
c yieldunit     : unit for isotope yield: num (number),
c                 mug (micro-gram), mg, g, or kg
c Ibeam         : beam current in mA for isotope production
c Area          : target area in cm^2 for isotope production
c Tirrad        : irradiation time per unit
c unitTirrad    : irradiation time unit (y,d,h,m,s)
c Tcool         : cooling time per unit
c unitTcool     : cooling time unit (y,d,h,m,s)
c rhotarget     : target material density
c
      if (k0.eq.1.and.flagendf) then
        eninclow=0.
      else
        eninclow=1.e-6
      endif
c
c Advice of S. Goriely: no gamma normalization for A < 40.
c
      if (k0.ne.1.or.Atarget.lt.40.or.strength.eq.3.or.strength.eq.4)
     +  then
        gnorm=1.
      else
        gnorm=-1.
      endif
      Rspincut=1.
      spincutmodel=1
      shellmodel=1
      kvibmodel=2
      gammashell2=0.
      pairconstant=12.
      Ufermi=30.
      cfermi=5.
      Ufermibf=45.
      cfermibf=5.
      Kph=15.
      M2constant=1.
      M2limit=1.
      M2shift=1.
      Rpipi=1.
      Rnunu=1.5
      Rpinu=1.
      Rnupi=1.
      Esurf0=-1.
      Rgamma=2.
      elwidth=0.5
      xscaptherm=0.
      xsptherm=0.
      xsalphatherm=0.
      alphaldall=-99.
      betaldall=-99.
      gammashell1all=-99.
      Pshiftconstantall=-99.
      do 10 Zix=0,numZ
        do 20 Nix=0,numN
          if (ldmodel(Zix,Nix).eq.1.or.ldmodel(Zix,Nix).ge.4) then
            if (flagcol(Zix,Nix)) then
              alphald(Zix,Nix)=0.0207305
              betald(Zix,Nix)=0.229537
              gammashell1(Zix,Nix)=0.473625
              Pshiftconstant(Zix,Nix)=0.
            else
              alphald(Zix,Nix)=0.0692559
              betald(Zix,Nix)=0.282769
              gammashell1(Zix,Nix)=0.433090
              Pshiftconstant(Zix,Nix)=0.
            endif
            if (flagcolldamp) then
              alphald(Zix,Nix)=0.0666
              betald(Zix,Nix)=0.258
              gammashell1(Zix,Nix)=0.459
              Pshiftconstant(Zix,Nix)=0.
            endif
          endif
          if (ldmodel(Zix,Nix).eq.2) then
            if (flagcol(Zix,Nix)) then
              alphald(Zix,Nix)=0.0381563
              betald(Zix,Nix)=0.105378
              gammashell1(Zix,Nix)=0.546474
              Pshiftconstant(Zix,Nix)=0.743229
            else
              alphald(Zix,Nix)=0.0722396
              betald(Zix,Nix)=0.195267
              gammashell1(Zix,Nix)=0.410289
              Pshiftconstant(Zix,Nix)=0.173015
            endif
          endif
          if (ldmodel(Zix,Nix).eq.3) then
            if (flagcol(Zix,Nix)) then
              alphald(Zix,Nix)=0.0357750
              betald(Zix,Nix)=0.135307
              gammashell1(Zix,Nix)=0.699663
              Pshiftconstant(Zix,Nix)=-0.149106
            else
              alphald(Zix,Nix)=0.110575
              betald(Zix,Nix)=0.0313662
              gammashell1(Zix,Nix)=0.648723
              Pshiftconstant(Zix,Nix)=1.13208
            endif
          endif
          alev(Zix,Nix)=0.
          alimit(Zix,Nix)=0.
          gammald(Zix,Nix)=-1.
          pair(Zix,Nix)=1.e-20
          do 30 ibar=0,numbar
            Pshift(Zix,Nix,ibar)=1.e-20
            deltaW(Zix,Nix,ibar)=0.
            Exmatch(Zix,Nix,ibar)=0.
            T(Zix,Nix,ibar)=0.
            E0(Zix,Nix,ibar)=1.e-20
            Tadjust(Zix,Nix,ibar)=1.
            E0adjust(Zix,Nix,ibar)=1.
            Exmatchadjust(Zix,Nix,ibar)=1.
            Nlow(Zix,Nix,ibar)=-1
            Ntop(Zix,Nix,ibar)=-1
            s2adjust(Zix,Nix,ibar)=1.
            Krotconstant(Zix,Nix,ibar)=1.
            ctable(Zix,Nix,ibar)=cglobal
            ptable(Zix,Nix,ibar)=pglobal
   30     continue
          g(Zix,Nix)=0.
          gp(Zix,Nix)=0.
          gn(Zix,Nix)=0.
          gamgam(Zix,Nix)=0.
          D0(Zix,Nix)=0.
          etable(Zix,Nix)=0.
          ftable(Zix,Nix)=1.
          ldadjust(Zix,Nix)=.false.
          gamadjust(Zix,Nix)=.false.
          fisadjust(Zix,Nix)=.false.
          do 40 irad=0,1
            do 40 lval=1,numgam
              do 40 igr=1,2
                egr(Zix,Nix,irad,lval,igr)=0.
                ggr(Zix,Nix,irad,lval,igr)=0.
                sgr(Zix,Nix,irad,lval,igr)=0.
                epr(Zix,Nix,irad,lval,igr)=0.
                gpr(Zix,Nix,irad,lval,igr)=0.
                tpr(Zix,Nix,irad,lval,igr)=0.
                egradjust(Zix,Nix,irad,lval,igr)=1.
                ggradjust(Zix,Nix,irad,lval,igr)=1.
                sgradjust(Zix,Nix,irad,lval,igr)=1.
                epradjust(Zix,Nix,irad,lval,igr)=1.
                gpradjust(Zix,Nix,irad,lval,igr)=1.
                tpradjust(Zix,Nix,irad,lval,igr)=1.
   40     continue
          do type=0,6
            fiso(Zix,Nix,type)=1.
          enddo
          if (Zix.eq.0.and.Nix.eq.0) then
            if (Zinit.eq.Ninit) then
              if (k0.eq.1) fiso(Zix,Nix,k0)=2.
              if (k0.eq.2) fiso(Zix,Nix,k0)=2.
              if (k0.eq.6) fiso(Zix,Nix,k0)=5.
            endif
            if (Zinit.eq.Ninit-1.or.Zinit.eq.Ninit+1) then
              if (k0.eq.1) fiso(Zix,Nix,k0)=1.5
              if (k0.eq.2) fiso(Zix,Nix,k0)=1.5
              if (k0.eq.6) fiso(Zix,Nix,k0)=1.5
            endif
          endif
          aadjust(Zix,Nix)=1.
          gnadjust(Zix,Nix)=1.
          gpadjust(Zix,Nix)=1.
          gamgamadjust(Zix,Nix)=1.
          do 50 ibar=1,numbar
            axtype(Zix,Nix,ibar)=1
            fbarrier(Zix,Nix,ibar)=0.
            fwidth(Zix,Nix,ibar)=0.
            fbaradjust(Zix,Nix,ibar)=1.
            fwidthadjust(Zix,Nix,ibar)=1.
            Rtransmom(Zix,Nix,ibar)=1.
            Rclass2mom(Zix,Nix,ibar)=1.
            widthc2(Zix,Nix,ibar)=0.2
   50     continue
          betafiscor(Zix,Nix)=1.
          Z=Zinit-Zix
          N=Ninit-Nix
          oddZ=mod(Z,2)
          oddN=mod(N,2)
          vfiscor(Zix,Nix)=1.
          if (oddZ.eq.0.and.oddN.eq.0) vfiscor(Zix,Nix)=0.86
          if (oddZ.eq.1.and.oddN.eq.0) vfiscor(Zix,Nix)=0.94
          if (oddZ.eq.0.and.oddN.eq.1) vfiscor(Zix,Nix)=0.89
          if (oddZ.eq.1.and.oddN.eq.1) vfiscor(Zix,Nix)=1.02
          fismodelx(Zix,Nix)=fismodel
          if (Ninit-Nix.gt.144.or.fismodel.eq.5) axtype(Zix,Nix,1)=3
          if (fismodel.ne.5) axtype(Zix,Nix,2)=2
          Rtransmom(Zix,Nix,1)=0.6
          hbtransfile(Zix,Nix)='                                       '
          clas2file(Zix,Nix)='                                         '
   20   continue
        levelfile(Zix)='                                               '
        deformfile(Zix)='                                              '
        optmodfileN(Zix)='                                            '
        optmodfileP(Zix)='                                            '
        radialfile(Zix)='                                              '
   10 continue
      ompenergyfile='                                            '
      radialmodel=2
      do 60 Zix=0,numZ+4
        do 60 Nix=0,numN+4
          massnucleus(Zix,Nix)=0.
          massexcess(Zix,Nix)=0.
          beta2(Zix,Nix,0)=0.
          beta2(Zix,Nix,1)=0.6
          beta2(Zix,Nix,2)=0.8
          beta2(Zix,Nix,3)=1.
   60 continue
      msdbins=6
      Emsdmin=0.
      do 70 type=0,6
        Cstrip(type)=1.
        Cknock(type)=1.
c
c Extra setting for problems with Kalbach model for (alpha,p) reactions
c
        if (k0.eq.6) then
          Cbreak(type)=0.
        else
          Cbreak(type)=1.
        endif
   70 continue
      do 72 k=0,numlev
        do 72 i=0,numlev
          do 72 Nix=0,numN
            do 72 Zix=0,numZ
              branchlevel(Zix,Nix,i,k)=0
              branchratio(Zix,Nix,i,k)=0.
   72 continue
      do 74 i=0,numlev
        do 74 Nix=0,numN
          do 74 Zix=0,numZ
            nbranch(Zix,Nix,i)=0
   74 continue
      do 80 type=1,6
        v1adjust(type)=1.
        v2adjust(type)=1.
        v3adjust(type)=1.
        v4adjust(type)=1.
        rvadjust(type)=1.
        avadjust(type)=1.
        w1adjust(type)=1.
        w2adjust(type)=1.
        w3adjust(type)=1.
        w4adjust(type)=1.
        rwadjust(type)=-1.
        awadjust(type)=-1.
        rvdadjust(type)=-1.
        avdadjust(type)=-1.
        d1adjust(type)=1.
        d2adjust(type)=1.
        d3adjust(type)=1.
        rwdadjust(type)=1.
        awdadjust(type)=1.
        vso1adjust(type)=1.
        vso2adjust(type)=1.
        rvsoadjust(type)=1.
        avsoadjust(type)=1.
        wso1adjust(type)=1.
        wso2adjust(type)=1.
        rwsoadjust(type)=-1.
        awsoadjust(type)=-1.
        rcadjust(type)=1.
        Vinfadjust(type)=1.
        do 90 omptype=1,numompadj
          ompadjustN(type,omptype)=0
          do 92 nr=1,numrange
            ompadjustE1(type,omptype,nr)=0.
            ompadjustE2(type,omptype,nr)=0.
            ompadjustD(type,omptype,nr)=1.
            ompadjusts(type,omptype,nr)=1.
   92     continue
   90   continue
        ompadjustF(type)=.false.
        ompadjustp(type)=.false.
   80 continue
      do 81 type=1,6
        Ejoin(type)=200.
   81 continue
      preeqadjust=.false.
      Nadjust=0
      do 82 i=1,numadj
        adjustkey(i)='                                                 '
        adjustfile(i)='                                                '
        nenadjust(i)=0
        do 84 k=1,4
          adjustix(i,k)=0
          adjustpar(i,k)=0.
   84   continue
        do 85 k=1,numenadj
          Eadjust(i,k)=0.
          Dadjust(i,k)=0.
   85   continue
   82 continue
      do 94 type=-1,6
        tljadjust(type)=.false.
        do 96 l=0,numl
          tladjust(type,l)=1.
   96   continue
   94 continue
      flagrescue=.false.
      do 100 mt=1,nummt
        do 105 is=-1,1
          rescuefile(mt,is)='                                          '
          grescue(mt,is)=1.
  105   continue
  100 continue
      jlmmode=0
      lvadjust=1.
      lwadjust=1.
      lv1adjust=1.
      lw1adjust=1.
      lvsoadjust=1.
      lwsoadjust=1.
      aradialcor=1.
      adepthcor=1.
      if (flaglabddx) then
        nanglerec=numangrec
      else
        nanglerec=1
      endif
      RprimeU=0.
      nTmax=numT
      astroT9=0.
      astroE=0.
      Ebeam=-1.
      Eback=-1.
      Ibeam=1.
      radiounit='gbq'
      yieldunit='num'
      Area=1.
      Tirrad(1)=1
      unitTirrad(1)='d'
      Tcool(1)=1
      unitTcool(1)='d'
      do 108 k=2,5
        Tirrad(k)=0
        unitTirrad(k)=' '
        Tcool(k)=0
        unitTcool(k)=' '
  108 continue
      rhotarget=-1.
c
c **************** Read fifth set of input variables *******************
c
c nlines     : number of input lines
c getkeywords: subroutine to retrieve keywords and values from input
c              line
c inline     : input line
c word       : words on input line
c key        : keyword
c value      : value or string
c ch         : character
c
c The keyword is identified and the corresponding values are read.
c Erroneous input is immediately checked. The keywords and number of
c values on each line are retrieved from the input.
c
      do 110 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
        Zix=0
        Nix=0
        type=0
        lval=0
        ibar=0
        igr=1
c
c Test for keywords
c
c Here, the various model parameters can be set to overrule the default
c values. Most default values will be computed later on, since they
c require more computation (e.g. level density parameters).
c
c Each keyword is characterized by a certain order of parameter
c and value input. They are distinguished by different classes.
c
c Classes:
c
c 1: keyword Z A real-value [optional: local adjustment]
c 2: keyword Z A integer-value
c 3: keyword Z A real-value barrier [optional: local adjustment]
c 4: keyword Z A integer-value barrier [optional: local adjustment]
c 5: keyword Z A real-value rad-type l-val [optional: local adjustment]
c 6: keyword particle-type real-value [optional: local adjustment]
c 7: keyword particle-type integer-value
c 8: keyword value [optional: local adjustment]
c 9: keyword Z filename
c
c class     : input class
c getvalues : subroutine to assign values to keywords
c iz        : charge number
c ia        : mass number
c val       : real value
c ival      : integer value
c cval      : character value
c flagassign: flag to assign value or not
c
        if (key.eq.'elow') then
          read(value,*,end=1000,err=1000) eninclow
          goto 110
        endif
        if (key.eq.'massnucleus') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) massnucleus(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'massexcess') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) massexcess(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'a') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) alev(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'alimit') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) alimit(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'gammald') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gammald(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'pair') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) pair(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'pshift') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) pshift(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'deltaw') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) deltaW(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'exmatch') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Exmatch(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'t') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) T(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'e0') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) E0(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'nlow') then
          class=4
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Nlow(Zix,Nix,ibar)=ival
          goto 110
        endif
        if (key.eq.'ntop') then
          class=4
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Ntop(Zix,Nix,ibar)=ival
          goto 110
        endif
        if (key.eq.'s2adjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) s2adjust(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'krotconstant') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            Krotconstant(Zix,Nix,ibar)=val
            ldadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'ctable') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ctable(Zix,Nix,ibar)=val
            ldadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'ptable') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) ptable(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'g') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) g(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'gp') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gp(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'gn') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gp(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'egr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            egr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'ggr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ggr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'sgr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            sgr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'epr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            epr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'gpr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            gpr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'spr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            tpr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'egradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            egradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'ggradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ggradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'sgradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            sgradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'epradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            epradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'gpradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            gpradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'spradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            tpradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'fiso') then
          do 120 type=0,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 130
            endif
  120     continue
          goto 1000
  130     continue
          read(word(3),*,end=1000,err=1000) iz
          read(word(4),*,end=1000,err=1000) ia
          read(word(5),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            fiso(Zix,Nix,type2)=val
          endif
          goto 110
        endif
        if (key.eq.'gamgam') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gamgam(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'d0') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) D0(Zix,Nix)=val*1000.
          goto 110
        endif
        if (key.eq.'aadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) aadjust(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'tadjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Tadjust(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'e0adjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) E0adjust(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'exmatchadjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Exmatchadjust(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'gnadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gnadjust(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'gpadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gpadjust(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'gamgamadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gamgamadjust(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'etable') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            etable(Zix,Nix)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'ftable') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ftable(Zix,Nix)=val
            gamadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'fisbar') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fbarrier(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'fishw') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fwidth(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'fisbaradjust') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fbaradjust(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'fishwadjust') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fwidthadjust(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rtransmom') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Rtransmom(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'rclass2mom') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Rclass2mom(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'class2width') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) widthc2(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'levelfile') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) levelfile(Zix)=cval
          goto 110
        endif
        if (key.eq.'deformfile') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) deformfile(Zix)=cval
          goto 110
        endif
        if (key.eq.'hbtransfile') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) hbtransfile(Zix,Nix)=cval
          goto 110
        endif
        if (key.eq.'class2file') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) clas2file(Zix,Nix)=cval
          goto 110
        endif
        if (key.eq.'optmodfilen') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) optmodfilen(Zix)=cval
          goto 110
        endif
        if (key.eq.'optmodfilep') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) optmodfilep(Zix)=cval
          goto 110
        endif
        if (key.eq.'radialfile') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) radialfile(Zix)=cval
          goto 110
        endif
        if (key.eq.'ompenergyfile') then
          ompenergyfile=value
          goto 110
        endif
        if (key.eq.'beta2') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) beta2(Zix,Nix,ibar)=val
          goto 110
        endif
        if (key.eq.'axtype') then
          ibar=1
          class=4
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) axtype(Zix,Nix,ibar)=ival
          goto 110
        endif
        if (key.eq.'vfiscor') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) vfiscor(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'betafiscor') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) betafiscor(Zix,Nix)=val
          goto 110
        endif
        if (key.eq.'branch') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) ilev0
          read(word(5),*,end=1000,err=1000) nbr
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN)
     +      goto 1040
          if (ilev0.lt.0..or.ilev0.gt.numlev) goto 1060
          if (nbr.lt.0..or.nbr.gt.numlev) goto 1060
          iword=5
          sum=0.
          do 140 k=1,nbr
            iword=iword+1
            read(word(iword),*,end=1000,err=1000) ilev1
            if (ilev1.lt.0..or.ilev1.gt.numlev) goto 1060
            iword=iword+1
            read(word(iword),*,end=1000,err=1000) br
            if (br.lt.0.) goto 1070
            branchlevel(Zix,Nix,ilev0,k)=ilev1
            branchratio(Zix,Nix,ilev0,k)=br
            sum=sum+br
  140     continue
          if (sum.gt.0.) then
            do 150 k=1,nbr
              branchratio(Zix,Nix,ilev0,k)=branchratio(Zix,Nix,ilev0,k)
     +          /sum
  150       continue
          endif
          nbranch(Zix,Nix,ilev0)=nbr
          goto 110
        endif
        if (key.eq.'cstrip') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Cstrip(type)=val
          goto 110
        endif
        if (key.eq.'cknock') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Cknock(type)=val
          goto 110
        endif
        if (key.eq.'cbreak') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Cbreak(type)=val
          goto 110
        endif
        if (key.eq.'v1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'v2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'v3adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v3adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'v4adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v4adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rvadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rvadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'avadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            avadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'w1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'w2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'w3adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w3adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'w4adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w4adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rwadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rwadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'awadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            awadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rvdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rvdadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'avdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) avdadjust(type)=val
          ompadjustp(type)=.true.
          goto 110
        endif
        if (key.eq.'rwdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rwdadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'awdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            awdadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'d1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            d1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'d2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            d2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'d3adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            d3adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rvsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rvsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'avsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            avsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'vso1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            vso1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'vso2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            vso2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rwsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rwsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'awsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            awsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'wso1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            wso1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'wso2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            wso2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rcadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rcadjust(type)=val
            ompadjustp(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'ejoin') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Ejoin(type)=val
          goto 110
        endif
        if (key.eq.'vinfadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Vinfadjust(type)=val
          goto 110
        endif
        if (key.eq.'tljadjust') then
          class=8
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            tladjust(type,lval)=val
            tljadjust(type)=.true.
          endif
          goto 110
        endif
        if (key.eq.'rvadjustf'.or.key.eq.'avadjustf'.or.
     +    key.eq.'rwadjustf'.or.key.eq.'awadjustf'.or.
     +    key.eq.'rvdadjustf'.or.key.eq.'avdadjustf'.or.
     +    key.eq.'rwdadjustf'.or.key.eq.'awdadjustf'.or.
     +    key.eq.'rvsoadjustf'.or.key.eq.'avsoadjustf'.or.
     +    key.eq.'rwsoadjustf'.or.key.eq.'awsoadjustf') then
          if (key.eq.'rvadjustf') omptype=1
          if (key.eq.'avadjustf') omptype=2
          if (key.eq.'rwadjustf') omptype=3
          if (key.eq.'awadjustf') omptype=4
          if (key.eq.'rvdadjustf') omptype=5
          if (key.eq.'avdadjustf') omptype=6
          if (key.eq.'rwdadjustf') omptype=7
          if (key.eq.'awdadjustf') omptype=8
          if (key.eq.'rvsoadjustf') omptype=9
          if (key.eq.'avsoadjustf') omptype=10
          if (key.eq.'rwsoadjustf') omptype=11
          if (key.eq.'awsoadjustf') omptype=12
          if (key.eq.'rcadjustf') omptype=13
          do 160 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 170
            endif
  160     continue
          goto 1000
  170     ompadjustF(type2)=.true.
          ompadjustN(type2,omptype)=ompadjustN(type2,omptype)+1
          nr=ompadjustN(type2,omptype)
          read(word(3),*,end=1000,err=1000)
     +      ompadjustE1(type2,omptype,nr)
          read(word(4),*,end=1000,err=1000)
     +      ompadjustE2(type2,omptype,nr)
          read(word(5),*,end=1000,err=1000) ompadjustD(type2,omptype,nr)
          read(word(6),*,end=110,err=1000) ompadjusts(type2,omptype,nr)
          goto 110
        endif
        if (key.eq.'rescuefile') then
          read(word(2),*,end=1000,err=1000) mt
          if (mt.lt.1.or.mt.gt.nummt) goto 1050
          is=-1
          do 210 k=1,71
            if (word(3)(k:k+1).eq.'_g') then
              is=0
              goto 220
            endif
            if (word(3)(k:k+1).eq.'_m') then
              is=1
              goto 220
            endif
  210     continue
  220     rescuefile(mt,is)=word(3)
          val=1.
          read(word(4),*,end=230,err=1000) val
  230     grescue(mt,is)=val
          flagrescue=.true.
          goto 110
        endif
        if (key.eq.'jlmmode') then
          read(value,*,end=1000,err=1000) jlmmode
          goto 110
        endif
        if (key.eq.'lvadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lvadjust=val
            ompadjustp(1)=.true.
          endif
          goto 110
        endif
        if (key.eq.'lwadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lwadjust=val
            ompadjustp(1)=.true.
          endif
          goto 110
        endif
        if (key.eq.'lv1adjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lv1adjust=val
            ompadjustp(1)=.true.
          endif
          goto 110
        endif
        if (key.eq.'lw1adjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lw1adjust=val
            ompadjustp(1)=.true.
          endif
          goto 110
        endif
        if (key.eq.'lvsoadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lvsoadjust=val
            ompadjustp(1)=.true.
          endif
          goto 110
        endif
        if (key.eq.'lwsoadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lwsoadjust=val
            ompadjustp(1)=.true.
          endif
          goto 110
        endif
        if (key.eq.'aradialcor') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            aradialcor=val
            ompadjustp(6)=.true.
          endif
          goto 110
        endif
        if (key.eq.'adepthcor') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            adepthcor=val
            ompadjustp(6)=.true.
          endif
          goto 110
        endif
        if (key.eq.'gnorm') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gnorm=val
          goto 110
        endif
        if (key.eq.'spincutmodel') then
          read(value,*,end=1000,err=1000) spincutmodel
          goto 110
        endif
        if (key.eq.'shellmodel') then
          read(value,*,end=1000,err=1000) shellmodel
          goto 110
        endif
        if (key.eq.'kvibmodel') then
          read(value,*,end=1000,err=1000) kvibmodel
          goto 110
        endif
        if (key.eq.'radialmodel') then
          read(value,*,end=1000,err=1000) radialmodel
          goto 110
        endif
        if (key.eq.'rspincut') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Rspincut=val
          goto 110
        endif
        if (key.eq.'alphald') then
          read(value,*,end=1000,err=1000) alphaldall
          goto 110
        endif
        if (key.eq.'betald') then
          read(value,*,end=1000,err=1000) betaldall
          goto 110
        endif
        if (key.eq.'gammashell1') then
          read(value,*,end=1000,err=1000) gammashell1all
          goto 110
        endif
        if (key.eq.'gammashell2') then
          read(value,*,end=1000,err=1000) gammashell2
          goto 110
        endif
        if (key.eq.'pairconstant') then
          read(value,*,end=1000,err=1000) pairconstant
          goto 110
        endif
        if (key.eq.'pshiftconstant') then
          read(value,*,end=1000,err=1000) Pshiftconstantall
          goto 110
        endif
        if (key.eq.'ufermi') then
          read(value,*,end=1000,err=1000) Ufermi
          goto 110
        endif
        if (key.eq.'cfermi') then
          read(value,*,end=1000,err=1000) cfermi
          goto 110
        endif
        if (key.eq.'ufermibf') then
          read(value,*,end=1000,err=1000) Ufermibf
          goto 110
        endif
        if (key.eq.'cfermibf') then
          read(value,*,end=1000,err=1000) cfermibf
          goto 110
        endif
        if (key.eq.'kph') then
          read(value,*,end=1000,err=1000) Kph
          goto 110
        endif
        if (key.eq.'m2constant') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            M2constant=val
            preeqadjust=.true.
          endif
          goto 110
        endif
        if (key.eq.'m2limit') then
          read(value,*,end=1000,err=1000) M2limit
          goto 110
        endif
        if (key.eq.'m2shift') then
          read(value,*,end=1000,err=1000) M2shift
          goto 110
        endif
        if (key.eq.'rpipi') then
          read(value,*,end=1000,err=1000) Rpipi
          goto 110
        endif
        if (key.eq.'rnunu') then
          read(value,*,end=1000,err=1000) Rnunu
          goto 110
        endif
        if (key.eq.'rpinu') then
          read(value,*,end=1000,err=1000) Rpinu
          goto 110
        endif
        if (key.eq.'rnupi') then
          read(value,*,end=1000,err=1000) Rnupi
          goto 110
        endif
        if (key.eq.'esurf') then
          read(value,*,end=1000,err=1000) Esurf0
          goto 110
        endif
        if (key.eq.'rgamma') then
          read(value,*,end=1000,err=1000) Rgamma
          goto 110
        endif
        if (key.eq.'msdbins')  then
          read(value,*,end=1000,err=1000) msdbins
          goto 110
        endif
        if (key.eq.'emsdmin')  then
          read(value,*,end=1000,err=1000) Emsdmin
          goto 110
        endif
        if (key.eq.'elwidth')  then
          read(value,*,end=1000,err=1000) elwidth
          goto 110
        endif
        if (key.eq.'xscaptherm')  then
          read(value,*,end=1000,err=1000) xscaptherm
          goto 110
        endif
        if (key.eq.'xsptherm')  then
          read(value,*,end=1000,err=1000) xsptherm
          goto 110
        endif
        if (key.eq.'xsalphatherm')  then
          read(value,*,end=1000,err=1000) xsalphatherm
          goto 110
        endif
        if (key.eq.'anglesrec')  then
          read(value,*,end=1000,err=1000) nanglerec
          goto 110
        endif
        if (key.eq.'rprime')  then
          read(value,*,end=1000,err=1000) RprimeU
          goto 110
        endif
        if (key.eq.'astrot')  then
          read(value,*,end=1000,err=1000) astroT9
          nTmax=1
          flagastro=.true.
          goto 110
        endif
        if (key.eq.'astroe')  then
          read(value,*,end=1000,err=1000) astroE
          nTmax=1
          flagastro=.true.
          flagastrogs=.true.
          goto 110
        endif
        if (key.eq.'ebeam') then
          read(value,*,end=1000,err=1000) Ebeam
          goto 110
        endif
        if (key.eq.'eback') then
          read(value,*,end=1000,err=1000) Eback
          goto 110
        endif
        if (key.eq.'radiounit') then
          read(value,*,end=1000,err=1000) radiounit
          goto 110
        endif
        if (key.eq.'yieldunit') then
          read(value,*,end=1000,err=1000) yieldunit
          goto 110
        endif
        if (key.eq.'ibeam') then
          read(value,*,end=1000,err=1000) Ibeam
          goto 110
        endif
        if (key.eq.'area') then
          read(value,*,end=1000,err=1000) Area
          goto 110
        endif
        if (key.eq.'tirrad') then
          Tirrad(1)=0
          unitTirrad(1)=' '
          do 310 k=1,5
            read(word(2*k),'(i9)',end=110,err=1000) Tirrad(k)
            read(word(2*k+1),'(a1)',end=1000,err=1000) unitTirrad(k)
  310     continue
          goto 110
        endif
        if (key.eq.'rho') then
          read(value,*,end=1000,err=1000) rhotarget
          goto 110
        endif
        if (key.eq.'tcool') then
          Tcool(1)=0
          unitTcool(1)=' '
          do 320 k=1,5
            read(word(2*k),*,end=110,err=1000) Tcool(k)
            read(word(2*k+1),*,end=1000,err=1000) unitTcool(k)
  320     continue
          goto 110
        endif
        goto 110
 1040   write(*,'(" TALYS-warning: Z,N index out of range,",
     +    " keyword ignored: ",a80)') inline(i)
  110 continue
c
c Set level density parameters per nucleus
c
      if (alphaldall.ne.-99.) then
        do 910 Nix=0,numN
          do 910 Zix=0,numZ
            alphald(Zix,Nix)=alphaldall
  910   continue
      endif
      if (betaldall.ne.-99.) then
        do 920 Nix=0,numN
          do 920 Zix=0,numZ
            betald(Zix,Nix)=betaldall
  920   continue
      endif
      if (gammashell1all.ne.-99.) then
        do 930 Nix=0,numN
          do 930 Zix=0,numZ
            gammashell1(Zix,Nix)=gammashell1all
  930   continue
      endif
      if (Pshiftconstantall.ne.-99.) then
        do 940 Nix=0,numN
          do 940 Zix=0,numZ
            Pshiftconstant(Zix,Nix)=Pshiftconstantall
  940   continue
      endif
c
c Apply consistent OMP adjustment factors for Koning-Delaroche and other
c potentials.
c
      do 950 type=1,6
        if ((type.eq.3.and.deuteronomp.ge.2).or.(type.eq.6.and.
     +    alphaomp.ge.2)) then
          if (rwadjust(type).eq.-1.) rwadjust(type)=1.
          if (awadjust(type).eq.-1.) awadjust(type)=1.
          if (rvdadjust(type).eq.-1.) rvdadjust(type)=1.
          if (avdadjust(type).eq.-1.) avdadjust(type)=1.
          if (rwsoadjust(type).eq.-1.) rwsoadjust(type)=1.
          if (awsoadjust(type).eq.-1.) awsoadjust(type)=1.
        else
          if (rwadjust(type).eq.-1.) rwadjust(type)=rvadjust(type)
          if (awadjust(type).eq.-1.) awadjust(type)=avadjust(type)
          if (rvdadjust(type).eq.-1.) rvdadjust(type)=rwdadjust(type)
          if (avdadjust(type).eq.-1.) avdadjust(type)=awdadjust(type)
          if (rwsoadjust(type).eq.-1.) rwsoadjust(type)=rvsoadjust(type)
          if (awsoadjust(type).eq.-1.) awsoadjust(type)=avsoadjust(type)
        endif
  950 continue
      return

c Error and warning messages
c
 1000 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
 1050 write(*,'(" TALYS-error: 1 <= MT number <= ",i4,
     +  ", MT index out of range: ",a80)') nummt,inline(i)
      stop
 1060 write(*,'(" TALYS-error: 0 <= level number <= ",i4,
     +  ", ilev0, ilev1 or nbr index out of range: ",a80)')
     +  numlev,inline(i)
 1070 write(*,'(" TALYS-error: 0 <= branching ratio",
     +  ", br index out of range: ",a80)') inline(i)
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
