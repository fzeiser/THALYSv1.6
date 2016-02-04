      subroutine input3
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Marieke Duijvestijn
c | Date  : December 2, 2013
c | Task  : Read input for third set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ch
      character*80 word(40),key,value
      integer      Zix,Nix,type,i,i2,iz,ia
c
c ************** Defaults for third set of input variables *************
c
c Some defaults are set on the basis of the input variables specified
c before. This can of course be overruled in the input file.
c
c flageciscalc: flag for new ECIS calculation for outgoing particles
c               and energy grid
c flaginccalc : flag for new ECIS calculation for incident channel
c flagendfecis: flag for new ECIS calculation for ENDF-6 files
c flagrel     : flag for relativistic kinematics
c flagcomp    : flag for compound nucleus calculation
c flagfullhf  : flag for full spin dependence of transmission
c               coefficients
c ewfc        : off-set incident energy for width fluctuation
c               calculation
c ptype0      : type of incident particle
c epreeq      : on-set incident energy for preequilibrium calculation
c Emaxtalys   : maximum acceptable energy for TALYS
c emulpre     : on-set incident energy for multiple preequilibrium
c               calculation
c pespinmodel : model for pre-equilibrium spin distribution or compound
c               spin distribution for pre-equilibrium cross section
c maxband     : highest vibrational level added to rotational model
c maxrot      : number of included excited rotational levels
c k0          : index for incident particle
c strength    : model for E1 gamma-ray strength function
c strengthM1  : model for M1 gamma-ray strength function
c flagpecomp  : flag for Kalbach complex particle emission model
c flagsurface : flag for surface effects in exciton model
c flaggiant0  : flag for collective contribution from giant resonances
c flag2comp   : flag for two-component pre-equilibrium model
c flagchannels: flag for exclusive channels calculation
c flagfission : flag for fission
c Atarget     : mass number of target nucleus
c fislim      : mass above which nuclide fissions
c ldmodelall  : level density model for all nuclides
c flagparity  : flag for non-equal parity distribution
c flaghbstate : flag for head band states in fission
c flagclass2  : flag for class2 states in fission
c flagbasic   : flag for output of basic information and results
c flageciscomp: flag for compound nucleus calculation by ECIS
c flagcpang   : flag for compound angular distribution calculation for
c               incident charged particles
c flagecisdwba: flag for new ECIS calculation for DWBA for MSD
c flagonestep : flag for continuum one-step direct only
c flaglocalomp: flag for local (y) or global (n) optical model
c flagdisp    : flag for dispersive optical model
c flagompall  : flag for new optical model calculation for all
c               residual nuclei
c flagautorot : flag for automatic rotational coupled channels
c               calculations for A > 150
c flagstate   : flag for optical model potential for each excited state
c Zix         : charge number index for residual nucleus
c numZph      : maximal number of protons away from the initial
c               compound nucleus for multiple pre-equilibrium emission
c Nix         : neutron number index for residual nucleus
c numNph      : maximal number of neutrons away from the initial
c               compound nucleus for multiple pre-equilibrium emission
c optmod      : file with optical model parameters
c flagsys     : flag for reaction cross section from systematics
c flagrot     : flag for use of rotational optical model per
c               outgoing particle, if available
c flagasys    : flag for all level density parameters a from systematics
c flaggshell  : flag for energy dependence of single particle level
c               density parameter g
c flagmassdis : flag for calculation of fission fragment mass yields
c flagffevap  : flag for calculation of particle evaporation from
c               fission fragment mass yields
c fymodel     : fission yield model, 1: Brosa 2: GEF
c flagfisfeed : flag for output of fission per excitation bin
c flagendf    : flag for information for ENDF-6 file
c flagendfdet : flag for detailed ENDF-6 information per channel
c flagrecoil  : flag for calculation of recoils
c flaglabddx  : flag for calculation of DDX in LAB system
c flagrecoilav: flag for average velocity in recoil calculation
c flagEchannel: flag for channel energy for emission spectrum
c flagreaction: flag for calculation of nuclear reactions
c flagastrogs : flag for calculation of astrophysics reaction rate
c               with target in ground state only
c flagexpmass : flag for using experimental nuclear mass if available
c flagjlm     : flag for using semi-microscopic JLM OMP
c flagomponly : flag to execute ONLY an optical model calculation
c flagmicro   : flag for completely microscopic Talys calculation
c
      flageciscalc=.true.
      flaginccalc=.true.
      flagendfecis=.true.
      flagrel=.true.
      flagcomp=.true.
      flagfullhf=.false.
      ewfc=-1.
      if (ptype0.eq.'0') then
        epreeq=Emaxtalys
      else
        epreeq=-1.
      endif
      emulpre=20.
      if (k0.le.1) then
        pespinmodel=1
      else
        pespinmodel=2
      endif
      maxband=0
      maxrot=2
      if (k0.ge.1) then
        strength=1
        flagpecomp=.true.
        flagsurface=.true.
      else
        strength=2
        flagpecomp=.false.
        flagsurface=.false.
      endif
      strengthM1=2
      if (k0.eq.1.or.k0.eq.2) then
        flaggiant0=.true.
      else
        flaggiant0=.false.
      endif
      flag2comp=.true.
      flagchannels=.false.
      flagfission=.false.
      if (Atarget.gt.fislim) flagfission=.true.
c
c ldmodel=5 (Combinatorial HFB model) has parity-dependent
c level densities.
c
      if (ldmodelall.ge.5) then
        flagparity=.true.
      else
        flagparity=.false.
      endif
      flaghbstate=.true.
      flagclass2=.true.
      flagbasic=.false.
      flageciscomp=.false.
      flagcpang=.false.
      flagecisdwba=.true.
      flagonestep=.false.
      flaglocalomp=.true.
      flagdisp=.false.
      flagompall=.false.
      flagautorot=.false.
      flagstate=.false.
      do 10 Zix=0,numZph
        do 10 Nix=0,numNph
          do 10 type=1,6
            optmod(Zix,Nix,type)=
     +        '                                                        '
   10 continue
      do 20 type=0,6
        flagsys(type)=.false.
   20 continue
      do 30 type=0,6
        flagrot(type)=.false.
   30 continue
      if (k0.eq.1.or.k0.eq.2) flagrot(k0)=.true.
      flagasys=.false.
      flaggshell=.false.
      flagmassdis=.false.
      flagffevap=.true.
      fymodel=2
      flagfisfeed=.false.
      flagendf=.false.
      if (k0.le.1) then
        flagendfdet=.true.
      else
        flagendfdet=.false.
      endif
      flagrecoil=.false.
      flaglabddx=.false.
      flagrecoilav=.false.
      flagEchannel=.false.
      flagreaction=.true.
      flagastrogs=.false.
      flagexpmass=.true.
      flagjlm=.false.
      if (flagomponly) then
        flagcomp=.false.
        epreeq=Emaxtalys
        emulpre=Emaxtalys
        flaggiant0=.false.
      endif
      if (flagmicro) then
        strength=4
        flagautorot=.true.
        flagjlm=.true.
      endif
c
c **************** Read third set of input variables *******************
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
c
c Test for keywords
c
c Zinit : charge number of initial compound nucleus
c Ninit : neutron number of initial compound nucleus
c parsym: symbol of particle
c
        if (key.eq.'maxband') then
          read(value,*,end=300,err=300) maxband
          goto 110
        endif
        if (key.eq.'maxrot')  then
          read(value,*,end=300,err=300) maxrot
          goto 110
        endif
        if (key.eq.'strength') then
          read(value,*,end=300,err=300) strength
          goto 110
        endif
        if (key.eq.'strengthM1') then
          read(value,*,end=300,err=300) strengthM1
          goto 110
        endif
        if (key.eq.'eciscalc') then
          if (ch.eq.'n') flageciscalc=.false.
          if (ch.eq.'y') flageciscalc=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'inccalc') then
          if (ch.eq.'n') flaginccalc=.false.
          if (ch.eq.'y') flaginccalc=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'endfecis') then
          if (ch.eq.'n') flagendfecis=.false.
          if (ch.eq.'y') flagendfecis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'relativistic') then
          if (ch.eq.'n') flagrel=.false.
          if (ch.eq.'y') flagrel=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'compound') then
          if (ch.eq.'n') flagcomp=.false.
          if (ch.eq.'y') flagcomp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'widthfluc') then
          if (ch.eq.'y') then
            if (k0.gt.1) ewfc=10.
            goto 110
          endif
          if (ch.eq.'n') then
            ewfc=0.
            goto 110
          endif
          read(value,*,end=300,err=300) ewfc
          goto 110
        endif
        if (key.eq.'preequilibrium') then
          if (ch.eq.'y') then
            epreeq=0.
            goto 110
          endif
          if (ch.eq.'n') then
            epreeq=Emaxtalys
            goto 110
          endif
          read(value,*,end=300,err=300) epreeq
          goto 110
        endif
        if (key.eq.'multipreeq') then
          if (ch.eq.'y') then
            emulpre=0.
            goto 110
          endif
          if (ch.eq.'n') then
            emulpre=Emaxtalys
            goto 110
          endif
          read(value,*,end=300,err=300) emulpre
          goto 110
        endif
        if (key.eq.'preeqspin') then
          if (ch.eq.'n') then
            pespinmodel=1
            goto 110
          endif
          if (ch.eq.'y') then
            pespinmodel=3
            goto 110
          endif
          read(value,*,end=300,err=300) pespinmodel
          goto 110
        endif
        if (key.eq.'giantresonance') then
          if (ch.eq.'n') flaggiant0=.false.
          if (ch.eq.'y') flaggiant0=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'preeqsurface') then
          if (ch.eq.'n') flagsurface=.false.
          if (ch.eq.'y') flagsurface=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'preeqcomplex') then
          if (ch.eq.'n') flagpecomp=.false.
          if (ch.eq.'y') flagpecomp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'twocomponent') then
          if (ch.eq.'n') flag2comp=.false.
          if (ch.eq.'y') flag2comp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'channels') then
          if (ch.eq.'n') flagchannels=.false.
          if (ch.eq.'y') flagchannels=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'fission') then
          if (ch.eq.'n') flagfission=.false.
          if (ch.eq.'y') flagfission=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'hbstate') then
          if (ch.eq.'n') flaghbstate=.false.
          if (ch.eq.'y') flaghbstate=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'class2') then
          if (ch.eq.'n') flagclass2=.false.
          if (ch.eq.'y') flagclass2=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'outbasic') then
          if (ch.eq.'n') flagbasic=.false.
          if (ch.eq.'y') flagbasic=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'eciscompound') then
          if (ch.eq.'n') flageciscomp=.false.
          if (ch.eq.'y') flageciscomp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'cpang') then
          if (ch.eq.'n') flagcpang=.false.
          if (ch.eq.'y') flagcpang=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'ecisdwba') then
          if (ch.eq.'n') flagecisdwba=.false.
          if (ch.eq.'y') flagecisdwba=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'onestep') then
          if (ch.eq.'n') flagonestep=.false.
          if (ch.eq.'y') flagonestep=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'localomp') then
          if (ch.eq.'n') flaglocalomp=.false.
          if (ch.eq.'y') flaglocalomp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'dispersion') then
          if (ch.eq.'n') flagdisp=.false.
          if (ch.eq.'y') flagdisp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'optmodall') then
          if (ch.eq.'n') flagompall=.false.
          if (ch.eq.'y') flagompall=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'statepot') then
          if (ch.eq.'n') flagstate=.false.
          if (ch.eq.'y') flagstate=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'optmod') then
          read(word(2),*,end=300,err=300) iz
          read(word(3),*,end=300,err=300) ia
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZph.or.
     +      Nix.lt.0.or.Nix.gt.numNph) then
            write(*,'(" TALYS-warning: Z,N index out of range,",
     +        " keyword ignored: ",a80)') inline(i)
            goto 110
          else
            ch=word(5)(1:1)
            if (ch.eq.' ') ch='n'
            do 210 type=1,6
              if (ch.eq.parsym(type)) then
                optmod(Zix,Nix,type)=word(4)
                goto 110
              endif
  210       continue
          endif
        endif
        if (key.eq.'sysreaction') then
          do 220 type=0,6
            flagsys(type)=.false.
  220     continue
          do 230 i2=2,40
            ch=word(i2)(1:1)
            do 240 type=0,6
              if (ch.eq.parsym(type)) then
                flagsys(type)=.true.
                goto 230
              endif
  240       continue
  230     continue
          goto 110
        endif
        if (key.eq.'rotational') then
          do 250 type=1,6
            flagrot(type)=.false.
  250     continue
          do 260 i2=2,40
            ch=word(i2)(1:1)
            do 270 type=1,6
              if (ch.eq.parsym(type)) then
                flagrot(type)=.true.
                goto 260
              endif
  270       continue
  260     continue
          goto 110
        endif
        if (key.eq.'asys') then
          if (ch.eq.'n') flagasys=.false.
          if (ch.eq.'y') flagasys=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'gshell') then
          if (ch.eq.'n') flaggshell=.false.
          if (ch.eq.'y') flaggshell=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'massdis') then
          if (ch.eq.'n') flagmassdis=.false.
          if (ch.eq.'y') flagmassdis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'ffevaporation') then
          if (ch.eq.'n') flagffevap=.false.
          if (ch.eq.'y') flagffevap=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'fisfeed') then
          if (ch.eq.'n') flagfisfeed=.false.
          if (ch.eq.'y') flagfisfeed=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'fymodel') then
          read(value,*,end=300,err=300) fymodel
          goto 110
        endif
        if (key.eq.'endf') then
          if (ch.eq.'n') flagendf=.false.
          if (ch.eq.'y') flagendf=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'endfdetail') then
          if (ch.eq.'n') flagendfdet=.false.
          if (ch.eq.'y') flagendfdet=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'recoil') then
          if (ch.eq.'n') flagrecoil=.false.
          if (ch.eq.'y') flagrecoil=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'labddx') then
          if (ch.eq.'n') flaglabddx=.false.
          if (ch.eq.'y') flaglabddx=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'recoilaverage') then
          if (ch.eq.'n') flagrecoilav=.false.
          if (ch.eq.'y') flagrecoilav=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'channelenergy') then
          if (ch.eq.'n') flagEchannel=.false.
          if (ch.eq.'y') flagEchannel=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'fullhf') then
          if (ch.eq.'n') flagfullhf=.false.
          if (ch.eq.'y') flagfullhf=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'autorot') then
          if (ch.eq.'n') flagautorot=.false.
          if (ch.eq.'y') flagautorot=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'reaction') then
          if (ch.eq.'n') flagreaction=.false.
          if (ch.eq.'y') flagreaction=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'parity') then
          if (ch.eq.'n') flagparity=.false.
          if (ch.eq.'y') flagparity=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'astrogs') then
          if (ch.eq.'n') flagastrogs=.false.
          if (ch.eq.'y') flagastrogs=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'expmass') then
          if (ch.eq.'n') flagexpmass=.false.
          if (ch.eq.'y') flagexpmass=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'jlmomp') then
          if (ch.eq.'n') flagjlm=.false.
          if (ch.eq.'y') flagjlm=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
  110 continue
      return
  300 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
