      subroutine input4
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 1, 2013
c | Task  : Read input for fourth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ch
      character*80 word(40),key,value
      integer      i,type
c
c ************** Defaults for fourth set of input variables ************
c
c flagmain    : flag for main output
c flagbasic   : flag for output of basic information and results
c flagpop     : flag for output of population
c flagcheck   : flag for output of numerical checks
c flagoutomp  : flag for output of optical model parameters
c flagdirect  : flag for output of direct reaction cross sections
c flaginverse : flag for output of transmission coefficients and
c               inverse reaction cross sections
c flaggamma   : flag for output of gamma-ray information
c flaglevels  : flag for output of discrete level information
c flagdensity : flag for output of level densities
c flagdisc    : flag for output of discrete state cross sections
c flagfisout  : flag for output of fission information
c flagfission : flag for fission
c flagtransen : flag for output of transmission coefficients per energy
c flagpeout   : flag for output of pre-equilibrium results
c flagang     : flag for output of angular distributions
c ddxmode     : mode for double-differential cross sections: 0: None,
c               1: Angular distributions, 2: Spectra per angle, 3: Both
c flaglegendre: flag for output of Legendre coefficients
c flagspec    : flag for output of spectra
c flagrecoil  : flag for calculation of recoils
c flagddx     : flag for output of double-differential cross sections
c flagoutdwba : flag for output of DWBA cross sections for MSD
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flageciscomp: flag for compound nucleus calculation by ECIS
c flagoutecis : flag for output of ECIS results
c flagompall  : flag for new optical model calculation for all
c               residual nuclei
c flagincadj  : flag for OMP adjustment on incident channel also
c flagecissave: flag for saving ECIS input and output files
c numinc      : number of incident energies
c flagexc     : flag for output of excitation functions
c flagnatural : flag for calculation of natural element
c eadd        : on-set incident energy for addition of discrete states
c               to spectra
c eaddel      : on-set incident energy for addition of elastic peak
c               to spectra
c Emaxtalys   : maximum acceptable energy for TALYS
c flagelectron: flag for application of electron conversion coefficient
c flagrot     : flag for use of rotational optical model per
c               outgoing particle, if available
c flagspher   : flag to force spherical optical model
c flagcoulomb : flag for Coulomb excitation calculation with ECIS
c flagcolldamp: flag for damping of collective effects in effective
c               level density (without explicit collective enhancement)
c               Only used for Bruyeres-le-Chatel (Pascal Romain) fission
c               model
c flagctmglob : flag for global CTM model (no discrete level info)
c cglobal     : global constant to adjust tabulated level densities
c pglobal     : global constant to adjust tabulated level densities
c alphaomp    : alpha optical model (1=normal, 2= McFadden-Satchler,
c               3-5= folding potential)
c deuteronomp : deuteron optical model (1=normal, 2=Daehnick,
c               3=Bojowald, 4=Han-Shi-Shen, 5=An-Cai)
c soswitch    : switch for deformed spin-orbit calculation and sequential
c               iterations in ECIS
c flagpartable: flag for output of model parameters on separate file
c maxchannel  : maximal number of outgoing particles in individual
c               channel description (e.g. this is 3 for (n,2np))
c Ztarget     : charge number of target nucleus
c massmodel   : model for theoretical nuclear mass
c pairmodel   : model for preequilibrium pairing energy
c flagmicro   : flag for completely microscopic Talys calculation
c fismodel    : fission model
c fismodelalt : alternative fission model for default barriers
c eurr        : off-set incident energy for URR calculation
c lurr        : maximal orbital angular momentum for URR
c flagurrnjoy : normalization of URR parameters with NJOY method
c flagendf    : flag for information for ENDF-6 file
c Atarget     : mass of target nucleus
c flagurr     : flag for output of unresolved resonance parameters
c k0          : index of incident particle
c flagcomp    : flag for compound nucleus calculation
c flagchannels: flag for exclusive channels calculation
c flagendfdet : flag for detailed ENDF-6 information per channel
c flagprod    : flag for isotope production
c flagoutfy   : flag for output detailed fission yield calculation
c gefran      : number of random events for GEF calculation
c
      flagmain=.true.
      flagpop=flagbasic
      flagcheck=flagbasic
      flagoutomp=flagbasic
      flagdirect=flagbasic
      flaginverse=flagbasic
      flaggamma=flagbasic
      flaglevels=flagbasic
      flagdensity=flagbasic
      flagdisc=flagbasic
      flagfisout=flagbasic
      if (.not.flagfission) flagfisout=.false.
      flagtransen=.true.
      flagpeout=.false.
      flagang=.false.
      ddxmode=0
      flaglegendre=.false.
      flagspec=.false.
      if (flagrecoil) flagspec=.true.
      flagddx=.false.
      flagoutdwba=.false.
      flaggamdis=.false.
      flagoutecis=flageciscomp
      if (flagompall) then
        flagecissave=.true.
      else
        flagecissave=.false.
      endif
      flagincadj=.true.
c
c By default, we assume that with more than one incident energy output
c of excitation functions (e.g. residual production cross sections as a
c function of incident energy) are wanted.
c
      if (numinc.eq.1) then
        flagexc=.false.
      else
        flagexc=.true.
      endif
      if (flagnatural) flagexc=.true.
      eadd=0.
      eaddel=0.
      flagelectron=.true.
      if (flagrot(k0)) then
        flagspher=.false.
      else
        flagspher=.true.
      endif
      flagcoulomb=.true.
      flagpartable=.false.
      maxchannel=4
      pairmodel=2
      if (flagmicro) then
        fismodel=5
      else
        fismodel=1
      endif
      fismodelalt=4
      flagcolldamp=.false.
      flagctmglob=.false.
      cglobal=1.e-20
      pglobal=1.e-20
      alphaomp=1
      deuteronomp=1
      soswitch=3.
      if (k0.ne.1.or..not.flagcomp) then
        eurr=0.
      else
        eurr=-1.
      endif
      flagurr=.false.
      lurr=2
      flagurrnjoy=.false.
c
c If the results of TALYS are used to create ENDF-6 data files,
c several output flags are automatically set.
c
      if (flagendf) then
        flagcheck=.true.
        flagdisc=.true.
        if (k0.eq.1) then
          if (Atarget.gt.20) flagurr=.true.
          eadd=20.
          eaddel=Emaxtalys
        endif
        if (flagfission) flagfisout=.true.
        if (k0.eq.3) then
          ddxmode=2
          flagddx=.true.
        endif
        flagang=.true.
        flaglegendre=.true.
        flagspec=.true.
        flagexc=.true.
        flagelectron=.true.
        if (flagendfdet) then
          flaggamdis=.true.
          flagchannels=.true.
        endif
      endif
      flagprod=.false.
      flagoutfy=.false.
      gefran=50000
c
c **************** Read fourth set of input variables ******************
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
      do 10 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
c
c Test for keywords
c
c
        if (key.eq.'outmain') then
          if (ch.eq.'n') flagmain=.false.
          if (ch.eq.'y') flagmain=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outpopulation') then
          if (ch.eq.'n') flagpop=.false.
          if (ch.eq.'y') flagpop=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outcheck') then
          if (ch.eq.'n') flagcheck=.false.
          if (ch.eq.'y') flagcheck=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outomp') then
          if (ch.eq.'n') flagoutomp=.false.
          if (ch.eq.'y') flagoutomp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outinverse') then
          if (ch.eq.'n') flaginverse=.false.
          if (ch.eq.'y') flaginverse=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outtransenergy') then
          if (ch.eq.'n') flagtransen=.false.
          if (ch.eq.'y') flagtransen=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outgamma') then
          if (ch.eq.'n') flaggamma=.false.
          if (ch.eq.'y') flaggamma=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outlevels') then
          if (ch.eq.'n') flaglevels=.false.
          if (ch.eq.'y') flaglevels=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outdensity') then
          if (ch.eq.'n') flagdensity=.false.
          if (ch.eq.'y') flagdensity=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outfission') then
          if (ch.eq.'n') flagfisout=.false.
          if (ch.eq.'y') flagfisout=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outpreequilibrium') then
          if (ch.eq.'n') flagpeout=.false.
          if (ch.eq.'y') flagpeout=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outdiscrete') then
          if (ch.eq.'n') flagdisc=.false.
          if (ch.eq.'y') flagdisc=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outspectra') then
          if (ch.eq.'n') flagspec=.false.
          if (ch.eq.'y') flagspec=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outangle') then
          if (ch.eq.'n') flagang=.false.
          if (ch.eq.'y') flagang=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outlegendre') then
          if (ch.eq.'n') flaglegendre=.false.
          if (ch.eq.'y') flaglegendre=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'ddxmode') then
          read(value,*,end=200,err=200) ddxmode
          if (ddxmode.eq.0) flagddx=.false.
          if (ddxmode.gt.0) then
            flagddx=.true.
            flagspec=.true.
          endif
          goto 10
        endif
        if (key.eq.'incadjust') then
          if (ch.eq.'n') flagincadj=.false.
          if (ch.eq.'y') flagincadj=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outdwba') then
          if (ch.eq.'n') flagoutdwba=.false.
          if (ch.eq.'y') flagoutdwba=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outgamdis') then
          if (ch.eq.'n') flaggamdis=.false.
          if (ch.eq.'y') flaggamdis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outexcitation') then
          if (ch.eq.'n') flagexc=.false.
          if (ch.eq.'y') flagexc=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outecis') then
          if (ch.eq.'n') flagoutecis=.false.
          if (ch.eq.'y') flagoutecis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'ecissave') then
          if (ch.eq.'n') flagecissave=.false.
          if (ch.eq.'y') flagecissave=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outdirect') then
          if (ch.eq.'n') flagdirect=.false.
          if (ch.eq.'y') flagdirect=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'adddiscrete') then
          if (ch.eq.'y') then
            eadd=0.
            goto 10
          endif
          if (ch.eq.'n') then
            eadd=Emaxtalys
            goto 10
          endif
          read(value,*,end=200,err=200) eadd
          goto 10
        endif
        if (key.eq.'addelastic') then
          if (ch.eq.'y') then
            eaddel=0.
            goto 10
          endif
          if (ch.eq.'n') then
            eaddel=Emaxtalys
            goto 10
          endif
          read(value,*,end=200,err=200) eaddel
          goto 10
        endif
        if (key.eq.'electronconv') then
          if (ch.eq.'n') flagelectron=.false.
          if (ch.eq.'y') flagelectron=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'spherical') then
          if (ch.eq.'n') flagspher=.false.
          if (ch.eq.'y') then
            flagspher=.true.
            do 110 type=1,6
              flagrot(type)=.false.
  110       continue
          endif
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'coulomb') then
          if (ch.eq.'n') flagcoulomb=.false.
          if (ch.eq.'y') flagcoulomb=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'colldamp') then
          if (ch.eq.'n') flagcolldamp=.false.
          if (ch.eq.'y') flagcolldamp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'ctmglobal') then
          if (ch.eq.'n') flagctmglob=.false.
          if (ch.eq.'y') flagctmglob=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'partable') then
          if (ch.eq.'n') flagpartable=.false.
          if (ch.eq.'y') flagpartable=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'urr') then
          flagurr=.true.
          if (ch.eq.'y') goto 10
          if (ch.eq.'n') then
            eurr=0.
            flagurr=.false.
            goto 10
          endif
          read(value,*,end=200,err=200) eurr
          goto 10
        endif
        if (key.eq.'lurr') then
          read(value,*,end=200,err=200) lurr
          goto 10
        endif
        if (key.eq.'urrnjoy') then
          if (ch.eq.'n') flagurrnjoy=.false.
          if (ch.eq.'y') flagurrnjoy=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'production') then
          if (ch.eq.'n') flagprod=.false.
          if (ch.eq.'y') flagprod=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outfy') then
          if (ch.eq.'n') flagoutfy=.false.
          if (ch.eq.'y') flagoutfy=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'maxchannel') then
          read(value,*,end=200,err=200) maxchannel
          goto 10
        endif
        if (key.eq.'pairmodel') then
          read(value,*,end=200,err=200) pairmodel
          goto 10
        endif
        if (key.eq.'fismodel') then
          read(value,*,end=200,err=200) fismodel
          goto 10
        endif
        if (key.eq.'fismodelalt') then
          read(value,*,end=200,err=200) fismodelalt
          goto 10
        endif
        if (key.eq.'cglobal') then
          read(value,*,end=200,err=200) cglobal
          goto 10
        endif
        if (key.eq.'pglobal') then
          read(value,*,end=200,err=200) pglobal
          goto 10
        endif
        if (key.eq.'alphaomp') then
          read(value,*,end=200,err=200) alphaomp
          goto 10
        endif
        if (key.eq.'deuteronomp') then
          read(value,*,end=200,err=200) deuteronomp
          goto 10
        endif
        if (key.eq.'soswitch') then
          read(value,*,end=200,err=200) soswitch
          goto 10
        endif
        if (key.eq.'gefran') then
          read(value,*,end=200,err=200) gefran
          goto 10
        endif
   10 continue
      return
  200 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
