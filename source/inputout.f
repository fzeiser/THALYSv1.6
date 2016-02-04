      subroutine inputout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 1, 2013
c | Task  : Write input parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  yesno
      character*12 sysstring,rotstring
      integer      i,type
c
c ************************** User input file ***************************
c
c nlines: number of input lines
c inline: input line
c
      write(*,'(/" ########## USER INPUT ##########")')
      write(*,'(/" USER INPUT FILE"/)')
      do 10 i=1,nlines
        write(*,'(1x,a80)') inline(i)
   10 continue
c
c ********* All possible input parameters including defaults ***********
c
      write(*,'(/" USER INPUT FILE + DEFAULTS"/)')
      write(*,'(" Keyword           Value   Variable",
     +  "     Explanation"/)')
c
c 1. Four main keywords
c
c ptype0     : type of incident particle
c Starget    : symbol of target nucleus
c Atarget    : mass number of target nucleus
c numinc     : number of incident energies
c flaginitpop: flag for initial population distribution
c eninc      : incident energy in MeV
c energyfile : file with incident energies
c
      write(*,'(" #"/" # Four main keywords"/" #")')
      write(*,'(" projectile          ",a1,"     ptype0",
     +  "       type of incident particle")') ptype0
      write(*,'(" element            ",a2,"     Starget",
     +  "      symbol of target nucleus")') Starget
      write(*,'(" mass              ",i3,"     mass     ",
     +  "    mass number of target nucleus")') Atarget
      if (numinc.eq.1.and..not.flaginitpop) then
        write(*,'(" energy            ",f7.3," eninc",
     +    "        incident energy in MeV")') eninc(1)
      else
        write(*,'(" energy ",a14,"     energyfile",
     +    "   file with incident energies")') energyfile
      endif
c
c 2. Basic physical and numerical parameters
c
c outtype     : type of outgoing particles
c maxZ        : maximal number of protons away from the initial compound
c               nucleus
c maxN        : maximal number of neutrons away from the initial
c               compound nucleus
c nbins0      : number of continuum excitation energy bins
c flagequi    : flag to use equidistant excitation bins instead of
c               logarithmic bins
c segment     : number of segments to divide emission energy grid
c nlevmax     : maximum number of included discrete levels for target
c nlevmaxres  : maximum number of included discrete levels for residual
c               nucleus
c parsym      : symbol of particle
c nlevbin     : number of excited levels for binary nucleus
c parname     : name of particle
c Ltarget     : excited level of target
c isomer      : definition of isomer in seconds
c transpower  : power for transmission coefficient limit
c transeps    : absolute limit for transmission coefficient
c xseps       : limit for cross sections
c popeps      : limit for population cross section per nucleus
c Rfiseps     : ratio for limit for fission cross section per nucleus
c eninclow    : minimal incident energy for nuclear model calculations
c nangle      : number of angles
c nanglecont  : number of angles for continuum
c nanglerec   : number of recoil angles
c maxenrec    : number of recoil energies
c yesno       : function to assign y or n to logical value
c flagchannels: flag for exclusive channels calculation
c maxchannel  : maximal number of outgoing particles in individual
c               channel description (e.g. this is 3 for (n,2np))
c flagmicro   : flag for completely microscopic Talys calculation
c flagbest    : flag to use best set of adjusted parameters
c flagrel     : flag for relativistic kinematics
c flagrecoil  : flag for calculation of recoils
c flaglabddx  : flag for calculation of DDX in LAB system
c flagrecoilav: flag for average velocity in recoil calculation
c flagEchannel: flag for channel energy for emission spectrum
c flagreaction: flag for calculation of nuclear reactions
c flagastro   : flag for calculation of astrophysics reaction rate
c flagastrogs : flag for calculation of astrophysics reaction rate with
c               target in ground state only
c massmodel   : model for theoretical nuclear mass
c flagexpmass : flag for using experimental nuclear mass if available
c disctable   : table with discrete levels
c flagprod    : flag for isotope production
c flagoutfy   : flag for output detailed fission yield calculation
c gefran      : number of random events for GEF calculation
c
      write(*,'(" #"/" # Basic physical and numerical parameters")')
      write(*,'(" #")')
      write(*,'(" ejectiles ",7(1x,a1),"  outtype      ",
     +  "outgoing particles")') (outtype(type),type=0,6)
      write(*,'(" maxz              ",i3,"     maxZ         ",
     +  "maximal number of protons away from the initial",
     +  " compound nucleus")') maxZ
      write(*,'(" maxn              ",i3,"     maxN         ",
     +  "maximal number of neutrons away from the initial",
     +  " compound nucleus")') maxN
      write(*,'(" bins              ",i3,"     nbins        ",
     +  "number of continuum excitation energy bins")') nbins0
      write(*,'(" equidistant         ",a1,"     flagequi    ",
     +  " flag to use equidistant excitation bins instead of ",
     +  "logarithmic bins")') yesno(flagequi)
      write(*,'(" segment           ",i3,"     segment    ",
     +  "  number of segments to divide emission energy grid")') segment
      write(*,'(" maxlevelstar      ",i3,"     nlevmax",
     +  "      maximum number of included discrete levels for target")')
     +  nlevmax
      write(*,'(" maxlevelsres      ",i3,"     nlevmaxres",
     +  "   maximum number of included discrete levels",
     +  " for residual nucleus")') nlevmaxres
      do 20 type=0,6
        write(*,'(" maxlevelsbin ",a1,"    ",i3,"     nlevbin   ",
     +    "   maximum number of included discrete levels for ",
     +    a8," channel")') parsym(type),nlevbin(type),parname(type)
   20 continue
      write(*,'(" ltarget           ",i3,"     ltarget",
     +  "      excited level of target")') Ltarget
      write(*,'(" isomer          ",1p,e9.2," isomer ",
     +  "      definition of isomer in seconds")') isomer
      write(*,'(" transpower        ",i3,"     transpower",
     +  "   power for transmission coefficient limit")') transpower
      write(*,'(" transeps        ",1p,e9.2," transeps",
     +  "     limit for transmission coefficient")') transeps
      write(*,'(" xseps           ",1p,e9.2," xseps   ",
     +  "     limit for cross sections")') xseps
      write(*,'(" popeps          ",1p,e9.2," popeps  ",
     +  "     limit for population cross section per nucleus")') popeps
      write(*,'(" Rfiseps         ",1p,e9.2," Rfiseps      ratio",
     +  " for limit for fission cross section per nucleus")') Rfiseps
      write(*,'(" elow            ",1p,e9.2," elow         minimal",
     +  " incident energy for nuclear model calculations")') eninclow
      write(*,'(" angles            ",i3,"     nangle",
     +  "       number of angles")') nangle
      write(*,'(" anglescont        ",i3,"     nanglecont",
     +  "   number of angles for continuum")') nanglecont
      write(*,'(" anglesrec         ",i3,"     nanglerec ",
     +  "   number of recoil angles")') nanglerec
      write(*,'(" maxenrec          ",i3,"     maxenrec  ",
     +  "   number of recoil energies")') maxenrec
      write(*,'(" channels            ",a1,"     flagchannels flag",
     +  " for exclusive channels calculation")') yesno(flagchannels)
      write(*,'(" maxchannel         ",i2,"     maxchannel",
     +  "   maximal number of outgoing particles in",
     +  " individual channel description")') maxchannel
      write(*,'(" micro               ",a1,"     flagmicro",
     +  "    flag for completely microscopic Talys calculation")')
     +  yesno(flagmicro)
      write(*,'(" best                ",a1,"     flagbest ",
     +  "    flag to use best set of adjusted parameters")')
     +  yesno(flagbest)
      write(*,'(" relativistic        ",a1,"     flagrel",
     +  "      flag for relativistic kinematics")') yesno(flagrel)
      write(*,'(" recoil              ",a1,"     flagrecoil",
     +  "   flag for calculation of recoils")') yesno(flagrecoil)
      write(*,'(" labddx              ",a1,"     flaglabddx   flag",
     +  " for calculation of DDX in LAB system")') yesno(flaglabddx)
      write(*,'(" recoilaverage       ",a1,"     flagrecoilav flag for",
     +  " average velocity in recoil calculation")') yesno(flagrecoilav)
      write(*,'(" channelenergy       ",a1,"     flagEchannel flag for",
     +  " channel energy for emission spectrum")') yesno(flagEchannel)
      write(*,'(" reaction            ",a1,"     flagreaction flag",
     +  " for calculation of nuclear reactions")') yesno(flagreaction)
      write(*,'(" astro               ",a1,"     flagastro    flag for",
     +  " calculation of astrophysics reaction rate")') yesno(flagastro)
      write(*,'(" astrogs             ",a1,"     flagastrogs  flag for",
     +  " calculation of astrophysics reaction rate with target in",
     +  " ground state only")') yesno(flagastrogs)
      write(*,'(" massmodel          ",i2,"     massmodel",
     +  "    model for theoretical nuclear mass")') massmodel
      write(*,'(" expmass             ",a1,"     flagexpmass ",
     +  " flag for using experimental nuclear mass if available")')
     +  yesno(flagexpmass)
      write(*,'(" disctable          ",i2,"     disctable",
     +  "    table with discrete levels")') disctable
      write(*,'(" production          ",a1,"     flagprod    ",
     +  " flag for isotope production")') yesno(flagprod)
      write(*,'(" outfy               ",a1,"     flagoutfy   ",
     +  " flag for output detailed fission yield calculation")')
     +  yesno(flagoutfy)
      write(*,'(" gefran        ",i7,"     gefran       number ",
     +  "of random events for GEF calculation")') gefran
c
c Isotope production
c
c Ebeam     : incident energy in MeV for isotope production
c Eback     : lower end of energy range in MeV for isotope production
c radiounit : unit for radioactivity: Bq, kBq, MBq, Gbq,
c             mCi, Ci or kCi
c yieldunit : unit for isotope yield: num (number),
c             mug (micro-gram), mg, g, or kg
c Ibeam     : beam current in mA
c Tirrad    : irradiation time per unit
c unitTirrad: irradiation time unit (y,d,h,m,s)
c Area      : target area in cm^2
c Tcool     : cooling time per unit
c unitTcool : cooling time unit (y,d,h,m,s)
c rhotarget : target density
c
      if (flagprod) then
        write(*,'(" #"/" # Isotope production "/" #")')
        write(*,'(" Ebeam             ",f7.3," beam    ",
     +    "     incident energy in MeV for isotope production")') Ebeam
        write(*,'(" Eback             ",f7.3," Eback  ",
     +    "      lower end of energy range in MeV for isotope",
     +    "  production")') Eback
        write(*,'(" radiounit             ",a3," radiounit ",
     +    "   unit for radioactivity")') radiounit
        write(*,'(" yieldunit             ",a3," yieldunit ",
     +    "   unit for isotope yield")') yieldunit
        write(*,'(" Ibeam             ",f7.3," Ibeam ",
     +    "       beam current in mA")') Ibeam
        do 30 i=1,5
          if (Tirrad(i).gt.0) write(*,'(" Tirrad      ",i9,
     +      "     Tirrad       ",a1," of irradiation time")')
     +      Tirrad(i),unitTirrad(i)
   30   continue
        write(*,'(" Area              ",f7.3," Area  ",
     +    "       target area in cm^2")') Area
        do 40 i=1,5
          if (Tcool(i).gt.0) write(*,'(" Tcool       ",i9,"     Tcool ",
     +      "       ",a1," of cooling time")') Tcool(i),unitTcool(i)
   40   continue
        write(*,'(" rho               ",f7.3," rhotarget",
     +    "    target density [g/cm^3] ")') rhotarget
      endif
c
c 3. Optical model
c
c flaglocalomp: flag for local (y) or global (n) optical model
c flagdisp    : flag for dispersive optical model
c flagjlm     : flag for using semi-microscopic JLM OMP
c flagompall  : flag for new optical model calculation for all residual
c               nuclei
c flagincadj  : flag for OMP adjustment on incident channel also
c flagomponly : flag to execute ONLY an optical model calculation
c flagautorot : flag for automatic rotational coupled channels
c               calculations for A > 150
c flagspher   : flag to force spherical optical model
c flagcoulomb : flag for Coulomb excitation calculation with ECIS
c flagstate   : flag for optical model potential for each excited state
c maxband     : highest vibrational level added to rotational model
c maxrot      : number of included excited rotational levels
c sysstring   : help variable
c flagsys     : flag for reaction cross section from systematics
c rotstring   : help variable
c flagrot     : flag for use of rotational optical model per
c               outgoing particle, if available
c core        : even-even core for weakcoupling (-1 or 1)
c flagecissave: flag for saving ECIS input and output files
c flageciscalc: flag for new ECIS calculation for outgoing particles
c               and energy grid
c flaginccalc : flag for new ECIS calculation for incident channel
c flagendfecis: flag for new ECIS calculation for ENDF-6 files
c radialmodel : model for radial matter densities (JLM OMP only)
c jlmmode     : option for JLM imaginary potential normalization
c alphaomp    : alpha optical model (1=normal, 2= McFadden-Satchler,
c               3-5= folding potential)
c deuteronomp : deuteron optical model (1=normal, 2=Daehnick,
c               3=Bojowald, 4=Han-Shi-Shen, 5=An-Cai)
c
      write(*,'(" #"/" # Optical model"/" #")')
      write(*,'(" localomp            ",a1,"     flaglocalomp flag for",
     +  " local (y) or global (n) optical model")') yesno(flaglocalomp)
      write(*,'(" dispersion          ",a1,"     flagdisp    ",
     +  " flag for dispersive optical model")') yesno(flagdisp)
      write(*,'(" jlmomp              ",a1,"     flagjlm      flag for",
     +  " using semi-microscopic JLM OMP")') yesno(flagjlm)
      write(*,'(" optmodall           ",a1,"     flagompall   flag for",
     +  " new optical model calculation for all residual nuclei")')
     +  yesno(flagompall)
      write(*,'(" incadjust           ",a1,"     flagincadj   flag for",
     +  " OMP adjustment on incident channel also")') yesno(flagincadj)
      write(*,'(" omponly             ",a1,"     flagomponly ",
     +  " flag to execute ONLY an optical model calculation")')
     +  yesno(flagomponly)
      write(*,'(" autorot             ",a1,"     flagautorot ",
     +  " flag for automatic rotational coupled channels ",
     +  "calculations for A > 150")') yesno(flagautorot)
      write(*,'(" spherical           ",a1,"     flagspher   ",
     +  " flag to force spherical optical model")') yesno(flagspher)
      write(*,'(" coulomb             ",a1,"     flagcoulomb ",
     +  " flag for Coulomb excitation calculation with ECIS")')
     +  yesno(flagcoulomb)
      write(*,'(" statepot            ",a1,"     flagstate   ",
     +  " flag for optical model potential for each excited state")')
     +  yesno(flagstate)
      write(*,'(" maxband           ",i3,"     maxband      highest",
     +  " vibrational band added to rotational model")') maxband
      write(*,'(" maxrot            ",i3,"     maxrot      ",
     +  " number of included excited rotational levels")') maxrot
      sysstring='            '
      i=-1
      do 110 type=1,6
        if (flagsys(type)) then
          i=i+2
          write(sysstring(i:i),'(a1)') parsym(type)
        endif
  110 continue
      write(*,'(" sysreaction  ",a12," sysreaction  particles",
     +  " with reaction cross section from systematics")') sysstring
      rotstring='            '
      i=-1
      do 120 type=1,6
        if (flagrot(type)) then
          i=i+2
          write(rotstring(i:i),'(a1)') parsym(type)
        endif
  120 continue
      write(*,'(" rotational   ",a12," rotational   ",
     +  "particles with possible rotational optical model")') rotstring
      write(*,'(" core              ",i3,"     core   ",
     +  "      even-even core for weakcoupling (-1 or 1)")') core
      write(*,'(" ecissave            ",a1,"     flagecissave flag",
     +  " for saving ECIS input and output files")') yesno(flagecissave)
      write(*,'(" eciscalc            ",a1,"     flageciscalc",
     +  " flag for new ECIS calculation for outgoing particles and",
     +  " energy grid")') yesno(flageciscalc)
      write(*,'(" inccalc             ",a1,"     flaginccalc ",
     +  " flag for new ECIS calculation for incident channel")')
     +  yesno(flaginccalc)
      write(*,'(" endfecis            ",a1,"     flagendfecis",
     +  " flag for new ECIS calculation for ENDF-6 files")')
     +  yesno(flagendfecis)
      write(*,'(" radialmodel        ",i2,"     radialmodel ",
     +  " model for radial matter densities (JLM OMP only)")')
     +  radialmodel
      write(*,'(" jlmmode            ",i2,"     jlmmode     ",
     +  " option for JLM imaginary potential normalization")')
     +  jlmmode
      write(*,'(" alphaomp           ",i2,"     alphaomp    ",
     +  " alpha optical model (1=normal, 2= McFadden-Satchler,",
     +  " 3-5= folding potential)")') alphaomp
      write(*,'(" deuteronomp        ",i2,"     deuteronomp ",
     +  " deuteron optical model (1=normal, 2=Daehnick,",
     +  " 3=Bojowald, 4=Han-Shi-Shen, 5=An-Cai)")') deuteronomp
c
c 4. Compound nucleus
c
c enincmin    : minimum incident energy
c ewfc        : off-set incident energy for width fluctuation
c               calculation
c enincmax    : maximum incident energy
c flagwidth   : flag for width fluctuation calculation
c wmode       : designator for width fluctuation model
c flagcomp    : flag for compound nucleus calculation
c flagfullhf  : flag for full spin dependence of transmission
c               coefficients
c flageciscomp: flag for compound nucleus calculation by ECIS
c flagcpang   : flag for compound angular distribution calculation for
c               incident charged particles
c eurr        : off-set incident energy for URR calculation
c flagurr     : flag for output of unresolved resonance parameters
c lurr        : maximal orbital angular momentum for URR
c flagurrnjoy : normalization of URR parameters with NJOY method
c
      write(*,'(" #"/" # Compound nucleus"/" #")')
      if (numinc.gt.1.and.enincmin.lt.ewfc.and.enincmax.ge.ewfc) then
        write(*,'(" widthfluc         ",f7.3," ewfc         off-set",
     +    " incident energy for width fluctuation calculation")') ewfc
      else
        write(*,'(" widthfluc           ",a1,"     flagwidth  ",
     +    "  flag for width fluctuation calculation")') yesno(flagwidth)
      endif
      write(*,'(" widthmode          ",i2,"     wmode      ",
     +  "  designator for width fluctuation model")') wmode
      write(*,'(" compound            ",a1,"     flagcomp     ",
     +  "flag for compound nucleus model")') yesno(flagcomp)
      write(*,'(" fullhf              ",a1,"     flagfullhf   ",
     +  "flag for full spin dependence of transmission coefficients")')
     +  yesno(flagfullhf)
      write(*,'(" eciscompound        ",a1,"     flageciscomp flag for",
     +  " compound nucleus calculation by ECIS")') yesno(flageciscomp)
      write(*,'(" cpang               ",a1,"     flagcpang    flag for",
     +  " compound angular distribution calculation for incident",
     +  " charged particles")') yesno(flagcpang)
      if (numinc.gt.1.and.enincmin.lt.eurr.and.enincmax.ge.eurr) then
        write(*,'(" urr               ",f7.3," eurr         off-set",
     +    " incident energy for URR calculation")') eurr
      else
        write(*,'(" urr                 ",a1,"     flagurr    ",
     +    "  flag for URR calculation")') yesno(flagurr)
      endif
      write(*,'(" urrnjoy             ",a1,"     flagurrnjoy  flag for",
     +  " normalization of URR parameters with NJOY method")')
     +  yesno(flagurrnjoy)
      write(*,'(" lurr              ",i3,"     lurr       ",
     +  "  maximal orbital angular momentum for URR")') lurr
c
c 5. Gamma emission
c
c gammax      : number of l-values for gamma multipolarity
c strength    : model for E1 gamma-ray strength function
c strengthM1  : model for M1 gamma-ray strength function
c flagelectron: flag for application of electron conversion coefficient
c flagracap   : flag for radiative capture model
c ldmodelracap: level density model for direct radiative capture
c
      write(*,'(" #"/" # Gamma emission"/" #")')
      write(*,'(" gammax             ",i2,"     gammax",
     +  "       number of l-values for gamma multipolarity")') gammax
      write(*,'(" strength           ",i2,"     strength",
     +  "     model for E1 gamma-ray strength function")') strength
      write(*,'(" strengthM1         ",i2,"     strengthM1",
     +  "   model for M1 gamma-ray strength function")') strengthM1
      write(*,'(" electronconv        ",a1,"     flagelectron",
     +  " flag for application of electron conversion coefficient")')
     +  yesno(flagelectron)
      write(*,'(" racap               ",a1,"     flagracap   ",
     +  " flag for radiative capture model")') yesno(flagracap)
      write(*,'(" ldmodelracap       ",i2,"     ldmodelracap",
     +  " level density model for direct radiative capture")')
     +  ldmodelracap
c
c 6. Pre-equilibrium
c
c epreeq      : on-set incident energy for preequilibrium calculation
c flagpreeq   : flag for pre-equilibrium calculation
c preeqmode   : designator for pre-equilibrium model
c flagmulpre  : flag for multiple pre-equilibrium calculation
c mpreeqmode  : designator for multiple pre-equilibrium model
c emulpre     : on-set incident energy for multiple preequilibrium
c phmodel     : particle-hole state density model
c pairmodel   : model for preequilibrium pairing energy
c pespinmodel : model for pre-equilibrium spin distribution or compound
c               spin distribution for pre-equilibrium cross section
c flaggiant0  : flag for collective contribution from giant resonances
c flagsurface : flag for surface effects in exciton model
c flagpecomp  : flag for Kalbach complex particle emission model
c flag2comp   : flag for two-component pre-equilibrium model
c flagecisdwba: flag for new ECIS calculation for DWBA for MSD
c flagonestep : flag for continuum one-step direct only
c
      write(*,'(" #"/" # Pre-equilibrium"/" #")')
      if (numinc.gt.1.and.enincmin.lt.epreeq.and.enincmax.ge.epreeq)
     +  then
        write(*,'(" preequilibrium    ",f7.3," epreeq       on-set",
     +    " incident energy for preequilibrium calculation")') epreeq
      else
        write(*,'(" preequilibrium      ",a1,"     flagpreeq  ",
     +    "  flag for pre-equilibrium calculation")') yesno(flagpreeq)
      endif
      write(*,'(" preeqmode          ",i2,"     preeqmode",
     +  "    designator for pre-equilibrium model")') preeqmode
      if (numinc.gt.1.and.enincmin.lt.emulpre.and.enincmax.ge.emulpre)
     +  then
        write(*,'(" multipreeq        ",f7.3," emulpre      on-set",
     +    " incident energy for multiple preequilibrium")') emulpre
      else
        write(*,'(" multipreeq          ",a1,"     flagmulpre   ",
     +    "flag for multiple pre-equilibrium calculation")')
     +    yesno(flagmulpre)
      endif
      write(*,'(" mpreeqmode         ",i2,"     mpreeqmode",
     +  "   designator for multiple pre-equilibrium model")')
     +  mpreeqmode
      write(*,'(" phmodel            ",i2,"     phmodel    ",
     +  "  particle-hole state density model")') phmodel
      write(*,'(" pairmodel          ",i2,"     pairmodel",
     +  "    designator for pre-equilibrium pairing model")') pairmodel
      write(*,'(" preeqspin           ",i1,"     pespinmodel",
     +  "  model for pre-equilibrium spin distribution")') pespinmodel
      write(*,'(" giantresonance      ",a1,"     flaggiant    ",
     +  "flag for collective contribution from giant resonances")')
     +  yesno(flaggiant0)
      write(*,'(" preeqsurface        ",a1,"     flagsurface  flag",
     +  " for surface effects in exciton model")') yesno(flagsurface)
      write(*,'(" preeqcomplex        ",a1,"     flagpecomp",
     +  "   flag for Kalbach complex particle emission model")')
     +  yesno(flagpecomp)
      write(*,'(" twocomponent        ",a1,"     flag2comp",
     +  "    flag for two-component pre-equilibrium model")')
     +  yesno(flag2comp)
      write(*,'(" ecisdwba            ",a1,"     flagecisdwba",
     +  " flag for new ECIS calculation for DWBA for MSD")')
     +  yesno(flagecisdwba)
      write(*,'(" onestep             ",a1,"     flagonestep  flag",
     +  " for continuum one-step direct only")') yesno(flagonestep)
c
c 7. Level densities
c
c ldmodelall  : level density model for all nuclides
c spincutmodel: model for spin cutoff factor for ground state
c shellmodel  : model for shell correction energies
c kvibmodel   : model for vibrational enhancement
c flagasys    : flag for all level density parameters a from
c               systematics
c flagparity  : flag for non-equal parity distribution
c flagcolall  : flag for collective enhancement of level density
c               for all nuclides
c flagctmglob : flag for global CTM model (no discrete level info)
c flaggshell  : flag for energy dependence of single particle level
c               density parameter g
c
      write(*,'(" #"/" # Level densities"/" #")')
      write(*,'(" ldmodel            ",i2,"     ldmodelall ",
     +  "  level density model")') ldmodelall
      write(*,'(" shellmodel         ",i2,"     shellmodel  ",
     +  " model for shell correction energies")') shellmodel
      write(*,'(" kvibmodel          ",i2,"     kvibmodel   ",
     +  " model for vibrational enhancement")') kvibmodel
      write(*,'(" spincutmodel       ",i2,"     spincutmodel",
     +  " model for spin cutoff factor for ground state")') spincutmodel
      write(*,'(" asys                ",a1,"     flagasys     flag",
     +  " for all level density parameters a from systematics")')
     +  yesno(flagasys)
      write(*,'(" parity              ",a1,"     flagparity  ",
     +  " flag for non-equal parity distribution")') yesno(flagparity)
      write(*,'(" colenhance          ",a1,"     flagcolall  ",
     +  " flag for collective enahncement of level density",
     +  " for all nuclides")') yesno(flagcolall)
      write(*,'(" ctmglobal           ",a1,"     flagctmglob ",
     +  " flag for global CTM model (no discrete level info)")')
     +  yesno(flagctmglob)
      write(*,'(" gshell              ",a1,"     flaggshell  ",
     +  " flag for energy dependence of single particle",
     +  " level density parameter g")') yesno(flaggshell)
c
c 8. Fission
c
c flagfission: flag for fission
c fismodel   : fission model
c fismodelalt: alternative fission model for default barriers
c flaghbstate: flag for head band states in fission
c flagclass2 : flag for class2 states in fission
c flagmassdis: flag for calculation of fission fragment mass yields
c flagffevap : flag for calculation of particle evaporation from
c              fission fragment mass yields
c flagfisfeed: flag for output of fission per excitation bin
c fymodel    : fission yield model, 1: Brosa 2: GEF
c
      write(*,'(" #"/" # Fission"/" #")')
      write(*,'(" fission             ",a1,"     flagfission",
     +  "  flag for fission")') yesno(flagfission)
      write(*,'(" fismodel           ",i2,"     fismodel   ",
     +  "  fission model")') fismodel
      write(*,'(" fismodelalt        ",i2,"     fismodelalt ",
     +  " alternative fission model for default barriers")') fismodelalt
      write(*,'(" hbstate             ",a1,"     flaghbstate ",
     +  " flag for head band states in fission")') yesno(flaghbstate)
      write(*,'(" class2              ",a1,"     flagclass2 ",
     +  "  flag for class2 states in fission")') yesno(flagclass2)
      write(*,'(" massdis             ",a1,"     flagmassdis",
     +  "  flag for calculation of fission fragment mass yields")')
     +  yesno(flagmassdis)
      write(*,'(" ffevaporation       ",a1,"     flagffevap  ",
     +  " flag for calculation of particle evaporation",
     +  " from fission fragment mass yields")') yesno(flagffevap)
      write(*,'(" fisfeed             ",a1,"     flagfisfeed",
     +  "  flag for output of fission per excitation bin")')
     +  yesno(flagfisfeed)
      write(*,'(" fymodel             ",i1,"     fymodel    ",
     +  "  fission yield model, 1: Brosa 2: GEF")') fymodel
c
c 9. Output
c
c flagmain    : flag for main output
c flagbasic   : flag for output of basic information and results
c flagpop     : flag for output of population
c flagcheck   : flag for output of numerical checks
c flaglevels  : flag for output of discrete level information
c flagdensity : flag for output of level densities
c flagoutomp  : flag for output of optical model parameters
c flagdirect  : flag for output of direct reaction results
c flaginverse : flag for output of transmission coefficients and
c               inverse reaction cross sections
c flagtransen : flag for output of transmission coefficients per energy
c flagoutecis : flag for output of ECIS results
c flagurr     : flag for output of unresolved resonance parameters
c flaggamma   : flag for output of gamma-ray information
c flagpeout   : flag for output of pre-equilibrium results
c flagfisout  : flag for output of fission information
c flagdisc    : flag for output of discrete state cross sections
c flagspec    : flag for output of spectra
c eadd        : on-set incident energy for addition of discrete states
c               to spectra
c eaddel      : on-set incident energy for addition of elastic peak
c               to spectra
c flagadd     : flag for addition of discrete states to spectra
c flagaddel   : flag for addition of elastic peak to spectra
c flagang     : flag for output of angular distributions
c flaglegendre: flag for output of Legendre coefficients
c ddxmode     : mode for double-differential cross sections: 0: None,
c               1: Angular distributions, 2: Spectra per angle, 3: Both
c flagoutdwba : flag for output of DWBA cross sections for MSD
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flagexc     : flag for output of excitation functions
c flagcompo   : flag for output of cross section components
c flagendf    : flag for information for ENDF-6 file
c flagendfdet : flag for detailed ENDF-6 information per channel
c flagpartable: flag for output of model parameters on separate file
c
      write(*,'(" #"/" # Output"/" #")')
      write(*,'(" outmain             ",a1,"     flagmain     ",
     +  "flag for main output")') yesno(flagmain)
      write(*,'(" outbasic            ",a1,"     flagbasic    flag for",
     +  " output of basic information and results")') yesno(flagbasic)
      write(*,'(" outpopulation       ",a1,"     flagpop",
     +  "      flag for output of population")') yesno(flagpop)
      write(*,'(" outcheck            ",a1,"     flagcheck",
     +  "    flag for output of numerical checks")') yesno(flagcheck)
      write(*,'(" outlevels           ",a1,"     flaglevels   flag",
     +  " for output of discrete level information")') yesno(flaglevels)
      write(*,'(" outdensity          ",a1,"     flagdensity",
     +  "  flag for output of level densities")') yesno(flagdensity)
      write(*,'(" outomp              ",a1,"     flagoutomp   flag",
     +  " for output of optical model parameters")') yesno(flagoutomp)
      write(*,'(" outdirect           ",a1,"     flagdirect   flag",
     +  " for output of direct reaction results")') yesno(flagdirect)
      write(*,'(" outinverse          ",a1,"     flaginverse",
     +  "  flag for output of transmission coefficients",
     +  " and inverse reaction cross sections")') yesno(flaginverse)
      write(*,'(" outtransenergy      ",a1,"     flagtransen",
     +  "  flag for output of transmission coefficients per energy")')
     +  yesno(flagtransen)
      write(*,'(" outecis             ",a1,"     flagoutecis ",
     +  " flag for output of ECIS results")') yesno(flagoutecis)
      write(*,'(" outgamma            ",a1,"     flaggamma    flag",
     +  " for output of gamma-ray information")') yesno(flaggamma)
      write(*,'(" outpreequilibrium   ",a1,"     flagpeout",
     +  "    flag for output of pre-equilibrium results ")')
     +  yesno(flagpeout)
      write(*,'(" outfission          ",a1,"     flagfisout   flag",
     +  " for output of fission information")') yesno(flagfisout)
      write(*,'(" outdiscrete         ",a1,"     flagdisc",
     +  "     flag for output of discrete state cross sections")')
     +  yesno(flagdisc)
      write(*,'(" outspectra          ",a1,"     flagspec",
     +  "     flag for output of double-differential cross",
     +  " sections")') yesno(flagspec)
      if (numinc.gt.1.and.enincmin.lt.eadd.and.enincmax.ge.eadd)
     +  then
        write(*,'(" adddiscrete       ",f7.3," eadd         on-set ",
     +    "incident energy for addition of discrete peaks to spectra")')
     +    eadd
      else
        write(*,'(" addiscrete          ",a1,"     flagadd      ",
     +    "flag for addition of discrete states to spectra")')
     +    yesno(flagadd)
      endif
      if (numinc.gt.1.and.enincmin.lt.eaddel.and.enincmax.ge.eaddel)
     +  then
        write(*,'(" addelastic        ",f7.3," eaddel       on-set",
     +    " incident energy addition of elastic peak to spectra")')
     +    eaddel
      else
        write(*,'(" addelastic          ",a1,"     flagaddel    ",
     +    "flag for addition of elastic peak to spectra")')
     +    yesno(flagaddel)
      endif
      write(*,'(" outangle            ",a1,"     flagang",
     +  "      flag for output of angular distributions")')
     +  yesno(flagang)
      write(*,'(" outlegendre         ",a1,"     flaglegendre flag",
     +  " for output of Legendre coefficients")') yesno(flaglegendre)
      write(*,'(" ddxmode             ",i1,"     ddxmode      ",
     +  "mode for double-differential cross sections")') ddxmode
      write(*,'(" outdwba             ",a1,"     flagoutdwba  ",
     +  "flag for output of DWBA cross sections for MSD")')
     +  yesno(flagoutdwba)
      write(*,'(" outgamdis           ",a1,"     flaggamdis  ",
     +  " flag for output of discrete gamma-ray intensities")')
     +  yesno(flaggamdis)
      write(*,'(" outexcitation       ",a1,"     flagexc",
     +  "      flag for output of excitation functions")')
     +  yesno(flagexc)
      write(*,'(" components          ",a1,"     flagcompo",
     +  "    flag for output of cross section components")')
     +  yesno(flagcompo)
      write(*,'(" endf                ",a1,"     flagendf",
     +  "     flag for information for ENDF-6 file")') yesno(flagendf)
      write(*,'(" endfdetail          ",a1,"     flagendfdet",
     +  "  flag for detailed ENDF-6 information per channel")')
     +  yesno(flagendfdet)
      write(*,'(" partable            ",a1,"     flagpartable",
     +  " flag for output of model parameters on separate file")')
     +  yesno(flagpartable)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
