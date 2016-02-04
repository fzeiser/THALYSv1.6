      subroutine checkvalue
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 3, 2013
c | Task  : Check for errors in values
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical lexist
      integer type,i,Zix,Nix,mt,is,l,omptype,nr,nr2,A,irad,lval,igr,
     +        ibar,fax0,n,m,k
      real    tl0,egr0,ggr0,sgr0,epr0,gpr0,tpr0,value,fbar0,fhw0,fR0,
     +        Ea,Eb,Em,D,Ea2,Eb2
c
c All parameters need to fall within certain ranges. These ranges are
c specified in this subroutine and in the manual.
c
c ******************* Check for wrong input variables ******************
c
c 1. Check of values for four main keywords.
c
c ptype0   : type of incident particle
c parsym   : symbol of particle
c numelem  : number of elements
c Starget  : symbol of target nucleus
c nuc      : symbol of nucleus
c Atarget  : mass number of target nucleus
c nummass  : number of masses
c Zinit    : charge number of initial compound nucleus
c Ninit    : neutron number of initial compound nucleus
c enincmin : minimum incident energy
c enincmax : maximum incident energy
c Emaxtalys: maximum acceptable energy for TALYS
c
      do 10 type=0,6
        if (ptype0.eq.parsym(type)) goto 20
   10 continue
      if (ptype0.eq.'0') goto 20
      write(*,'(" TALYS-error: Wrong symbol for projectile: ",a1)')
     +  ptype0
      stop
   20 do 30 i=3,numelem
        if (Starget.eq.nuc(i)) goto 40
   30 continue
      write(*,'(" TALYS-error: Wrong symbol for element: ",a2)') Starget
      stop
   40 if (Atarget.le.5.or.Atarget.gt.nummass) then
        write(*,'(" TALYS-error: 5 < Target mass <= ",i3)') nummass
        stop
      endif
      if (Zinit.le.2) then
        write(*,'(" TALYS-error: Target Z > 2")')
        stop
      endif
      if (Ninit.le.2) then
        write(*,'(" TALYS-error: Target N > 2")')
        stop
      endif
      if (Zinit.gt.numelem) then
        write(*,'(" TALYS-error: combination of projectile ",
     +    "+ element > ",i3)') numelem
        stop
      endif
      if (enincmin.lt.1.e-11.or.enincmax.gt.Emaxtalys) then
        write(*,'(" TALYS-error: 1.e-5 eV <= Incident energy",
     +    " < ",f10.3," MeV")') Emaxtalys
        stop
      endif
c
c 2. Check of values for basic physical and numerical parameters
c
c maxZ,numZ    : maximal number of protons away from the initial
c                compound nucleus
c maxN,numN    : maximal number of neutrons away from the initial
c                compound nucleus
c nbins,numbins: number of continuum excitation energy bins
c numbins      : maximal number of continuum excitation energy bins
c ptype0       : type of incident particle
c segment      : number of segments to divide emission energy grid
c flagastro    : flag for calculation of astrophysics reaction rate
c nlevmax      : maximum number of included discrete levels for target
c nlevmaxres   : maximum number of included discrete levels for residual
c                nucleus
c numlev       : maximum number of included discrete levels
c nlevbin      : number of excited levels for binary nucleus
c nlev         : number of excited levels for nucleus
c Zix          : charge number index for residual nucleus
c numZ         : maximal number of protons away from initial compound
c                nucleus
c Nix          : neutron number index for residual nucleus
c numN         : maximal number of neutrons away from initial compound
c                nucleus
c A            : mass number
c Ainit        : mass number of initial compound nucleus
c massnucleus  : mass of nucleus in amu as read from user input file
c massexcess   : mass excess in MeV as read from user input file
c Ltarget      : excited level of target
c Lisoinp      : user assignment of target isomer number
c isomer       : definition of isomer in seconds
c core         : even-even core for weakcoupling (-1 or 1)
c transpower   : power for transmission coefficient limit
c transeps     : absolute limit for transmission coefficient
c xseps        : limit for cross sections
c popeps       : limit for population cross section per nucleus
c Rfiseps      : ratio for limit for fission cross section per nucleus
c eninclow     : minimal incident energy for nuclear model calculations
c nangle       : number of angles
c numang       : maximum number of angles
c nanglecont   : number of angles for continuum
c numangcont   : maximum number of angles for continuum
c nanglerec    : number of recoil angles
c numangrec    : maximum number of recoil angles
c maxenrec     : number of recoil energies
c numenrec     : maximum number of recoil energies
c maxchannel   : maximal number of outgoing particles in individual
c                channel description (e.g. this is 3 for (n,2np))
c massmodel    : model for theoretical nuclear mass
c disctable    : table with discrete levels
c astroT9      : temperature, in 10^9 K, for Maxwellian average
c astroE       : energy, in MeV, for Maxwellian average
c kT           : energy kT expressed in MeV corresponding to a
c                temperature T9=1
c flagprod     : flag for isotope production
c Ebeam        : incident energy in MeV for isotope production
c Eback        : lower end of energy range in MeV for isotope
c                production
c eninc        : incident energy in MeV
c radiounit    : unit for radioactivity: Bq, kBq, MBq, Gbq,
c                mCi, Ci or kCi
c yieldunit    : unit for isotope yield: num (number),
c                mug (micro-gram), mg, g, or kg
c Ibeam        : beam current in mA for isotope production
c Area         : target area in cm^2 for isotope production
c Tirrad       : irradiation time per unit
c unitTirrad   : irradiation time unit (y,d,h,m,s)
c Tcool        : cooling time per unit
c unitTcool    : cooling time unit (y,d,h,m,s)
c rhotarget    : target material density
c
      if (maxZ.lt.0.or.maxZ.gt.numZ-2) then
        write(*,'(" TALYS-error: 0 <= maxZ <=",i3)') numZ-2
        stop
      endif
      if (maxN.lt.0.or.maxN.gt.numN-2) then
        write(*,'(" TALYS-error: 0 <= maxN <=",i3)') numN-2
        stop
      endif
      if (nbins0.ne.0.and.(nbins0.lt.2.or.nbins0.gt.numbins)) then
        write(*,'(" TALYS-error: 2 <= bins <=",i3)') numbins
        stop
      endif
      if (segment.le.0.or.segment.gt.4) then
        write(*,'(" TALYS-error: 1 <= segment <= 4")')
        stop
      else
        if (segment.gt.1.and.enincmax.gt.100.) then
          write(*,'(" TALYS-error: segment = 1",
     +      " for incident energy of ",f7.3," MeV")') enincmax
          stop
        endif
        if (segment.gt.2.and.enincmax.gt.40.) then
          write(*,'(" TALYS-error: 1 <= segment = 2",
     +      " for incident energy of ",f7.3," MeV")') enincmax
          stop
        endif
        if (segment.gt.3.and.enincmax.gt.20.) then
          write(*,'(" TALYS-error: 1 <= segment = 3",
     +      " for incident energy of ",f7.3," MeV")') enincmax
          stop
        endif
        if (segment.gt.1.and.flagastro) then
          write(*,'(" TALYS-error: segment = 1",
     +      " for astrophysical calculations")')
          stop
        endif
      endif
      if (nlevmax.lt.0.or.nlevmax.gt.numlev) then
        write(*,'(" TALYS-error: 0 <= maxlevelstar <=",i3)') numlev
        stop
      endif
      if (nlevmaxres.lt.0.or.nlevmaxres.gt.numlev) then
        write(*,'(" TALYS-error: 0 <= maxlevelsres <=",i3)') numlev
        stop
      endif
      do 50 type=0,6
        if (nlevbin(type).lt.0.or.nlevbin(type).gt.numlev) then
          write(*,'(" TALYS-error: 0 <= maxlevelsbin <=",i3)') numlev
          stop
        endif
   50 continue
      do 60 Zix=0,numZ
        do 70 Nix=0,numN
          if (nlev(Zix,Nix).lt.0.or.nlev(Zix,Nix).gt.numlev) then
            write(*,'(" TALYS-error: 0 <= nlevels <=",i3)') numlev
            stop
          endif
          A=Ainit-Zix-Nix
          if (massnucleus(Zix,Nix).ne.0.and.
     +      (massnucleus(Zix,Nix).lt.real(A)-1..or.
     +      massnucleus(Zix,Nix).gt.real(A)+1.)) then
            write(*,'(" TALYS-error: A-1. <= massnucleus <= A+1.")')
            stop
          endif
          if (massexcess(Zix,Nix).ne.0.and.
     +      (massexcess(Zix,Nix).lt.-600..or.
     +      massexcess(Zix,Nix).gt.600.)) then
            write(*,'(" TALYS-error: -600. <= massexcess <= 600.")')
            stop
          endif
   70   continue
   60 continue
      if (Ltarget.lt.0.or.Ltarget.gt.numlev) then
        write(*,'(" TALYS-error: 0 <= Ltarget <=",i3)') numlev
        stop
      endif
      if (Lisoinp.ne.-1.and.(Lisoinp.lt.0.or.Lisoinp.gt.9)) then
        write(*,'(" TALYS-error: 0 <= Liso <= 9")')
        stop
      endif
      if (isomer.lt.0..or.isomer.gt.1.e38) then
        write(*,'(" TALYS-error: 0. <= isomer <= 1.e38")')
        stop
      endif
      if (core.ne.-1.and.core.ne.1) then
        write(*,'(" TALYS-error: core = -1 or 1")')
        stop
      endif
      if (transpower.lt.2.or.transpower.gt.20) then
        write(*,'(" TALYS-error: 2 <= transpower <= 20")')
        stop
      endif
      if (transeps.lt.0..or.transeps.gt.1.) then
        write(*,'(" TALYS-error: 0. <= transeps <= 1.")')
        stop
      endif
      if (xseps.lt.0..or.xseps.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= xseps <= 1000.")')
        stop
      endif
      if (popeps.lt.0..or.popeps.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= popeps <= 1000.")')
        stop
      endif
      if (Rfiseps.lt.0..or.Rfiseps.gt.1.) then
        write(*,'(" TALYS-error: 0. <= Rfiseps <= 1.")')
        stop
      endif
      if (eninclow.ne.0..and.(eninclow.lt.1.e-6.or.eninclow.gt.1.))
     +  then
        write(*,'(" TALYS-error: 1.e-6 <= Elow <= 1.")')
        stop
      endif
      if (nangle.lt.1..or.nangle.gt.numang) then
        write(*,'(" TALYS-error: 1 <= angles <=",i3)') numang
        stop
      endif
      if (nanglecont.lt.1..or.nanglecont.gt.numangcont) then
        write(*,'(" TALYS-error: 1 <= anglescont <=",i3)') numangcont
        stop
      endif
      if (nanglerec.lt.1..or.nanglerec.gt.numangrec) then
        write(*,'(" TALYS-error: 1 <= anglesrec <=",i3)') numangrec
        stop
      endif
      if (maxenrec.lt.1..or.maxenrec.gt.numenrec) then
        write(*,'(" TALYS-error: 1 <= maxenrec <=",i3)') numenrec
        stop
      endif
      if (maxchannel.lt.1.or.maxchannel.gt.8) then
        write(*,'(" TALYS-error: 1 <= maxchannel <= 8")')
        stop
      endif
      if (massmodel.lt.0.or.massmodel.gt.3) then
        write(*,'(" TALYS-error: 0 <= massmodel <= 3")')
        stop
      endif
      if (disctable.lt.1.or.disctable.gt.3) then
        write(*,'(" TALYS-error: 1 <= disctable <= 3")')
        stop
      endif
      if (astroT9.ne.0..and.(astroT9.lt.0.0001.or.astroT9.gt.10.)) then
        write(*,'(" TALYS-error: 0.0001 <= astroT <= 10")')
        stop
      endif
      if (astroE.ne.0..and.(astroE.lt.0.00001.or.astroE.gt.1.)) then
        write(*,'(" TALYS-error: 0.00001 <= astroE <= 1")')
        stop
      endif
      if (astroE.ne.0..and.astroT9.ne.0.) then
        write(*,'(" TALYS-error: Only astroE OR astroT can be given")')
        stop
      endif
      if (astroE.ne.0.) astroT9=astroE/kT
      if (astroT9.ne.0.) astroE=astroT9*kT
      if (flagprod) then
        if (Ebeam.eq.-1.) then
          write(*,'(" TALYS-error: accelerator energy Ebeam must be ",
     +      "given for isotope production (production y)")')
          stop
        endif
        if (Ebeam.le.0..or.Ebeam.gt.Emaxtalys) then
          write(*,'(" TALYS-error: 0 < Ebeam < ",f8.3," MeV")')
     +      Emaxtalys
          stop
        endif
        if (Eback.eq.-1.) then
          Eback=max(Ebeam-5.,0.1)
        else
          if (Eback.le.0..or.Eback.gt.Emaxtalys) then
            write(*,'(" TALYS-error: 0 < Eback < ",f8.3," MeV")')
     +        Emaxtalys
            stop
          endif
        endif
        if (Eback.ge.Ebeam) then
          write(*,'(" TALYS-error: Ebeam must be larger than Eback")')
          stop
        endif
        if (Ebeam.gt.enincmax+1.e-4) then
          write(*,'(" TALYS-error: Ebeam is not in the energy range ",
     +      "with TALYS results, Ebeam=",f10.5," Ein(max)=",f10.5,
     +      ". Rerun with wider energy grid")') Ebeam,enincmax
          stop
        endif
        if (Eback.lt.eninc(1)-1.e-4) then
          write(*,'(" TALYS-error: Eback is not in the energy range ",
     +      "with TALYS results, Eback=",f10.5," Ein(1)=",f10.5,
     +      ". Rerun with wider energy grid")') Eback,eninc(1)
          stop
        endif
        if (radiounit.ne.'bq'.and.radiounit.ne.'kbq'.and.
     +    radiounit.ne.'mbq'.and.radiounit.ne.'gbq'.and.
     +    radiounit.ne.'mci'.and.radiounit.ne.'ci'.and.
     +    radiounit.ne.'kci') then
          write(*,'(" TALYS-error: radiounit should be equal to ",
     +      "Bq, kBq, MBq, Gbq, mCi, Ci or kCi")')
          stop
        endif
        if (yieldunit.ne.'num'.and.yieldunit.ne.'mug'.and.
     +    yieldunit.ne.'mg'.and.yieldunit.ne.'g'.and.yieldunit.ne.'kg')
     +    then
          write(*,'(" TALYS-error: yieldunit should be equal to ",
     +      "num (number), mug (micro-gram), mg, g, or kg")')
          stop
        endif
        if (Ibeam.le.0..or.Ibeam.gt.10000.) then
          write(*,'(" TALYS-error: 0 <= Ibeam < 1000 mA")')
          stop
        endif
        if (Area.le.0..or.Area.gt.10000.) then
          write(*,'(" TALYS-error: 0 <= Area < 10000 cm^2")')
          stop
        endif
        do 80 k=1,5
          if (Tirrad(k).lt.0.or.Tirrad(k).ge.1000000) then
            write(*,'(" TALYS-error: 0 <= Tirrad < 1.e6")')
            stop
          endif
          if (Tcool(k).lt.0.or.Tcool(k).ge.1000000) then
            write(*,'(" TALYS-error: 0 <= Tcool < 1.e6")')
            stop
          endif
   80   continue
        do 90 k=1,5
          if (unitTirrad(k).ne.' '.and.unitTirrad(k).ne.'y'.and.
     +      unitTirrad(k).ne.'d'.and.unitTirrad(k).ne.'h'.and.
     +      unitTirrad(k).ne.'m'.and.unitTirrad(k).ne.'s') then
            write(*,'(" TALYS-error: wrong unit for Tirrad= ",i9)')
     +        Tirrad(k)
            stop
          endif
          if (unitTcool(k).ne.' '.and.unitTcool(k).ne.'y'.and.
     +      unitTcool(k).ne.'d'.and.unitTcool(k).ne.'h'.and.
     +      unitTcool(k).ne.'m'.and.unitTcool(k).ne.'s') then
            write(*,'(" TALYS-error: wrong unit for Tcool= ",i9)')
     +        Tcool(k)
            stop
          endif
   90   continue
        if (rhotarget.ne.-1..and.(rhotarget.le.0..or.rhotarget.gt.100.))
     +    then
          write(*,'(" TALYS-error: 0 < rhotarget <= 100.")')
          stop
        endif
      endif
c
c 3. Check of values of optical model
c
c numNph     : maximal number of neutrons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c numZph     : maximal number of protons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c optmod     : file with optical model parameters
c optmodfileN: optical model parameter file for neutrons
c optmodfileP: optical model parameter file for protons
c radialfile : radial matter density file
c
      do 110 Zix=0,numZph
        do 120 Nix=0,numNph
          do 130 type=1,6
            if (optmod(Zix,Nix,type)(1:1).eq.' ') goto 130
            inquire (file=optmod(Zix,Nix,type),exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Non-existent optical model",
     +          " file: ",a72)') optmod(Zix,Nix,type)
              stop
            endif
  130     continue
  120   continue
        if (optmodfileN(Zix)(1:1).ne.' ') then
          inquire (file=optmodfileN(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent optical model",
     +        " file: ",a72)') optmodfileN(Zix)
            stop
          endif
        endif
        if (optmodfileP(Zix)(1:1).ne.' ') then
          inquire (file=optmodfileP(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent optical model",
     +        " file: ",a72)') optmodfileP(Zix)
            stop
          endif
        endif
        if (radialfile(Zix)(1:1).ne.' ') then
          inquire (file=radialfile(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent radial file: ",a72)')
     +        radialfile(Zix)
            stop
          endif
        endif
c
c Check other parameter input files
c
c levelfile    : discrete level file
c deformfile   : deformation parameter file
c hbtransfile  : file with head band transition states
c clas2file    : file with class 2 transition states
c ompenergyfile: file with energies for OMP calculation (ENDF files
c                only)
c rescuefile   : file with incident energy dependent adjustment factors
c grescue      : global multiplication factor for incident energy
c                dependent adjustment factors
c alphaomp     : alpha optical model (1=normal, 2= McFadden-Satchler,
c                3-5= folding potential)
c deuteronomp  : deuteron optical model (1=normal, 2=Daehnick,
c                3=Bojowald, 4=Han-Shi-Shen, 5=An-Cai)
c radialmodel  : model for radial matter densities (JLM OMP only)
c
        if (levelfile(Zix)(1:1).ne.' ') then
          inquire (file=levelfile(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent level file: ",a72)')
     +        levelfile(Zix)
            stop
          endif
        endif
        if (deformfile(Zix)(1:1).ne.' ') then
          inquire (file=deformfile(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent deformation ",
     +        "parameter file: ",a72)') deformfile(Zix)
            stop
          endif
        endif
        do 140 Nix=0,numN
          if (hbtransfile(Zix,Nix)(1:1).ne.' ') then
            inquire (file=hbtransfile(Zix,Nix),exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Non-existent head band ",
     +          "transition state file: ",a72)')
     +          hbtransfile(Zix,Nix)
              stop
            endif
          endif
          if (clas2file(Zix,Nix)(1:1).ne.' ') then
            inquire (file=clas2file(Zix,Nix),exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Non-existent class 2 ",
     +          "transition state file: ",a72)') clas2file(Zix,Nix)
              stop
            endif
          endif
  140   continue
  110 continue
      if (ompenergyfile(1:1).ne.' ') then
        inquire (file=ompenergyfile,exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: Non-existent ompenergyfile: ",a72)')
     +      ompenergyfile
          stop
        endif
      endif
      do 150 mt=1,nummt
        do 155 is=-1,1
          if (rescuefile(mt,is)(1:1).ne.' ') then
            inquire (file=rescuefile(mt,is),exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Non-existent rescue file: ",
     +          a72)') rescuefile(mt,is)
              stop
            endif
          endif
          if (grescue(mt,is).lt.0.001.or.grescue(mt,is).gt.1000.) then
            write(*,'(" TALYS-error: 0.001 <= grescue <= 1000.")')
            stop
          endif
  155   continue
  150 continue
      if (alphaomp.lt.1.or.alphaomp.gt.5) then
        write(*,'(" TALYS-error: 1 <= alphaomp <= 5")')
        stop
      endif
      if (deuteronomp.lt.1.or.deuteronomp.gt.5) then
        write(*,'(" TALYS-error: 1 <= deuteronomp <= 5")')
        stop
      endif
      if (radialmodel.lt.1.or.radialmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= radialmodel <= 2")')
        stop
      endif
c
c Check adjustable OMP parameters
c
c v1adjust.. : adjustable factors for OMP (default 1.)
c jlmmode    : option for JLM imaginary potential normalization
c soswitch   : switch for deformed spin-orbit calculation and
c              sequential iterations in ECIS
c ompadjustN : number of energy ranges for local OMP adjustment
c ompadjustE1: start energy of local OMP adjustment
c ompadjustE2: end energy of local OMP adjustment
c ompadjustD : depth of local OMP adjustment
c ompadjusts : variance of local OMP adjustment
c tladjust   : adjustable factor for Tlj (default 1.)
c aradialcor : adjustable parameter for shape of DF alpha potential
c adepthcor  : adjustable parameter for depth of DF alpha potential
c Ejoin      : joining energy for high energy OMP
c Vinfadjust : adjustable factor for high energy limit of
c              real central potential
c
      do 160 type=1,6
        if (v1adjust(type).lt.0.1.or.v1adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= v1adjust <= 10.")')
          stop
        endif
        if (v2adjust(type).lt.0.1.or.v2adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= v2adjust <= 10.")')
          stop
        endif
        if (v3adjust(type).lt.0.1.or.v3adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= v3adjust <= 10.")')
          stop
        endif
        if (v4adjust(type).lt.0.1.or.v4adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= v4adjust <= 10.")')
          stop
        endif
        if (rvadjust(type).lt.0.1.or.rvadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= rvadjust <= 10.")')
          stop
        endif
        if (avadjust(type).lt.0.1.or.avadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= avadjust <= 10.")')
          stop
        endif
        if (rwadjust(type).lt.0.1.or.rwadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= rwadjust <= 10.")')
          stop
        endif
        if (awadjust(type).lt.0.1.or.awadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= awadjust <= 10.")')
          stop
        endif
        if (w1adjust(type).lt.0.1.or.w1adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= w1adjust <= 10.")')
          stop
        endif
        if (w2adjust(type).lt.0.1.or.w2adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= w2adjust <= 10.")')
          stop
        endif
        if (w3adjust(type).lt.0.1.or.w3adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= w3adjust <= 10.")')
          stop
        endif
        if (w4adjust(type).lt.0.1.or.w4adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= w4adjust <= 10.")')
          stop
        endif
        if (rvdadjust(type).lt.0.1.or.rvdadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= rvdadjust <= 10.")')
          stop
        endif
        if (avdadjust(type).lt.0.1.or.avdadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= avdadjust <= 10.")')
          stop
        endif
        if (d1adjust(type).lt.0.1.or.d1adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= d1adjust <= 10.")')
          stop
        endif
        if (d2adjust(type).lt.0.1.or.d2adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= d2adjust <= 10.")')
          stop
        endif
        if (d3adjust(type).lt.0.1.or.d3adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= d3adjust <= 10.")')
          stop
        endif
        if (rwdadjust(type).lt.0.1.or.rwdadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= rwdadjust <= 10.")')
          stop
        endif
        if (awdadjust(type).lt.0.1.or.awdadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= awdadjust <= 10.")')
          stop
        endif
        if (vso1adjust(type).lt.0.1.or.vso1adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= vso1adjust <= 10.")')
          stop
        endif
        if (vso2adjust(type).lt.0.1.or.vso2adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= vso2adjust <= 10.")')
          stop
        endif
        if (wso1adjust(type).lt.0.1.or.wso1adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= wso1adjust <= 10.")')
          stop
        endif
        if (wso2adjust(type).lt.0.1.or.wso2adjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= wso2adjust <= 10.")')
          stop
        endif
        if (rvsoadjust(type).lt.0.1.or.rvsoadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= rvsoadjust <= 10.")')
          stop
        endif
        if (avsoadjust(type).lt.0.1.or.avsoadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= avsoadjust <= 10.")')
          stop
        endif
        if (rwsoadjust(type).lt.0.1.or.rwsoadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= rwsoadjust <= 10.")')
          stop
        endif
        if (awsoadjust(type).lt.0.1.or.awsoadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= awsoadjust <= 10.")')
          stop
        endif
        if (rcadjust(type).lt.0.1.or.rcadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.1 <= rcadjust <= 10.")')
          stop
        endif
        do 170 omptype=1,numompadj
          do 180 nr=1,ompadjustN(type,omptype)
            if (ompadjustE1(type,omptype,nr).lt.0.or.
     +        ompadjustE1(type,omptype,nr).gt.Emaxtalys) then
              write(*,'(" TALYS-error: 0. <= ompadjustE1 <= ",f10.3)')
     +          Emaxtalys
              stop
            endif
            if (ompadjustE2(type,omptype,nr).lt.0.or.
     +        ompadjustE2(type,omptype,nr).gt.Emaxtalys) then
              write(*,'(" TALYS-error: 0. <= ompadjustE2 <= ",f10.3)')
     +          Emaxtalys
              stop
            endif
            if (ompadjustE2(type,omptype,nr).lt.
     +        ompadjustE1(type,omptype,nr)) then
              write(*,'(" TALYS-error: ompadjustE1 <= ompadjustE2 ")')
              stop
            endif
            do 190 nr2=1,ompadjustN(type,omptype)
              if (nr.eq.nr2) goto 190
              if (ompadjustE1(type,omptype,nr).gt.
     +          ompadjustE1(type,omptype,nr2).and.
     +          ompadjustE1(type,omptype,nr).lt.
     +          ompadjustE2(type,omptype,nr2)) then
                write(*,'(" TALYS-error: ompadjustE1 and ",
     +            "ompadjustE2 overlapping")')
                stop
              endif
  190       continue
            if (ompadjustD(type,omptype,nr).lt.-100.or.
     +        ompadjustD(type,omptype,nr).gt.100.) then
              write(*,'(" TALYS-error: -100. <= ompadjustD <= 100.")')
              stop
            endif
            if (ompadjusts(type,omptype,nr).lt.0..or.
     +        ompadjusts(type,omptype,nr).gt.100.) then
              write(*,'(" TALYS-error: 0. <= ompadjusts <= 100.")')
              stop
            endif
  180     continue
  170   continue
  160 continue
      do 192 type=-1,6
        do 194 l=0,numl
          tl0=tladjust(type,l)
          if (tl0.lt.0.001.or.tl0.gt.1000.) then
            write(*,'(" TALYS-error: 0.001 <= tladjust <= 1000.")')
            stop
          endif
  194   continue
  192 continue
      if (jlmmode.lt.0.or.jlmmode.gt.3) then
        write(*,'(" TALYS-error: 0 <= jlmmode <= 3")')
        stop
      endif
      if (lvadjust.lt.0.5.or.lvadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.2 <= lvadjust <= 5.")')
        stop
      endif
      if (lwadjust.lt.0.5.or.lwadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lwadjust <= 1.5")')
        stop
      endif
      if (lv1adjust.lt.0.5.or.lv1adjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lv1adjust <= 1.5")')
        stop
      endif
      if (lw1adjust.lt.0.5.or.lw1adjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lw1adjust <= 1.5")')
        stop
      endif
      if (lvsoadjust.lt.0.5.or.lvsoadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lvsoadjust <= 1.5")')
        stop
      endif
      if (lwsoadjust.lt.0.5.or.lwsoadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lwsoadjust <= 1.5")')
        stop
      endif
      if (aradialcor.lt.0.5.or.aradialcor.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= aradialcor <= 1.5")')
        stop
      endif
      if (adepthcor.lt.0.5.or.adepthcor.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= adepthcor <= 1.5")')
        stop
      endif
      if (soswitch.lt.0.1.or.soswitch.gt.10.) then
        write(*,'(" TALYS-error: 0.1 <= soswitch <= 10.")')
        stop
      endif
      do 196 type=1,2
        if (Ejoin(type).le.0..or.Ejoin(type).gt.Emaxtalys) then
          write(*,'(" TALYS-error: 0. < Ejoin <= ",f10.5)')
     +      Emaxtalys
          stop
        endif
        if (Vinfadjust(type).lt.0.01.or.Vinfadjust(type).gt.10.) then
          write(*,'(" TALYS-error: 0.01 <= Vinfadjust <= 10.")')
          stop
        endif
  196 continue
c
c Check direct reaction parameters
c
c maxband   : highest vibrational band added to rotational model
c maxrot    : number of included excited rotational levels
c k0        : index for incident particle
c flaggiant0: flag for collective contribution from giant resonances
c
      if (maxband.lt.0.or.maxband.gt.10) then
        write(*,'(" TALYS-error: 0 <= maxband <= 10")')
        stop
      endif
      if (maxrot.lt.1.or.maxrot.gt.20) then
        write(*,'(" TALYS-error: 1 <= maxrot <= 20")')
        stop
      endif
      if (k0.eq.0.and.flaggiant0) then
        write(*,'(" TALYS-error: No giant resonance sumrules for",
     +    " photonuclear reactions")')
        stop
      endif
c
c 4. Check of values for compound nucleus
c
c ewfc        : off-set incident energy for width fluctuation
c               calculation
c eurr        : off-set incident energy for URR calculation
c flagurr     : flag for output of unresolved resonance parameters
c wmode       : designator for width fluctuation model
c lurr        : maximal orbital angular momentum for URR
c numl        : maximal number of l-values
c flageciscomp: flag for compound nucleus calculation by ECIS
c
      if ((ewfc.lt.0..or.ewfc.gt.20.).and.ewfc.ne.-1.) then
        write(*,'(" TALYS-error: 0. <= ewfc <= 20.")')
        stop
      endif
      if ((eurr.lt.0..or.eurr.gt.20.).and.eurr.ne.-1.) then
        write(*,'(" TALYS-error: 0. <= eurr <= 20.")')
        stop
      endif
      if (wmode.lt.0.or.wmode.gt.3) then
        write(*,'(" TALYS-error: 0 <= wmode <= 3")')
        stop
      endif
      if (k0.eq.0.and.ewfc.gt.0.) then
        write(*,'(" TALYS-error: No width fluctuations for",
     +    " photonuclear reactions")')
        stop
      endif
      if (k0.ne.1.and.(eurr.gt.0..or.flagurr)) then
        write(*,'(" TALYS-error: URR calculation only possible",
     +    " for incident neutrons")')
        stop
      endif
      if (.not.flagcomp.and.eurr.gt.0.) then
        write(*,'(" TALYS-error: URR calculation only possible",
     +    " if compound nucleus model enabled")')
        stop
      endif
      if (lurr.lt.0.or.lurr.gt.numl) then
        write(*,'(" TALYS-error: 0 <= lurr <= ",i3)') numl
        stop
      endif
      if (k0.eq.0..and.flageciscomp) then
        write(*,'(" TALYS-error: No compound calculation by",
     +    " ECIS for incident photons")')
        stop
      endif
      if (enincmax.gt.20..and.flageciscomp) then
        write(*,'(" TALYS-error: No compound calculation by",
     +    " ECIS for E > 20 MeV")')
        stop
      endif
c
c 5. Check of values for gamma emission
c
c gammax      : number of l-values for gamma multipolarity
c strength    : model for E1 gamma-ray strength function
c strengthM1  : model for M1 gamma-ray strength function
c egr,egr0    : energy of GR
c irad        : variable to indicate M(=0) or E(=1) radiation
c lval        : multipolarity
c igr         : giant resonance
c ggr,ggr0    : width of GR
c sgr,sgr0    : strength of GR
c epr         : energy of PR
c gpr         : width of PR
c tpr         : strength of PR
c egradjust.  : adjustable factors for giant resonance parameters
c               (default 1.)
c gamgam      : total radiative width in eV
c D0          : s-wave resonance spacing in eV
c fiso        : correction factor for isospin forbidden transitions
c gnorm       : gamma normalization factor
c etable      : constant to adjust tabulated strength functions
c ftable      : constant to adjust tabulated strength functions
c RprimeU     : potential scattering radius
c flagracap   : flag for radiative capture model
c ldmodelracap: level density model for direct radiative capture
c spectfacexp : experimental spectroscopic factor
c spectfacth  : theoretical spectroscopic factor
c
      if (gammax.lt.1.or.gammax.gt.6) then
        write(*,'(" TALYS-error: 1 <= gammax <= 6")')
        stop
      endif
      if (strength.lt.1.or.strength.gt.5) then
        write(*,'(" TALYS-error: 1 <= strength <= 5")')
        stop
      endif
      if (strength.eq.3.or.strength.eq.4) then
        inquire (file=path(1:lenpath)//'gamma/hfb/z050',exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: Microscopic HFB tables are not ",
     +    " installed: download the full TALYS package from ",
     +    "www.talys.eu")')
          stop
        endif
      endif
      if (strengthM1.lt.1.or.strengthM1.gt.2) then
        write(*,'(" TALYS-error: 1 <= strengthM1 <= 2")')
        stop
      endif
      do 210 Zix=0,numZ
        do 220 Nix=0,numN
          do 230 irad=0,1
            do 240 lval=1,gammax
              do 250 igr=1,2
                egr0=egr(Zix,Nix,irad,lval,igr)
                if (egr0.ne.0..and.(egr0.lt.1..or.egr0.gt.100.)) then
                  write(*,'(" TALYS-error: 1.<= energy of GR <=100.")')
                  stop
                endif
                ggr0=ggr(Zix,Nix,irad,lval,igr)
                if (ggr0.ne.0..and.(ggr0.lt.1..or.ggr0.gt.100.)) then
                  write(*,'(" TALYS-error: 1.<= width of GR <=100.")')
                  stop
                endif
                sgr0=sgr(Zix,Nix,irad,lval,igr)
                if (sgr0.lt.0..or.sgr0.gt.10000.) then
                  write(*,'(" TALYS-error: 0.<= strength of GR",
     +              " <=10000.")')
                  stop
                endif
                epr0=epr(Zix,Nix,irad,lval,igr)
                if (epr0.ne.0..and.(epr0.lt.1..or.epr0.gt.100.)) then
                  write(*,'(" TALYS-error: 1.<= energy of PR <=100.")')
                  stop
                endif
                gpr0=gpr(Zix,Nix,irad,lval,igr)
                if (gpr0.ne.0..and.(gpr0.lt.0.1.or.gpr0.gt.100.)) then
                  write(*,'(" TALYS-error: 0.1<= width of PR <=100.")')
                  stop
                endif
                tpr0=tpr(Zix,Nix,irad,lval,igr)
                if (tpr0.lt.0..or.tpr0.gt.10000.) then
                  write(*,'(" TALYS-error: 0.<= strength of PR",
     +              " <=10000.")')
                  stop
                endif
                egr0=egradjust(Zix,Nix,irad,lval,igr)
                if (egr0.lt.0.1.or.egr0.gt.10.) then
                  write(*,'(" TALYS-error: 0.1 <= egradjust <= 10.")')
                  stop
                endif
                ggr0=ggradjust(Zix,Nix,irad,lval,igr)
                if (ggr0.lt.0.1.or.ggr0.gt.10.) then
                  write(*,'(" TALYS-error: 0.1 <= ggradjust <= 10.")')
                  stop
                endif
                sgr0=sgradjust(Zix,Nix,irad,lval,igr)
                if (sgr0.lt.0.1.or.sgr0.gt.10.) then
                  write(*,'(" TALYS-error: 0.1 <= sgradjust <= 10.")')
                  stop
                endif
                epr0=epradjust(Zix,Nix,irad,lval,igr)
                if (epr0.lt.0.1.or.epr0.gt.10.) then
                  write(*,'(" TALYS-error: 0.1 <= epradjust <= 10.")')
                  stop
                endif
                gpr0=gpradjust(Zix,Nix,irad,lval,igr)
                if (gpr0.lt.0.1.or.gpr0.gt.10.) then
                  write(*,'(" TALYS-error: 0.1 <= gpradjust <= 10.")')
                  stop
                endif
                tpr0=tpradjust(Zix,Nix,irad,lval,igr)
                if (tpr0.lt.0.1.or.tpr0.gt.10.) then
                  write(*,'(" TALYS-error: 0.1 <= spradjust <= 10.")')
                  stop
                endif
  250         continue
  240       continue
  230     continue
          if (gamgam(Zix,Nix).ne.0.) then
            if (gamgam(Zix,Nix).lt.0..or.gamgam(Zix,Nix).gt.10.) then
              write(*,'(" TALYS-error: 0. < gamgam <= 10.")')
              stop
            endif
          endif
          if (D0(Zix,Nix).ne.0.) then
            if (D0(Zix,Nix).lt.1.e-3.or.D0(Zix,Nix).gt.1.e7) then
              write(*,'(" TALYS-error: 1.e-6 <= D0 <= 10000. keV")')
              stop
            endif
          endif
          if (etable(Zix,Nix).lt.-10..or.etable(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: -10. <= etable <= 10.")')
            stop
          endif
          if (ftable(Zix,Nix).lt.0.1.or.ftable(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: 0.1 <= ftable <= 10.")')
            stop
          endif
          do 260 type=0,6
            if (fiso(Zix,Nix,type).lt.0.01.or.fiso(Zix,Nix,type).gt.100)
     +        then
              write(*,'(" TALYS-error: 0.01 <= fiso <= 100.")')
              stop
            endif
  260     continue
          if (gamgamadjust(Zix,Nix).lt.0.1.or.
     +      gamgamadjust(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: 0.1 <= gamgamadjust <= 10.")')
            stop
          endif
  220   continue
  210 continue
      if ((gnorm.le.0..or.gnorm.gt.1000.).and.gnorm.ne.-1.) then
        write(*,'(" TALYS-error: 0. < gnorm <= 100.")')
        stop
      endif
      if (RprimeU.lt.0..or.RprimeU.gt.10.) then
        write(*,'(" TALYS-error: 0. <= Rprime <= 10.")')
        stop
      endif
      if (flagracap.and.k0.eq.0) then
        write(*,'(" TALYS-error: Radiative capture model not possible",
     +    " for incident photons")')
        stop
      endif
      if (ldmodelracap.lt.1.or.ldmodelracap.gt.3) then
        write(*,'(" TALYS-error: 1 <= ldmodelracap <= 3")')
        stop
      endif
      do 270 Zix=0,numZ
        do 270 Nix=0,numN
          if (spectfacth(Zix,Nix).lt.0..or.
     +      spectfacth(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: 0. <= sfth <= 10.")')
            stop
          endif
          do 270 i=0,numlev
            if (spectfacexp(Zix,Nix,i).lt.0..or.
     +        spectfacexp(Zix,Nix,i).gt.10.) then
              write(*,'(" TALYS-error: 0. <= sfexp <= 10.")')
              stop
            endif
  270 continue
c
c 6. Check of values for pre-equilibrium
c
c epreeq     : on-set incident energy for preequilibrium calculation
c preeqmode  : designator for pre-equilibrium model
c mpreeqmode : designator for multiple pre-equilibrium model
c emulpre    : on-set incident energy for multiple preequilibrium
c phmodel    : particle-hole state density model
c pairmodel  : model for preequilibrium pairing energy
c pespinmodel: model for pre-equilibrium spin distribution or compound
c              spin distribution for pre-equilibrium cross section
c M2constant : constant for matrix element in exciton model
c M2limit    : constant for asymptotic value for matrix element
c M2shift    : constant for energy shift for matrix element
c Rpinu,.... : ratio for two-component matrix element
c Esurf0     : well depth for surface interaction
c Rgamma     : adjustable parameter for pre-equilibrium gamma decay
c msdbins    : number of energy points for DWBA calculation for MSD
c numenmsd   : maximum number of energy points for DWBA calculation for
c              MSD
c Emsdmin    : minimal outgoing energy for MSD calculation
c elwidth    : width of elastic peak in MeV
c flagpecomp : flag for Kalbach complex particle emission model
c Cstrip     : adjustable parameter for stripping/pick-up reactions
c Cknock     : adjustable parameter for knockout reactions
c Cbreak     : adjustable parameter for breakup reactions
c
      if ((epreeq.lt.0..or.epreeq.gt.Emaxtalys).and.epreeq.ne.-1.) then
        write(*,'(" TALYS-error: 0. <= epreeq < ",f10.3)') Emaxtalys
        stop
      endif
      if (preeqmode.lt.1.or.preeqmode.gt.4) then
        write(*,'(" TALYS-error: 1 <= preeqmode <= 4")')
        stop
      endif
      if (preeqmode.gt.3.and.k0.eq.0) then
        write(*,'(" TALYS-error: preeqmode <= 3 for incident photons")')
        stop
      endif
      if (mpreeqmode.lt.1.or.mpreeqmode.gt.2) then
        write(*,'(" TALYS-error: 1 <= mpreeqmode <= 2")')
        stop
      endif
      if (emulpre.lt.0..or.emulpre.gt.Emaxtalys) then
        write(*,'(" TALYS-error: 0. <= emulpre < ",f10.3)')  Emaxtalys
        stop
      endif
      if (phmodel.lt.1.or.phmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= phmodel <= 2")')
        stop
      endif
      if (pairmodel.lt.1.or.pairmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= pairmodel <= 2")')
        stop
      endif
      if (pespinmodel.lt.1.or.pespinmodel.gt.3) then
        write(*,'(" TALYS-error: 1 <= pespinmodel <= 3")')
        stop
      endif
      if (M2constant.lt.0..or.M2constant.gt.100.) then
        write(*,'(" TALYS-error: 0. <= M2constant <= 100.")')
        stop
      endif
      if (M2limit.lt.0..or.M2limit.gt.100.) then
        write(*,'(" TALYS-error: 0. <= M2limit <= 100.")')
        stop
      endif
      if (M2shift.lt.0..or.M2shift.gt.100.) then
        write(*,'(" TALYS-error: 0. <= M2shift <= 100.")')
        stop
      endif
      if (Rpipi.lt.0..or.Rpipi.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rpipi <= 100.")')
        stop
      endif
      if (Rnunu.lt.0..or.Rnunu.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rnunu <= 100.")')
        stop
      endif
      if (Rpinu.lt.0..or.Rpinu.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rpinu <= 100.")')
        stop
      endif
      if (Rnupi.lt.0..or.Rnupi.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rnupi <= 100.")')
        stop
      endif
      if (Esurf0.ne.-1..and.(Esurf0.lt.0..or.Esurf0.gt.38.)) then
        write(*,'(" TALYS-error: 0. <= Esurf <= 38.")')
        stop
      endif
      if (Rgamma.lt.0..or.Rgamma.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rgamma <= 100.")')
        stop
      endif
      if ((msdbins.lt.2.or.msdbins.gt.numenmsd/2-1).and.msdbins.ne.0)
     +  then
        write(*,'(" TALYS-error: 2 <= msdbins <=",i3)') numenmsd/2-1
        stop
      endif
      if (Emsdmin.lt.0.) then
        write(*,'(" TALYS-error: 0. <= E-in")')
        stop
      endif
      if (elwidth.lt.1.e-6.or.elwidth.gt.100.) then
        write(*,'(" TALYS-error: 1.e-6 <= elwidth <= 100.")')
        stop
      endif
      if (xscaptherm.ne.0..and.
     +  (xscaptherm.lt.1.e-20.or.xscaptherm.gt.1.e10)) then
        write(*,'(" TALYS-error: 1.e-20 <= xscaptherm <= 1.e10")')
        stop
      endif
      if (xsptherm.ne.0..and.
     +  (xsptherm.lt.1.e-20.or.xsptherm.gt.1.e10)) then
        write(*,'(" TALYS-error: 1.e-20 <= xsptherm <= 1.e10")')
        stop
      endif
      if (xsalphatherm.ne.0..and.
     +  (xsalphatherm.lt.1.e-20.or.xsalphatherm.gt.1.e10)) then
        write(*,'(" TALYS-error: 1.e-20 <= xsalphatherm <= 1.e10")')
        stop
      endif
      if (k0.eq.0.and.flagpecomp) then
        write(*,'(" TALYS-error: No pick-up and knock-out ",
     +    "mechanism for photonuclear reactions")')
        stop
      endif
      do 280 type=0,6
        if (Cstrip(type).lt.0..or.Cstrip(type).gt.100.) then
          write(*,'(" TALYS-error: 0. <= Cstrip <= 100.")')
          stop
        endif
        if (Cknock(type).lt.0..or.Cknock(type).gt.100.) then
          write(*,'(" TALYS-error: 0. <= Cknock <= 100.")')
          stop
        endif
        if (Cbreak(type).lt.0..or.Cbreak(type).gt.100.) then
          write(*,'(" TALYS-error: 0. <= Cbreak <= 100.")')
          stop
        endif
  280 continue
c
c 7. Check of values for level densities
c
c spincutmodel  : model for spin cutoff factor for ground state
c shellmodel    : model for shell correction energies
c kvibmodel     : model for vibrational enhancement
c ldmodel       : level density model
c alev          : level density parameter
c alimit        : asymptotic level density parameter
c gammald       : gamma-constant for asymptotic level density parameter
c ibar          : fission barrier
c numbar        : number of fission barriers
c deltaW        : shell correction in nuclear mass
c Nlow          : lowest discrete level for temperature matching
c Ntop          : highest discrete level for temperature matching
c E0            : constant of temperature formula
c beta2         : deformation parameter
c s2adjust      : adjustable constant (Z,A,barrier-dependent) for spin
c                 cutoff parameter
c Krotconstant  : normalization constant for rotational enhancement
c T             : nuclear temperature
c Exmatch       : matching point for Ex
c Pshift        : adjustable pairing shift
c pair          : pairing energy
c ctable        : constant to adjust tabulated level densities
c ptable        : constant to adjust tabulated level densities
c aadjust....   : adjustable factors for level density parameters
c                 (default 1.)
c cglobal       : global constant to adjust tabulated level densities
c pglobal       : global constant to adjust tabulated level densities
c g             : single-particle level density parameter
c gp            : single-particle proton level density parameter
c gn            : single-particle neutron level density parameter
c alphad        : alpha-constant for asymptotic level density parameter
c betald        : beta-constant for asymptotic level density parameter
c gammashell1   : gamma-constant for asymptotic level density parameter
c Pshiftconstant: global constant for pairing shift
c
      if (spincutmodel.lt.1.or.spincutmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= spincutmodel <= 2")')
        stop
      endif
      if (shellmodel.lt.1.or.shellmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= shellmodel <= 2")')
        stop
      endif
      if (kvibmodel.lt.1.or.kvibmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= kvibmodel <= 2")')
        stop
      endif
      do 310 Zix=0,numZ
        do 320 Nix=0,numN
          if (ldmodel(Zix,Nix).lt.1.or.ldmodel(Zix,Nix).gt.6) then
            write(*,'(" TALYS-error: 1 <= ldmodel <= 6")')
            stop
          endif
          if (ldmodel(Zix,Nix).ge.4) then
            inquire (file=path(1:lenpath)//
     +        'density/ground/hilaire/z050.tab',exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Microscopic HFB tables are not ",
     +        " installed: download the full TALYS package from ",
     +        "www.talys.eu")')
              stop
            endif
          endif
          if (alev(Zix,Nix).ne.0.) then
            if (alev(Zix,Nix).lt.1.or.alev(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 1. <= a <= 100.")')
              stop
            endif
          endif
          if (alimit(Zix,Nix).ne.0.) then
            if (alimit(Zix,Nix).lt.1.or.alimit(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 1. <= alimit <= 100.")')
              stop
            endif
          endif
          if (gammald(Zix,Nix).ne.-1.) then
            if (gammald(Zix,Nix).lt.0..or.gammald(Zix,Nix).gt.1.) then
              write(*,'(" TALYS-error: 0. <= gammald <= 1.")')
              stop
            endif
          endif
          do 330 ibar=0,numbar
            if (deltaW(Zix,Nix,ibar).ne.0.) then
              if (deltaW(Zix,Nix,ibar).lt.-20..or.
     +          deltaW(Zix,Nix,ibar).gt.20.) then
                write(*,'(" TALYS-error: -20. <= deltaW <= 20.")')
                stop
              endif
            endif
            if (Nlow(Zix,Nix,ibar).ne.-1) then
              if (Nlow(Zix,Nix,ibar).lt.0.or.
     +          Nlow(Zix,Nix,ibar).gt.200) then
                write(*,'(" TALYS-error: 0 <= Nlow <= 200")')
                stop
              endif
            endif
            if (Ntop(Zix,Nix,ibar).ne.-1) then
              if (Ntop(Zix,Nix,ibar).lt.0.or.
     +          Ntop(Zix,Nix,ibar).gt.200) then
                write(*,'(" TALYS-error: 0 <= Ntop <= 200")')
                stop
              endif
            endif
            if (Nlow(Zix,Nix,ibar).ne.-1.and.Ntop(Zix,Nix,ibar).ne.-1.
     +        and.Nlow(Zix,Nix,ibar).ge.Ntop(Zix,Nix,ibar)) then
              write(*,'(" TALYS-error: Ntop <= Nlow")')
              stop
            endif
            if (E0(Zix,Nix,ibar).ne.0.) then
              if (E0(Zix,Nix,ibar).lt.-10..or.
     +          E0(Zix,Nix,ibar).gt.10.) then
                write(*,'(" TALYS-error: -10. <= E0 <= 10.")')
                stop
              endif
            endif
            if (beta2(Zix,Nix,ibar).lt.-0.5.or.
     +        beta2(Zix,Nix,ibar).ge.1.5) then
              write(*,'(" TALYS-error: -0.5 <= beta2 < 1.5")')
              stop
            endif
            if (s2adjust(Zix,Nix,ibar).lt.0.01.or.
     +        s2adjust(Zix,Nix,ibar).gt.10.) then
              write(*,'(" TALYS-error: 0.01 <= s2adjust <= 10.")')
              stop
            endif
            if (Krotconstant(Zix,Nix,ibar).lt.0.01.or.
     +        Krotconstant(Zix,Nix,ibar).gt.100.) then
              write(*,'(" TALYS-error: 0.01 <= Krotconstant <= 100.")')
              stop
            endif
            if (T(Zix,Nix,ibar).ne.0.) then
              if (T(Zix,Nix,ibar).lt.1.e-3.or.
     +          T(Zix,Nix,ibar).gt.10.) then
                write(*,'(" TALYS-error: 1.e-3 <= T <= 10.")')
                stop
              endif
            endif
            if (Exmatch(Zix,Nix,ibar).ne.0.) then
              if (Exmatch(Zix,Nix,ibar).lt.0.1.or.
     +          Exmatch(Zix,Nix,ibar).gt.20.) then
                write(*,'(" TALYS-error: 0.1 <= Exmatch <= 20.")')
                stop
              endif
            endif
            if (Tadjust(Zix,Nix,ibar).lt.0.1.or.
     +        Tadjust(Zix,Nix,ibar).gt.10.) then
              write(*,'(" TALYS-error: 0.1 <= Tadjust <= 10.")')
              stop
            endif
            if (E0adjust(Zix,Nix,ibar).lt.0.1.or.
     +        E0adjust(Zix,Nix,ibar).gt.10.) then
              write(*,'(" TALYS-error: 0.1 <= E0adjust <= 10.")')
              stop
            endif
            if (Exmatchadjust(Zix,Nix,ibar).lt.0.2.or.
     +        Exmatchadjust(Zix,Nix,ibar).gt.2.) then
              write(*,'(" TALYS-error: 0.2 <= Exmatchadjust <= 2.")')
              stop
            endif
            if (Pshift(Zix,Nix,ibar).lt.-10..or.
     +        Pshift(Zix,Nix,ibar).gt.10.) then
              write(*,'(" TALYS-error: -10. <= Pshift <= 10.")')
              stop
            endif
            if (ctable(Zix,Nix,ibar).ne.0.) then
              if (ctable(Zix,Nix,ibar).lt.-10..or.ctable(Zix,Nix,ibar).
     +          gt.10.) then
                write(*,'(" TALYS-error: -10. <= ctable <= 10.")')
                stop
              endif
            endif
            if (ptable(Zix,Nix,ibar).ne.0.) then
              if (ptable(Zix,Nix,ibar).lt.-10..or.ptable(Zix,Nix,ibar).
     +          gt.10.) then
                write(*,'(" TALYS-error: -10. <= ptable <= 10.")')
                stop
              endif
            endif
            if (aadjust(Zix,Nix).lt.0.5.or.aadjust(Zix,Nix).gt.2.)
     +        then
              write(*,'(" TALYS-error: 0.5 <= aadjust <= 2.")')
              stop
            endif
            if (gnadjust(Zix,Nix).lt.0.5.or.gnadjust(Zix,Nix).gt.2.)
     +        then
              write(*,'(" TALYS-error: 0.5 <= gnadjust <= 2.")')
              stop
            endif
            if (gpadjust(Zix,Nix).lt.0.5.or.gpadjust(Zix,Nix).gt.2.)
     +        then
              write(*,'(" TALYS-error: 0.5 <= gpadjust <= 2.")')
              stop
            endif
  330     continue
          if (pair(Zix,Nix).lt.-10..or.pair(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: -10. <= pair <= 10.")')
            stop
          endif
          if (cglobal.ne.0.) then
            if (cglobal.lt.-10..or.cglobal.gt.10.) then
              write(*,'(" TALYS-error: -10. <= cglobal <= 10.")')
              stop
            endif
          endif
          if (pglobal.ne.0.) then
            if (pglobal.lt.-10..or.pglobal.gt.10.) then
              write(*,'(" TALYS-error: -10. <= pglobal <= 10.")')
              stop
            endif
          endif
          if (g(Zix,Nix).ne.0.) then
            if (g(Zix,Nix).lt.0.1.or.g(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 0.1 <= g <= 100.")')
              stop
            endif
          endif
          if (gp(Zix,Nix).ne.0.) then
            if (gp(Zix,Nix).lt.0.1.or.gp(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 0.1 <= gp <= 100.")')
              stop
            endif
          endif
          if (gn(Zix,Nix).ne.0.) then
            if (gn(Zix,Nix).lt.0.1.or.gn(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 0.1 <= gn <= 100.")')
              stop
            endif
          endif
          if (alphald(Zix,Nix).lt.0.01.or.alphald(Zix,Nix).gt.0.2) then
            write(*,'(" TALYS-error: 0.01 <= alphald <= 0.2")')
            stop
          endif
          if (betald(Zix,Nix).lt.-0.5.or.betald(Zix,Nix).gt.0.5) then
            write(*,'(" TALYS-error: -0.5 <= betald <= 0.5")')
            stop
          endif
          if (betald(Zix,Nix).lt.0..and.
     +      abs(betald(Zix,Nix)).gt.alphald(Zix,Nix)) then
            write(*,'(" TALYS-error: if betald<0, |betald|<alphald")')
            stop
          endif
          if (gammashell1(Zix,Nix).lt.0..or.gammashell1(Zix,Nix).gt.1.)
     +      then
            write(*,'(" TALYS-error: 0. <= gammashell1 <= 1.")')
            stop
          endif
          if (Pshiftconstant(Zix,Nix).lt.-5..or.
     +      Pshiftconstant(Zix,Nix).gt.5.) then
            write(*,'(" TALYS-error: -5. <= Pshiftconstant <= 5.")')
            stop
          endif
  320   continue
  310 continue
c
c There are many input possibilities for the energy dependent level
c density parameter of the Ignatyuk formula. The required parameters
c are alev, alimit, gammald and deltaW. The Ignatyuk formula implies
c that they can not all be given at the same time in the input file.
c
c gammashell2   : gamma-constant for asymptotic level density parameter
c pairconstant  : constant for pairing energy systematics
c Ufermi        : energy of Fermi distribution for damping of
c                 ground-state rotational effects
c cfermi        : width of Fermi distribution for damping of
c               : ground-state rotational effects
c Ufermibf      : energy of Fermi distribution for damping of barrier
c               : rotational effects
c cfermibf      : width of Fermi distribution for damping of barrier
c               : rotational effects
c Kph           : constant for single-particle level density parameter
c                 (g=A/Kph)
c Rspincut      : adjustable constant (global) for spin cutoff factor
c
      do 340 Zix=0,numZ
        do 350 Nix=0,numN
          do 360 ibar=0,numbar
            if (alev(Zix,Nix).ne.0.and.deltaW(Zix,Nix,ibar).ne.0..and.
     +        alimit(Zix,Nix).ne.0.and.gammald(Zix,Nix).ne.-1.) then
              write(*,'(" TALYS-error: Level density conflict - a,",
     +          " deltaW, alimit and gammald are ALL given",
     +          " in the input for Z=",i3," A=",i3,
     +          " fission barrier=",i3)') Zinit-Zix,Ainit-Zix-Nix,ibar
              stop
            endif
  360     continue
  350   continue
  340 continue
      if (gammashell2.lt.0..or.gammashell2.gt.0.2) then
        write(*,'(" TALYS-error: 0. <= gammashell2 <= 0.2")')
        stop
      endif
      if (pairconstant.lt.0..or.pairconstant.gt.30.) then
        write(*,'(" TALYS-error: 0. <= pairconstant <= 30.")')
        stop
      endif
      if (Ufermi.lt.0..or.Ufermi.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= Ufermi <= 1000.")')
        stop
      endif
      if (cfermi.le.0..or.cfermi.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= cfermi <= 1000.")')
        stop
      endif
      if (Ufermibf.lt.0..or.Ufermibf.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= Ufermibf <= 1000.")')
        stop
      endif
      if (cfermibf.le.0..or.cfermibf.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= cfermibf <= 1000.")')
        stop
      endif
      if (Kph.lt.1..or.Kph.gt.100.) then
        write(*,'(" TALYS-error: 1. <= Kph <= 100.")')
        stop
      endif
      if (Rspincut.le.0..or.Rspincut.gt.10.) then
        write(*,'(" TALYS-error: 0. < Rspincut <= 10.")')
        stop
      endif
c
c 8. Check of values for fission
c
c flagfission: flag for fission
c flagfisout : flag for output of fission information
c fismodel   : fission model
c fismodelalt: alternative fission model for default barriers
c gefran     : number of random events for GEF calculation
c axtype     : type of axiality of barrier
c                 1: axial symmetry
c                 2: left-right asymmetry
c                 3: triaxial and left-right symmetry
c                 4: triaxial no left-right symmetry
c                 5: no symmetry
c fbarrier   : height of fission barrier
c fwidth     : width of fission barrier
c fbaradjust.: adjustable factors for fission parameters
c              (default 1.)
c Rtransmom  : normalization constant for moment of inertia for
c              transition states
c Rclass2mom : normalization constant for moment of inertia for
c              class 2 states
c widthc2    : width of class2 states
c betafiscor : adjustable factor for fission path width
c vfiscor    : adjustable factor for fission path height
c
      if ((flagfission.or.flagfisout).and.Atarget.le.120) then
        write(*,'(" TALYS-error: Fission not allowed for A <= 120")')
        stop
      endif
      if (fismodel.lt.1.or.fismodel.gt.5) then
        write(*,'(" TALYS-error: 1 <= fismodel <= 5")')
        stop
      endif
      if (fismodelalt.lt.3.or.fismodelalt.gt.4) then
        write(*,'(" TALYS-error: 3 <= fismodelalt <= 4")')
        stop
      endif
      if (gefran.lt.1000.or.gefran.gt.1000000) then
        write(*,'(" TALYS-error: 1000 <= gefran <= 1000000")')
        stop
      endif
      do 410 Zix=0,numZ
        do 420 Nix=0,numN
          do 430 ibar=1,numbar
            fax0=axtype(Zix,Nix,ibar)
            if (fax0.ne.0.and.(fax0.lt.1.or.fax0.gt.5)) then
              write(*,'(" TALYS-error: 1 <= type of axiality <= 5")')
              stop
            endif
            fbar0=fbarrier(Zix,Nix,ibar)
            if (fbar0.ne.0..and.(fbar0.lt.0..or.fbar0.gt.100.)) then
              write(*,'(" TALYS-error: 0.<= fission barrrier <=100.")')
              stop
            endif
            fbar0=fbaradjust(Zix,Nix,ibar)
            if (fbar0.ne.0..and.(fbar0.lt.0.1.or.fbar0.gt.10.)) then
              write(*,'(" TALYS-error: 0.1 <= fbaradjust <= 10.")')
              stop
            endif
            fhw0=fwidth(Zix,Nix,ibar)
            if (fhw0.ne.0..and.(fhw0.lt.0.01.or.fhw0.gt.10.)) then
              write(*,'(" TALYS-error: 0.01 <= fission width <= 10.")')
              stop
            endif
            fhw0=fwidthadjust(Zix,Nix,ibar)
            if (fhw0.ne.0..and.(fhw0.lt.0.1.or.fhw0.gt.10.)) then
              write(*,'(" TALYS-error: 0.1 <= fwidthadjust <= 10.")')
              stop
            endif
            fR0=Rtransmom(Zix,Nix,ibar)
            if (fR0.ne.0..and.(fR0.lt.0.1.or.fR0.gt.10.)) then
              write(*,'(" TALYS-error: 0.1 <= Rtransmom <= 10.")')
              stop
            endif
            fR0=Rclass2mom(Zix,Nix,ibar)
            if (fR0.ne.0..and.(fR0.lt.0.1.or.fR0.gt.10.)) then
              write(*,'(" TALYS-error: 0.1 <= Rclass2mom <= 10.")')
              stop
            endif
  430     continue
          fhw0=betafiscor(Zix,Nix)
          if (fhw0.lt.0.1.or.fhw0.gt.10.) then
            write(*,'(" TALYS-error: 0.1 <= betafiscor <=10.")')
            stop
          endif
          fbar0=vfiscor(Zix,Nix)
          if (fbar0.lt.0.1.or.fbar0.gt.10.) then
            write(*,'(" TALYS-error: 0.1 <= vfiscor <=10.")')
            stop
          endif
  420   continue
  410 continue
c
c 9. Check of values for output
c
c eadd     : on-set incident energy for addition of discrete states
c            to spectra
c eaddel   : on-set incident energy for addition of elastic peak
c            to spectra
c ddxmode  : mode for double-differential cross sections: 0: None,
c            1: Angular distributions, 2: Spectra per angle, 3: Both
c ddxecount: counter for double-differential cross section files
c fileddxe : designator for double-differential cross sections on
c            separate file: angular distribution
c enincmax : maximum incident energy
c ddxacount: counter for double-differential cross section files
c fileddxa : designator for double-differential cross sections on
c            separate file: spectrum per angle
c
      if (eadd.lt.0..or.eadd.gt.Emaxtalys) then
        write(*,'(" TALYS-error: 0. <= eadd < ",f10.3)')  Emaxtalys
        stop
      endif
      if (eaddel.lt.0..or.eaddel.gt.Emaxtalys) then
        write(*,'(" TALYS-error: 0. <= eaddel < ",f10.3)')  Emaxtalys
        stop
      endif
      if (ddxmode.lt.0.or.ddxmode.gt.3) then
        write(*,'(" TALYS-error: 0 <= ddxmode <= 3")')
        stop
      endif
      do 510 type=0,6
        do 520 i=1,ddxecount(type)
          value=fileddxe(type,i)
          if (value.lt.0.or.value.gt.enincmax) then
            write(*,'(" TALYS-error: 0. <= fileddxe <=",f7.3)') enincmax
            stop
          endif
  520   continue
        do 530 i=1,ddxacount(type)
          value=fileddxa(type,i)
          if (value.lt.0.or.value.gt.180.) then
            write(*,'(" TALYS-error: 0. <= fileddxa <= 180.")')
            stop
          endif
  530   continue
  510 continue
c
c 10. Check of values energy-dependent parameter adjustment
c
c Nadjust  : number of adjustable parameters
c adjustpar: local adjustment parameters
c Ea       : start energy of local adjustment
c Eb       : end energy of local adjustment
c Em       : intermediate energy of local adjustment
c D        : depth of local adjustment
c adjustkey: keyword for local adjustment
c
      do 610 n=1,Nadjust
        if (adjustfile(i)(1:1).ne.' ') goto 610
        Ea=adjustpar(n,1)
        Eb=adjustpar(n,2)
        Em=adjustpar(n,3)
        D=adjustpar(n,4)
        if ((Ea.ge.Eb).or.(Ea.ge.Em).or.(Em.ge.Eb)) then
          write(*,'(" TALYS-error: energy range for adjustment should",
     +      " be given as follows: Ea Eb Em D, with Ea < Em < Eb",
     +      " for keyword ",a80)') adjustkey(n)
          stop
        endif
        if (D.le.0..or.D.ge.10.) then
          write(*,'(" TALYS-error: 0. < D <= 10. for keyword ",a80)')
     +      adjustkey(n)
          stop
        endif
        do 620 m=1,Nadjust
          if (m.eq.n) goto 620
          if (adjustkey(m).ne.adjustkey(n)) goto 620
          Ea2=adjustpar(m,1)
          Eb2=adjustpar(m,2)
          if ((Ea2.gt.Ea.and.Ea2.lt.Eb).or.(Eb2.gt.Ea.and.Eb2.lt.Eb))
     +      then
            write(*,'(" TALYS-error: overlapping energy ranges for ",
     +        "keyword ",a80)') adjustkey(n)
            stop
          endif
  620   continue
  610 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
