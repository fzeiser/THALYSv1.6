      subroutine input1
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Read input for first set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      projexist,massexist,elemexist,enerexist,lexist,fexist
      character*1  ch
      character*15 bestchar
      character*40 bestpath
      character*80 word(40),key,value,bestfile
      integer      i,i2,inull,type,iz,nbest,J,k,pbeg,parity,inum,negrid,
     +             lenbest
      real         Ein,etmp,enincF,deninc,E
c
c ************ Read first set of variables from input lines ************
c
c energyfile : file with incident energies
c projexist  : logical for existence of projectile
c massexist  : logical for existence of mass
c elemexist  : logical for existence of element
c enerexist  : logical for existence of energy
c flagnatural: flag for calculation of natural element
c flagmicro  : flag for completely microscopic Talys calculation
c flagastro  : flag for calculation of astrophysics reaction rate
c flagbest   : flag to use best set of adjusted parameters
c bestpath   : alternative directory for best values
c eninc      : incident energy in MeV
c numJ       : maximal J-value
c Exdist     : excitation energy of population distribution
c Pdistex    : population distribution, spin-independent
c Pdist      : population distribution per spin and parity
c Ztarget    : charge number of target nucleus
c Ltarget    : excited level of target
c Starget    : symbol of target nucleus
c enincF     : final incident energy
c deninc     : incident energy increment
c
c 1. Initializations
c
      energyfile='                                                     '
      projexist=.false.
      massexist=.false.
      elemexist=.false.
      enerexist=.false.
      flagnatural=.false.
      flagmicro=.false.
      flagastro=.false.
      flagbest=.false.
      bestpath='                                        '
      do 5 i=0,numen6+2
        eninc(i)=0.
    5 continue
      do 7 i=0,numex
        Exdist(i)=0.
        Pdistex(i)=0.
        do 8 parity=-1,1,2
          do 9 J=0,numJ
            Pdist(i,J,parity)=0.
    9     continue
    8   continue
    7 continue
      Ztarget=0
      Ltarget=0
      Starget='  '
      enincF=0.
      deninc=0.
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
c 2. The projectile is read
c
c ptype0: type of incident particle
c
        if (key.eq.'projectile') then
          projexist=.true.
          ptype0=ch
          goto 10
        endif
c
c 3. The target mass is read
c
c Atarget: mass number of target nucleus
c
        if (key.eq.'mass') then
          massexist=.true.
          read(value,*,end=500,err=500) Atarget
          goto 10
        endif
c
c 4. The nuclear symbol or charge number is read
c
c numelem: number of elements
c
        if (key.eq.'element') then
          elemexist=.true.
          if (ch.ge.'0'.and.ch.le.'9') then
            read(value,*,end=500,err=500) Ztarget
            if (Ztarget.lt.1.or.Ztarget.gt.numelem) goto 500
            goto 10
          else
            read(value,'(a2)',end=500,err=500) Starget
            Starget(1:1)=char(ichar(Starget(1:1))-32)
            goto 10
          endif
        endif
c
c 5. The level of the target is read
c
        if (key.eq.'ltarget') then
          read(value,*,end=500,err=500) Ltarget
          goto 10
        endif
c
c 6. The incident energy or file with incident energies is read
c
        if (key.eq.'energy') then
          enerexist=.true.
          if ((ch.ge.'0'.and.ch.le.'9').or.ch.eq.'.') then
            read(value,*,end=500,err=500) eninc(1)
            read(word(3),*,end=10,err=10) enincF
            read(word(4),*,end=10,err=10) deninc
            goto 10
          else
            eninc(1)=0.
            energyfile=value
            goto 10
          endif
        endif
c
c 7. Test for completely microscopic and/or astrophysical calculation
c
        if (key.eq.'micro') then
          if (ch.eq.'n') flagmicro=.false.
          if (ch.eq.'y') flagmicro=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
        if (key.eq.'astro') then
          if (ch.eq.'n') flagastro=.false.
          if (ch.eq.'y') flagastro=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
c
c 8. Possibility to use "best" adjusted parameter sets
c
        if (key.eq.'best') then
          if (ch.eq.'n') flagbest=.false.
          if (ch.eq.'y') flagbest=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
        if (key.eq.'bestpath') then
          bestpath=value
          goto 10
        endif
   10 continue
c
c The four main keywords MUST be present in the input file.
c
      if (.not.projexist) then
        write(*,'(" TALYS-error: projectile must be given")')
        stop
      endif
      if (.not.massexist) then
        write(*,'(" TALYS-error: mass must be given")')
        stop
      endif
      if (.not.elemexist) then
        write(*,'(" TALYS-error: element must be given")')
        stop
      endif
      if (.not.enerexist) then
        write(*,'(" TALYS-error: energy must be given")')
        stop
      endif
c
c Manual input of structure path and null device.
c
c path   : directory containing structure files to be read
c lenpath: length of pathname
c nulldev: null device
c
      do 50 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
        if (key.eq.'strucpath') then
          lenpath=0
          do 60 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              lenpath=lenpath+1
              path(lenpath:lenpath)=ch
            endif
   60     continue
        endif
        inull=0
        if (key.eq.'nulldev') then
          nulldev='             '
          do 70 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              inull=inull+1
              nulldev(inull:inull)=ch
            endif
   70     continue
        endif
   50 continue
c
c Test to check accessibility of structure files and null device
c
      if (path(lenpath:lenpath).ne.'/') then
        lenpath=lenpath+1
        path(lenpath:lenpath)='/'
      endif
      if (lenpath.gt.60) then
        write(*,'(" TALYS-warning: path name should contain 60",
     +    " characters or less")')
      endif
      inquire (file=path(1:lenpath)//'abundance/z001',exist=lexist)
      if (.not.lexist) then
        write(*,'(" TALYS-error: Structure database not installed:",
     +    " change path in machine.f or strucpath keyword",
     +    " in input file")')
        stop
      endif
      if (inull.gt.13) then
        write(*,'(" TALYS-warning: null device should contain 13",
     +    " characters or less")')
      endif
c
c ************* Process first set of input variables *******************
c
c 1. Assignment of index k0 to incident particle
c
c flaginitpop: flag for initial population distribution
c parsym     : symbol of particle
c k0         : index of incident particle
c
c Throughout TALYS, the initial particle can always be identified
c by the index k0.
c
      flaginitpop=.false.
      k0=0
      do 110 type=0,6
        if (ptype0.eq.parsym(type)) then
          k0=type
          goto 200
        endif
  110 continue
c
c It is also possible to define a population distribution as the
c initial state, through the keyword projectile 0. In that case,
c we assign a photon projectile.
c
      if (ptype0.eq.'0') flaginitpop=.true.
c
c 2. Identification of target and initial compound nucleus
c
c nuc: symbol of nucleus
c
  200 if (Ztarget.eq.0) then
        do 210 iz=1,numelem
          if (nuc(iz).eq.Starget) then
            Ztarget=iz
            goto 220
          endif
  210   continue
      else
        Starget=nuc(Ztarget)
      endif
c
c A calculation for a natural element is specified by target mass 0
c
c abundance: subroutine for natural abundances
c iso      : counter for isotope
c nbest    : number of lines in best file
c isotope  : isotope number of residual nucleus
c Ntarget  : neutron number of target nucleus
c Zinit    : charge number of initial compound nucleus
c parZ     : charge number of particle
c Ninit    : neutron number of initial compound nucleus
c parN     : neutron number of particle
c Ainit    : mass number of initial compound nucleus
c
  220 if (Atarget.eq.0) then
        flagnatural=.true.
        if (iso.eq.1) then
          call abundance
        else
          if (flagbest) then
            nbest=nlines-nlines0
            do 230 i=nbest+1,nlines
              inline(i-nbest)=inline(i)
  230       continue
            nlines=nlines0
          endif
        endif
        Atarget=isotope(iso)
      endif
      Ntarget=Atarget-Ztarget
      Zinit=Ztarget+parZ(k0)
      Ninit=Ntarget+parN(k0)
      Ainit=Zinit+Ninit
c
c 3. Determine incident energies
c
c Ein: incident energy
c nin: counter for incident energies
c
c 1. Normal case of a projectile with a target nucleus
c
      Ein=0.
      fexist=.false.
      if (.not.flaginitpop) then
c
c A. If no incident energy is given in the input file, incident energies
c    should be read from a file.
c
c incidentgrid: subroutine with predefined incident energy grids
c
        eninc(0)=0.
        if (eninc(1).eq.0.) then
          inquire (file=energyfile,exist=lexist)
          if (.not.lexist) then
            call incidentgrid(energyfile(1:14),fexist)
            if (fexist) then
              goto 300
            else
              write(*,'(" TALYS-error: give a single incident energy ",
     +          "in the input file using the energy keyword ")')
              write(*,'(14x,"or specify a range of incident energies ",
     +          "in a file ")')
              write(*,'(14x,"or give a correct name for a pre-defined",
     +          " energy grid ",a73)') energyfile
              stop
            endif
          endif
          nin=0
          open (unit=2,status='old',file=energyfile)
  310     read(2,*,end=320,err=510) Ein
          if (Ein.ne.0.) then
            nin=nin+1
c
c There is a maximum number of incident energies
c
c numenin : maximal number of incident energies
c
            if (nin.gt.numenin) then
              write(*,'(" TALYS-error: there are more than",i4,
     +          " incident energies in file ",a73)') numenin,energyfile
              write(*,'(" numenin in talys.cmb should be increased")')
              stop
            endif
            eninc(nin)=Ein
          endif
          goto 310
  320     close (unit=2)
          if (nin.eq.0) then
            write(*,'(" TALYS-error: there are no",
     +        " incident energies in file ",a73)') energyfile
            stop
          endif
c
c Sort incident energies in ascending order and remove double points
c
c numinc: number of incident energies
c
          do 322 i=1,nin
            do 324 k=1,i
              if (eninc(i).ge.eninc(k)) goto 324
              etmp=eninc(i)
              eninc(i)=eninc(k)
              eninc(k)=etmp
  324       continue
  322     continue
          numinc=nin
          do 326 i=1,nin-1
            if (eninc(i).eq.eninc(i+1)) then
              do 328 k=i+1,nin
                eninc(k)=eninc(k+1)
  328         continue
              numinc=numinc-1
            endif
  326     continue
c
c The minimum and maximum value of all the incident energies is
c determined.
c
c enincmin: minimum incident energy
c enincmax: maximum incident energy
c
          enincmin=eninc(1)
          enincmax=eninc(1)
          do 330 nin=2,numinc
            enincmin=min(enincmin,eninc(nin))
            enincmax=max(enincmax,eninc(nin))
  330     continue
        else
          if (enincF.eq.0.) then
c
c B1. Single value given in the user input file
c
            numinc=1
            enincmin=eninc(1)
            enincmax=eninc(1)
          else
c
c B2. Energy grid based on input values E1, E2, dE
c
            if (enincF.le.eninc(1)) then
              write(*,'(" TALYS-error: final incident energy should",
     +          " be larger than the first ",a80)') inline(i)
              stop
            endif
            if (deninc.lt.0.) then
              write(*,'(" TALYS-error: energy increment should",
     +          " be larger than zero ",a80)') inline(i)
              stop
            endif
            if (deninc.eq.0.) then
              negrid=10
              deninc=(enincF-eninc(1))/(negrid-1)
            endif
            nin=1
  335       nin=nin+1
            if (nin.gt.numenin) then
              write(*,'(" TALYS-error: number of incident energies ",
     +          " greater than ",i4)') numenin
              stop
            endif
            E=eninc(nin-1)+deninc
            if (E.lt.enincF-1.e-4) then
              eninc(nin)=E
              goto 335
            else
              eninc(nin)=enincF
              numinc=nin
              enincmin=eninc(1)
              enincmax=eninc(nin)
            endif
          endif
        endif
      else
c
c 2. Population distribution as the initial state
c
c npopbins: number of excitation energy bins for population distribution
c npopJ   : number of spins for population distribution
c npopP   : number of parities for population distribution
c pbeg    : help variable
c numbins : maximal number of continuum excitation energy bins
c
        inquire (file=energyfile,exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: if projectile 0, specify a range",
     +        " of excitation energies in a file ",a73)') energyfile
          stop
        endif
        open (unit=2,status='old',file=energyfile)
        read(2,*,end=510,err=510) npopbins,npopJ,npopP,eninc(1)
        if (npopbins.lt.2.or.npopbins.gt.numbins) then
          write(*,'(" TALYS-error: 2 <= bins <=",i3," in population ",
     +      "distribution file")') numbins
          stop
        endif
        if (npopJ.lt.0.or.npopJ.gt.numJ+1) then
          write(*,'(" TALYS-error: 0 <= number of spins <=",i3,
     +      " + 1 in population distribution file")') numJ
          stop
        endif
        if (npopJ.gt.0.and.(npopP.lt.1.or.npopP.gt.2)) then
          write(*,'(" TALYS-error: 1 <= number of parities <= 2",
     +      " in population distribution file")')
          stop
        endif
c
c Only excitation energy distribution (no spins)
c
        if (npopJ.eq.0) then
          do 350 nin=1,npopbins
            read(2,*,end=510,err=510) Exdist(nin),Pdistex(nin)
  350     continue
        else
c
c Spin-dependent excitation energy distribution (no total)
c
          if (npopP.eq.1) then
            pbeg=1
          else
            pbeg=-1
          endif
          do 360 parity=pbeg,1,2
            do 370 nin=1,npopbins
                read(2,*,end=510,err=510) Exdist(nin),
     +            (Pdist(nin,J,parity),J=0,npopJ-1)
  370       continue
  360     continue
          do 375 nin=2,npopbins
            if (Exdist(nin).le.Exdist(nin-1)) then
              write(*,'(" TALYS-error: excitation energies must",
     +          " be given in ascending order, or the number",
     +          " of population bins is not correct")')
              stop
            endif
  375     continue
          if (npopP.eq.1) then
            do 380 nin=1,npopbins
              do 380 J=0,npopJ-1
                Pdist(nin,J,-1)=Pdist(nin,J,1)
  380       continue
          endif
        endif
        close (unit=2)
        numinc=1
        enincmin=eninc(1)
        enincmax=eninc(1)
      endif
c
c In case of built-in energy range, write an explicit 'energies' file
c
  300 if (enincF.gt.0..or.fexist) then
        open (unit=2,status='unknown',file='energies')
        do 385 nin=1,numinc
          write(2,'(1p,g12.4)') eninc(nin)
  385   continue
        close (unit=2)
      endif
c
c If requested by input: retrieve best set of adjusted input parameters
c
c bestchar: help variable
c bestfile: adjusted "best" parameter file
c convert : subroutine to convert input line from upper case to lowercase
c
      if (flagbest) then
        bestchar='z000a000x.talys'
        write(bestchar(2:4),'(i3.3)') Ztarget
        write(bestchar(6:8),'(i3.3)') Atarget
        write(bestchar(9:9),'(a1)') ptype0
        if (bestpath(1:1).eq.' ')
     +    bestpath='best/                                   '
        do 390 i=1,40
          if (bestpath(i:i).eq.' ') then
            if (bestpath(i-1:i-1).ne.'/') then
              bestpath(i:i)='/'
              lenbest=i
            else
              lenbest=i-1
            endif
            goto 400
          endif
  390   continue
  400   if (Starget(2:2).eq.' ') then
          if (Ltarget.eq.0) then
            write(bestpath(lenbest+1:lenbest+5),'(a1,i3.3,"/")')
     +        Starget(1:1),Atarget
            write(bestpath(lenbest+6:lenbest+20),'(a15)') bestchar
          else
            write(bestpath(lenbest+1:lenbest+6),'(a1,i3.3,"m/")')
     +        Starget(1:1),Atarget
            write(bestpath(lenbest+7:lenbest+21),'(a15)') bestchar
          endif
        else
          if (Ltarget.eq.0) then
            write(bestpath(lenbest+1:lenbest+6),'(a2,i3.3,"/")')
     +        Starget(1:2),Atarget
            write(bestpath(lenbest+7:lenbest+21),'(a15)') bestchar
          else
            write(bestpath(lenbest+1:lenbest+7),'(a2,i3.3,"m/")')
     +        Starget(1:2),Atarget
            write(bestpath(lenbest+8:lenbest+22),'(a15)') bestchar
          endif
        endif
        bestfile=path(1:lenpath)//bestpath
        inquire (file=bestfile,exist=lexist)
        if (.not.lexist) return
        open (unit=3,status='old',file=bestfile)
        inum=0
  410   read(3,'(a80)',end=420) key
        inum=inum+1
        i=numlines-inum
        inline(i)=key
        call convert(i)
        goto 410
  420   close (unit=3)
        if (inum.gt.0) then
          do 430 i=nlines,1,-1
            inline(i+inum)=inline(i)
  430     continue
          do 440 i=1,inum
            inline(i)=inline(numlines-i)
  440     continue
          nlines=nlines+inum
        endif
      endif
      return
  500 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
  510 write(*,'(" TALYS-error: Problem in file ",a73)') energyfile
      write(*,'(" after E= ",e12.5)') Ein
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
