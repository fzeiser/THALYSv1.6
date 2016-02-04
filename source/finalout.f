        subroutine finalout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 18, 2013
c | Task  : Output of final results
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*6   discfile,contfile
      character*9   totfile,reactionstring(0:6)
      character*10  binfile
      character*11  fisfile
      character*12  rpfile,isofile,xsfile
      character*19  gamfile
      character*80  string
      character*126 stringtot
      integer       istat,nen,type,l,Acomp,Zcomp,Ncomp,Z,A,nex,Zix,Nix,
     +              npart,ia,ih,it,id,ip,in,ident,idc,i1,i2
      real          Egam
c
c Write model parameters to separate file
c
c flagpartable  : flag for output of model parameters on separate file
c ldmodel       : level density model
c alphad        : alpha-constant for asymptotic level density parameter
c betald        : beta-constant for asymptotic level density parameter
c gammashell1   : gamma-constant for asymptotic level density parameter
c gammashell2   : gamma-constant for asymptotic level density parameter
c pairconstant  : constant for pairing energy systematics
c Pshiftconstant: global constant for pairing shift
c Rspincut      : adjustable constant (global) for spin cutoff factor
c cglobal,...   : global constant to adjust tabulated level densities
c flagcolall    : flag for collective enhancement of level density
c                 for all nuclides
c Ufermi        : energy of Fermi distribution for damping of
c               : ground-state rotational effects
c cfermi        : width of Fermi distribution for damping of
c               : ground-state rotational effects
c Ufermibf      : energy of Fermi distribution for damping of barrier
c               : rotational effects
c cfermibf      : width of Fermi distribution for damping of barrier
c               : rotational effects
c phmodel       : particle-hole state density model
c Kph           : constant for single-particle level density parameter
c                 (g=A/Kph)
c gnorm         : gamma normalization factor
c xscaptherm    : thermal capture cross section
c M2constant    : overall constant for matrix element in exciton model
c M2limit       : constant for asymptotic value for matrix element
c M2shift       : constant for energy shift for matrix element
c Rpinu,Rnupi...: ratio for two-component matrix element
c Rgamma        : adjustable parameter for pre-equilibrium gamma decay
c Esurf         : well depth for surface interaction
c parsym        : symbol of particle
c xspreeqtotps  : preequilibrium cross section per particle type for
c                 pickup and stripping
c xspreeqtotki  : preequilibrium cross section per particle type for
c                 knockout and inelastic
c xspreeqtotbu  : preequilibrium cross section per particle type for
c                 breakup
c Cstrip        : adjustable parameter for stripping/pick-up reactions
c Cknock        : adjustable parameter for knockout reactions
c Cbreak        : adjustable parameter for breakup reactions
c flagjlm       : flag for using semi-microscopic JLM OMP
c alphaomp      : alpha optical model (1=normal, 2= McFadden-Satchler,
c                 3-5= folding potential)
c aradialcor    : adjustable parameter for shape of DF alpha potential
c adepthcor     : adjustable parameter for depth of DF alpha potential
c v1adjust..    : adjustable factors for OMP (default 1.)
c tladjust      : adjustable factor for Tlj (default 1.)
c k0            : index of incident particle
c parinclude    : logical to include outgoing particle
c flagcomp      : flag for compound nucleus calculation
c Rprime        : potential scattering radius
c Sstrength     : s,p,d,etc-wave strength function
c D0theo        : theoretical s-wave resonance spacing
c
      if (flagpartable) then
        write(11,'("##")')
        write(11,'("## General parameters")')
        write(11,'("##")')
        write(11,'("## Level density")')
        write(11,'("##")')
        if (ldmodel(0,0).le.3) then
          write(11,'("alphald        ",f10.5)') alphald(0,0)
          write(11,'("betald         ",f10.5)') betald(0,0)
          write(11,'("gammashell1    ",f10.5)') gammashell1(0,0)
          write(11,'("gammashell2    ",f10.5)') gammashell2
          write(11,'("pairconstant   ",f10.5)') pairconstant
          write(11,'("pshiftconstant ",f10.5)') Pshiftconstant(0,0)
          write(11,'("Rspincut       ",f10.5)') Rspincut
        endif
        write(11,'("cglobal        ",1p,e12.5)') cglobal
        write(11,'("pglobal        ",1p,e12.5)') pglobal
        if (flagcolall) then
          write(11,'("Ufermi         ",f10.5)') Ufermi
          write(11,'("cfermi         ",f10.5)') cfermi
          write(11,'("Ufermibf       ",f10.5)') Ufermibf
          write(11,'("cfermibf       ",f10.5)') cfermibf
        endif
        if (phmodel.eq.1) write(11,'("Kph            ",f10.5)') Kph
        write(11,'("##")')
        write(11,'("## Gamma-ray")')
        write(11,'("##")')
        write(11,'("gnorm          ",f10.5)') gnorm
        write(11,'("xscaptherm     ",1p,e12.5)') xscaptherm
        write(11,'("##")')
        write(11,'("## Pre-equilibrium")')
        write(11,'("##")')
        write(11,'("M2constant     ",f10.5)') M2constant
        write(11,'("M2limit        ",f10.5)') M2limit
        write(11,'("M2shift        ",f10.5)') M2shift
        write(11,'("Rpipi          ",f10.5)') Rpipi
        write(11,'("Rnunu          ",f10.5)') Rnunu
        write(11,'("Rpinu          ",f10.5)') Rpinu
        write(11,'("Rnupi          ",f10.5)') Rnupi
        write(11,'("Rgamma         ",f10.5)') Rgamma
        write(11,'("Esurf          ",f10.5)') Esurf
        do type=1,6
          if (xspreeqtotps(type).gt.0.)
     +      write(11,'("Cstrip        ",a1,f10.5)') parsym(type),
     +      Cstrip(type)
          if (xspreeqtotki(type).gt.0.)
     +      write(11,'("Cknock        ",a1,f10.5)') parsym(type),
     +      Cknock(type)
          if (xspreeqtotbu(type).gt.0.)
     +      write(11,'("Cbreak        ",a1,f10.5)') parsym(type),
     +      Cbreak(type)
        enddo
        write(11,'("##")')
        write(11,'("## Optical model")')
        write(11,'("##")')
        if (flagjlm) then
          write(11,'("lvadjust       ",f10.5)') lvadjust
          write(11,'("lwadjust       ",f10.5)') lwadjust
          write(11,'("lv1adjust      ",f10.5)') lv1adjust
          write(11,'("lw1adjust      ",f10.5)') lw1adjust
          write(11,'("lvsoadjust     ",f10.5)') lvsoadjust
          write(11,'("lwsoadjust     ",f10.5)') lwsoadjust
        endif
        if (alphaomp.ge.3) then
          write(11,'("aradialcor     ",f10.5)') aradialcor
          write(11,'("adepthcor      ",f10.5)') adepthcor
        endif
        do type=1,6
          if (.not.flagjlm.or.type.gt.2) then
            write(11,'("v1adjust      ",a1,f10.5)') parsym(type),
     +        v1adjust(type)
            write(11,'("v2adjust      ",a1,f10.5)') parsym(type),
     +        v2adjust(type)
            write(11,'("v3adjust      ",a1,f10.5)') parsym(type),
     +        v3adjust(type)
            write(11,'("v4adjust      ",a1,f10.5)') parsym(type),
     +        v4adjust(type)
            write(11,'("rvadjust      ",a1,f10.5)') parsym(type),
     +        rvadjust(type)
            write(11,'("avadjust      ",a1,f10.5)') parsym(type),
     +        avadjust(type)
            write(11,'("w1adjust      ",a1,f10.5)') parsym(type),
     +        w1adjust(type)
            write(11,'("w2adjust      ",a1,f10.5)') parsym(type),
     +        w2adjust(type)
            write(11,'("w3adjust      ",a1,f10.5)') parsym(type),
     +        w3adjust(type)
            write(11,'("w4adjust      ",a1,f10.5)') parsym(type),
     +        w4adjust(type)
            write(11,'("d1adjust      ",a1,f10.5)') parsym(type),
     +        d1adjust(type)
            write(11,'("d2adjust      ",a1,f10.5)') parsym(type),
     +        d2adjust(type)
            write(11,'("d3adjust      ",a1,f10.5)') parsym(type),
     +        d3adjust(type)
            write(11,'("rvdadjust     ",a1,f10.5)') parsym(type),
     +        rvdadjust(type)
            write(11,'("avdadjust     ",a1,f10.5)') parsym(type),
     +        avdadjust(type)
            write(11,'("rvsoadjust    ",a1,f10.5)') parsym(type),
     +        rvsoadjust(type)
            write(11,'("avsoadjust    ",a1,f10.5)') parsym(type),
     +        avsoadjust(type)
            write(11,'("vso1adjust    ",a1,f10.5)') parsym(type),
     +        vso1adjust(type)
            write(11,'("vso2adjust    ",a1,f10.5)') parsym(type),
     +        vso2adjust(type)
            write(11,'("wso1adjust    ",a1,f10.5)') parsym(type),
     +        wso1adjust(type)
            write(11,'("wso2adjust    ",a1,f10.5)') parsym(type),
     +        wso2adjust(type)
            write(11,'("rcadjust      ",a1,f10.5)') parsym(type),
     +        rcadjust(type)
          endif
        enddo
        do type=1,6
          do l=0,4
            write(11,'("tljadjust     ",a1,f10.5,i2)') parsym(type),
     +        tladjust(type,l),l
          enddo
        enddo
        if (k0.eq.1.and.(parinclude(0).or.flagcomp).and.
     +    Rprime.ne.0.) then
          write(11,'("##")')
          write(11,'("## Resonance parameters")')
          write(11,'("## Z   A     S0        R      xs(therm)    D0",
     +      "         a         P        Sn")')
          write(11,'("##",2i4,1p,7e10.3)') Ztarget,Atarget,
     +      Sstrength(0)*1.e4,Rprime,xscaptherm,D0theo(0,0),
     +      alev(0,0),pair(0,0),S(0,0,1)
        endif
      endif
      close (unit=11)
c
c ****************** Integrated binary cross sections ******************
c
c flagomponly: flag to execute ONLY an optical model calculation
c flagexc    : flag for output of excitation functions
c numinc     : number of incident energies
c parname    : name of particle
c
      if (flagomponly) return
      if (flagnatural.or..not.flagexc) return
      write(*,'(/" ########## EXCITATION FUNCTIONS ###########"/)')
      totfile='total.tot'
      open (unit=1,status='old',file=totfile,iostat=istat)
      if (istat.eq.0) then
        write(*,'(" 1. Total (binary) cross sections"/)')
        write(*,'("   Energy    Non-elastic  Elastic     Total  ",
     +    "  Comp. el.  Shape el.  Reaction   Comp non-el",
     +    "  Direct   Pre-equil."/)')
        read(1,'(////)')
        do 10 nen=1,numinc
          read(1,'(a112)',end=10) stringtot
          write(*,'(1x,a112)') stringtot
   10   continue
        close (unit=1)
      endif
      binfile='binary.tot'
      open (unit=1,status='old',file=binfile,iostat=istat)
      if (istat.eq.0) then
        write(*,'(/" 2. Binary non-elastic cross sections",
     +    " (non-exclusive)"/)')
        write(*,'("   Energy   ",7(2x,a8,1x)/)')
     +    (parname(type),type=0,6)
        read(1,'(////)')
        do 20 nen=1,numinc
          read(1,'(a112)',end=20) stringtot
          write(*,'(1x,a112)') stringtot
   20   continue
        close (unit=1)
      endif
c
c ************** Total particle production cross sections **************
c
c parskip    : logical to skip outgoing particle
c flagfission: flag for fission
c
      write(*,'(/" 3. Total particle production cross sections")')
      do 110 type=0,6
        if (parskip(type)) goto 110
        totfile=' prod.tot'
        write(totfile(1:1),'(a1)') parsym(type)
        open (unit=1,status='old',file=totfile,iostat=istat)
        if (istat.eq.0) then
          write(*,'(/1x,a8, " production"/)') parname(type)
          write(*,'(" Energy   Cross section Multiplicity"/)')
          read(1,'(////)')
          do 120 nen=1,numinc
            read(1,'(a80)',end=120) string
            write(*,'(1x,a80)') string
  120     continue
          close (unit=1)
        endif
  110 continue
      if (flagfission) then
        fisfile='fission.tot'
        open (unit=1,status='old',file=fisfile,iostat=istat)
        if (istat.eq.0) then
          write(*,'(/" 3b. Total fission cross sections "/)')
          write(*,'(" Energy   Cross section "/)')
          read(1,'(////)')
          do 130 nen=1,numinc
            read(1,'(a80)',end=130) string
            write(*,'(1x,a80)') string
  130     continue
          close (unit=1)
        endif
      endif
c
c ******************** Residual production cross sections **************
c
c Acomp     : mass number index for compound nucleus
c maxA      : maximal number of nucleons away from the initial
c                  compound nucleus
c Zcomp     : charge number index for compound nucleus
c maxZ      : maximal number of protons away from the initial compound
c             nucleus
c Ncomp     : neutron number index for compound nucleus
c maxN      : maximal number of neutrons away from the initial
c             compound nucleus
c rpexist   : flag for existence of residual production cross section
c ZZ,Z      : charge number of residual nucleus
c AA,A      : mass number of residual nucleus
c Qres      : Q-value for residual nucleus
c Ethresh   : threshold incident energy for residual nucleus
c Nlast     : last discrete level
c rpisoexist: flag for existence of isomeric residual production cross
c             section
c nuc       : symbol of nucleus
c tau       : lifetime of state in seconds
c
      write(*,'(/" 4. Residual production cross sections")')
      do 210 Acomp=0,maxA
        do 210 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 210
          if (.not.rpexist(Zcomp,Ncomp)) goto 210
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          rpfile='rp000000.tot'
          write(rpfile(3:8),'(2i3.3)') Z,A
          open (unit=1,status='old',file=rpfile,iostat=istat)
          if (istat.eq.0) then
            write(*,'(/" Production of Z=",i3,
     +        " A=",i3," (",i3,a2,") - Total"/)') Z,A,A,nuc(Z)
            write(*,'(" Q-value    =",f11.6)') Qres(Zcomp,Ncomp,0)
            write(*,'(" E-threshold=",f11.6/)') Ethresh(Zcomp,Ncomp,0)
            write(*,'(" Energy   Cross section"/)')
            read(1,'(////)')
            do 220 nen=1,numinc
              read(1,'(a80)',end=220) string
              write(*,'(1x,a80)') string
  220       continue
            close (unit=1)
          endif
          do 230 nex=0,Nlast(Zcomp,Ncomp,0)
            if (.not.rpisoexist(Zcomp,Ncomp,nex)) goto 230
            isofile='rp000000.L00'
            write(isofile(3:8),'(2i3.3)') Z,A
            write(isofile(11:12),'(i2.2)') nex
            open (unit=1,status='old',file=isofile,iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Production of Z=",i3," A=",i3," (",i3,a2,
     +        ") - Level",i3,"  (lifetime:",1p,e12.5," sec.)"/)')
     +          Z,A,A,nuc(Z),nex,tau(Zcomp,Ncomp,nex)
              write(*,'(" Q-value    =",f11.6)') Qres(Zcomp,Ncomp,nex)
              write(*,'(" E-threshold=",f11.6/)')
     +          Ethresh(Zcomp,Ncomp,nex)
              write(*,'(" Energy   Cross section Branching ratio"/)')
              read(1,'(////)')
              do 240 nen=1,numinc
                read(1,'(a80)',end=240) string
                write(*,'(1x,a80)') string
  240         continue
              close (unit=1)
            endif
  230     continue
  210 continue
c
c ********************** Fission cross sections ************************
c
c fisexist: flag for existence of fission cross section
c
      if (flagfission) then
        write(*,'(/" 4b. Fission cross sections per fissioning",
     +    " nuclide")')
        do 310 Acomp=0,maxA
          do 310 Zcomp=0,maxZ
            Ncomp=Acomp-Zcomp
            if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 310
            if (.not.fisexist(Zcomp,Ncomp)) goto 310
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            rpfile='rp000000.fis'
            write(rpfile(3:8),'(2i3.3)') Z,A
            open (unit=1,status='old',file=rpfile,iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Fission cross section for Z=",i3,
     +          " A=",i3," (",i3,a2,")"/)') Z,A,A,nuc(Z)
              write(*,'("  Energy     Cross section"/)')
              read(1,'(////)')
              do 320 nen=1,numinc
                read(1,'(a80)',end=320) string
                write(*,'(1x,a80)') string
  320         continue
              close (unit=1)
            endif
  310   continue
      endif
c
c ******************** Reactions to discrete states ********************
c
c flagdisc    : flag for output of discrete state cross sections
c Zindex,Zix  : charge number index for residual nucleus
c Nindex,Nix  : neutron number index for residual nucleus
c flagchannels: flag for exclusive channels calculation
c Ltarget     : excited level of target
c jdis        : spin of level
c cparity     : parity of level (character)
c parlev      : parity of level
c edis        : energy of level
c
      if (flagdisc) then
        do 410 type=0,6
          if (type.eq.k0) then
            reactionstring(type)='Inelastic'
          else
            reactionstring(type)='  ( , )  '
            write(reactionstring(type)(4:4),'(a1)') parsym(k0)
            write(reactionstring(type)(6:6),'(a1)') parsym(type)
          endif
  410   continue
        write(*,'(/" 5. Binary reactions to discrete levels",
     +    " and continuum")')
        Zcomp=0
        Ncomp=0
        do 420 type=0,6
          if (parskip(type)) goto 420
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          Z=ZZ(Zcomp,Ncomp,type)
          A=AA(Zcomp,Ncomp,type)
          if (flagchannels) then
            contfile='  .tot'
            write(contfile(1:2),'(2a1)') parsym(k0),parsym(type)
            open (unit=1,status='old',file=contfile,iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Total exclusive ",a9," scattering")')
     +          reactionstring(type)
              write(*,'(/"  Energy     Total      Discrete   ",
     +          "Continuum     (",a1,",g",a1,")"/)') parsym(k0),
     +          parsym(type)
              read(1,'(////)')
              do 430 nen=1,numinc
                read(1,'(a80)',end=430) string
                write(*,'(1x,a80)') string
  430         continue
              close (unit=1)
            endif
          endif
          do 440 nex=0,Nlast(Zix,Nix,0)
            if (type.eq.k0.and.nex.eq.Ltarget) goto 440
            discfile='  .L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(type)
            write(discfile(5:6),'(i2.2)') nex
            open (unit=1,status='old',file=discfile,iostat=istat)
            if (istat.eq.0) then
              write(*,'(/1x,a9," to level",i3," of ",i3,a2,":",f4.1,a1,
     +          " at",f8.5," MeV"/)') reactionstring(type),nex,A,nuc(Z),
     +          jdis(Zix,Nix,nex),cparity(parlev(Zix,Nix,nex)),
     +          edis(Zix,Nix,nex)
              write(*,'("  Energy     Total      Direct",
     +          "      Compound"/)')
              read(1,'(////)')
              do 450 nen=1,numinc
                read(1,'(a80)',end=450) string
                write(*,'(1x,a80)') string
  450         continue
              close (unit=1)
            endif
  440     continue
          if (flagchannels) then
            contfile='  .con'
            write(contfile(1:2),'(2a1)') parsym(k0),parsym(type)
            open (unit=1,status='old',file=contfile,iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Continuum exclusive ",a9," scattering")')
     +          reactionstring(type)
              write(*,'(/"  Energy     Total      Continuum   ( ",a1,
     +          ",g",a1,")"/)') parsym(k0),parsym(type)
              read(1,'(////)')
              do 460 nen=1,numinc
                read(1,'(a80)',end=460) string
                write(*,'(1x,a80)') string
  460         continue
              close (unit=1)
            endif
          endif
  420   continue
      endif
c
c ******************** Exclusive channels cross sections ***************
c
c npart     : number of particles in outgoing channel
c maxchannel: maximal number of outgoing particles in individual
c             channel description (e.g. this is 3 for (n,2np))
c numia,....: maximal number of ejectile in channel description
c chanopen  : flag to open channel with first non-zero cross section
c idnum     : counter for exclusive channel
c idchannel : identifier for exclusive channel
c reacstring: string for exclusive reaction channel
c Qexcl     : Q-value for exclusive channel
c Ethrexc   : threshold incident energy for exclusive channel
c
c 1. Exclusive cross sections
c
      if (flagchannels) then
        write(*,'(/" 6. Exclusive cross sections")')
        do 510 npart=0,maxchannel
        do 510 ia=0,numia
        do 510 ih=0,numih
        do 510 it=0,numit
        do 510 id=0,numid
        do 510 ip=0,numip
        do 510 in=0,numin
          if (in+ip+id+it+ih+ia.ne.npart) goto 510
          if (.not.chanopen(in,ip,id,it,ih,ia)) goto 510
          ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
          do 520 idc=0,idnum
            if (idchannel(idc).eq.ident) then
              xsfile='xs000000.tot'
              write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
              open (unit=1,status='old',file=xsfile,iostat=istat)
              if (istat.eq.0) then
                write(*,'(/"     Emitted particles     reaction")')
                write(*,'("   n   p   d   t   h   a")')
                write(*,'(6i4,7x,a17)')  in,ip,id,it,ih,ia,
     +            reacstring(idc)
                write(*,'(/" Q-value    =",f11.6)') Qexcl(idc,0)
                write(*,'(" E-threshold=",f11.6/)') Ethrexcl(idc,0)
                write(*,'(" Energy   Cross section Gamma c.s. ",
     +            "c.s./res.prod.cs"/)')
                read(1,'(////)')
                do 530 nen=1,numinc
                  read(1,'(a100)',end=530) stringtot
                  write(*,'(1x,a100)') stringtot
  530           continue
                close (unit=1)
              endif
              Zcomp=ip+id+it+2*ih+2*ia
              Ncomp=in+id+2*it+ih+2*ia
              isofile='xs000000.tot'
              write(isofile(3:8),'(6i1)') in,ip,id,it,ih,ia
              do 540 nex=0,Nlast(Zcomp,Ncomp,0)
                write(isofile(11:12),'(i2.2)') nex
                open (unit=1,status='old',file=isofile,iostat=istat)
                if (istat.eq.0) then
                  write(*,'(/" Level",i3," (lifetime:",
     +              1p,e12.5," sec.)"/)') nex,tau(Zcomp,Ncomp,nex)
                  write(*,'(" Q-value    =",f11.6)') Qexcl(idc,nex)
                  write(*,'(" E-threshold=",f11.6/)') Ethrexcl(idc,nex)
                  write(*,'(" Energy   Cross section  Branching"/)')
                  read(1,'(////)')
                  do 550 nen=1,numinc
                    read(1,'(a80)',end=550) string
                    write(*,'(1x,a80)') string
  550             continue
                  close (unit=1)
                endif
  540         continue
            endif
  520     continue
  510   continue
      endif
c
c 2. Exclusive fission cross sections
c
        if (flagfission) then
          write(*,'(/" 6b. Exclusive fission cross sections")')
          do 610 npart=0,maxchannel
          do 610 ia=0,numia
          do 610 ih=0,numih
          do 610 it=0,numit
          do 610 id=0,numid
          do 610 ip=0,numip
          do 610 in=0,numin
            if (in+ip+id+it+ih+ia.ne.npart) goto 610
            if (.not.chanopen(in,ip,id,it,ih,ia)) goto 610
            ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
            do 620 idc=0,idnum
              if (idchannel(idc).eq.ident) then
                xsfile='xs000000.fis'
                write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
                open (unit=1,status='old',file=xsfile,iostat=istat)
                if (istat.eq.0) then
                  write(*,'(/"     Emitted particles     reaction")')
                  write(*,'("   n   p   d   t   h   a")')
                  write(*,'(6i4,7x,a17)') in,ip,id,it,ih,ia,
     +              reacstring(idc)
                  write(*,'(/" Q-value    =",f11.6)') Qexcl(idc,0)
                  write(*,'(" E-threshold=",f11.6/)') Ethrexcl(idc,0)
                  write(*,'(" Energy   Cross section c.s./res.pr.cs"/)')
                  read(1,'(////)')
                  do 630 nen=1,numinc
                    read(1,'(a80)',end=630) string
                    write(*,'(1x,a80)') string
  630             continue
                  close (unit=1)
                endif
              endif
  620       continue
  610     continue
        endif
c
c ************************* Gamma-ray intensities **********************
c
c flaggamdis: flag for output of discrete gamma-ray intensities
c numlev    : maximum number of included discrete levels
c gamexist  : flag for existence of gamma production cross section
c Egam      : gamma energy
c
      if (flaggamdis) then
        write(*,'(/" 7. Gamma-ray intensities")')
        do 710 Zcomp=0,maxZ
          do 710 Ncomp=0,maxN
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            do 720 i1=0,numlev
              do 720 i2=0,i1
                if (.not.gamexist(Zcomp,Ncomp,i1,i2)) goto 720
                gamfile='gam000000L00L00.tot'
                write(gamfile(4:9),'(2i3.3)') Z,A
                write(gamfile(11:12),'(i2.2)') i1
                write(gamfile(14:15),'(i2.2)') i2
                open (unit=1,status='old',file=gamfile,iostat=istat)
                if (istat.eq.0) then
                  Egam=edis(Zcomp,Ncomp,i1)-edis(Zcomp,Ncomp,i2)
                  write(*,'(/1x,i3,a2,": Initial state",i3," (",f4.1,a1,
     +              " at",f8.4,") ---> Final state",i3," (",f4.1,a1,
     +              " at",f8.4,")  Egamma= ",f8.4/)') A,nuc(Z),i1,
     +              jdis(Zcomp,Ncomp,i1),
     +              cparity(parlev(Zcomp,Ncomp,i1)),
     +              edis(Zcomp,Ncomp,i1),i2,jdis(Zcomp,Ncomp,i2),
     +              cparity(parlev(Zcomp,Ncomp,i2)),
     +              edis(Zcomp,Ncomp,i2),Egam
                  write(*,'(" Energy   Cross section"/)')
                  read(1,'(////)')
                  do 730 nen=1,numinc
                    read(1,'(a80)',end=730) string
                    write(*,'(1x,a80)') string
  730             continue
                  close (unit=1)
                endif
  720       continue
  710   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
