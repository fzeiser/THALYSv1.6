      subroutine discreteout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Output of cross sections for discrete states
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*6  discfile,contfile,totfile
      character*9  reactionstring(0:6)
      integer      type,Zcomp,Ncomp,Zix,Nix,nex,nen,NL
c
c ************************* Make reaction string ***********************
c
c reactionstring: reaction string
c k0            : index of incident particle
c parsym        : symbol of particle
c
      do 10 type=0,6
        if (type.eq.k0) then
          reactionstring(type)='Inelastic'
        else
          reactionstring(type)='  ( , )  '
          write(reactionstring(type)(4:4),'(a1)') parsym(k0)
          write(reactionstring(type)(6:6),'(a1)') parsym(type)
        endif
   10 continue
c
c ****************** Cross sections for discrete states ****************
c
c Zcomp        : charge number index for compound nucleus
c Ncomp        : neutron number index for compound nucleus
c parskip      : logical to skip outgoing particle
c xsdisctot    : total cross section summed over discrete states
c Zindex,Zix   : charge number index for residual nucleus
c Nindex,Nix   : neutron number index for residual nucleus
c Nlast        : last discrete level
c Ltarget      : excited level of target
c edis         : energy of level
c eoutdis      : outgoing energy of discrete state reaction
c jdis         : spin of level
c cparity      : parity of level (character)
c parlev       : parity of level
c xsdirdisc    : direct cross section for discrete state
c xscompdisc   : compound cross section for discrete state
c xsdisc       : total cross section for discrete state
c dorigin      : origin of direct cross section (Direct or Preeq)
c xsdirdisctot : direct cross section summed over discrete states
c xscompdisctot: compound cross section summed over discrete states
c xsdircont    : direct cross section for continuum
c xscompcont   : compound cross section for continuum
c xsconttot    : total cross section for continuum
c xsdirect     : total direct cross section
c xscompound   : total compound cross section
c xsbinary     : cross section from initial compound to residual nucleus
c xsngn        : total (projectile,gamma-ejectile) cross section
c flagchannels : flag for exclusive channels calculation
c xsexclcont   : exclusive single channel cross section for continuum
c xsexclusive  : exclusive single channel cross section
c
      Zcomp=0
      Ncomp=0
      write(*,'(/" 5. Binary reactions to discrete levels",
     +  " and continuum")')
      do 110 type=0,6
        if (parskip(type)) goto 110
        if (xsbinary(type).eq.0.) goto 110
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        write(*,'(/1x,a9," cross sections:"/)') reactionstring(type)
        write(*,'(" Inclusive:"/)')
        write(*,'(" Level Energy    E-out     J/P       Direct    ",
     +    "Compound      Total     Origin"/)')
        do 120 nex=0,Nlast(Zix,Nix,0)
          if (type.eq.k0.and.nex.eq.Ltarget) goto 120
          write(*,'(1x,i2,2f10.5,f7.1,a1,3f12.5,4x,a6)') nex,
     +      edis(Zix,Nix,nex),eoutdis(type,nex),jdis(Zix,Nix,nex),
     +      cparity(parlev(Zix,Nix,nex)),xsdirdisc(type,nex),
     +      xscompdisc(type,nex),xsdisc(type,nex),dorigin(type,nex)
  120   continue
        write(*,'(31x,3("   ---------"))')
        write(*,'(" Discrete  ",a9,":",10x,3f12.5)')
     +    reactionstring(type),xsdirdisctot(type),xscompdisctot(type),
     +    xsdisctot(type)
        write(*,'(" Continuum ",a9,":",10x,3f12.5)')
     +    reactionstring(type),xsdircont(type),xscompcont(type),
     +    xsconttot(type)
        write(*,'(31x,3("   ---------"))')
        write(*,'(" Total     ",a9,":",10x,3f12.5/)')
     +    reactionstring(type),xsdirect(type),xscompound(type),
     +    xsbinary(type)
        if (type.ne.0)
     +    write(*,'(" (",a1,",g",a1,") cross section:",f12.5)')
     +      parsym(k0),parsym(type),xsngn(type)
        if (flagchannels) then
          write(*,'(/" Exclusive"/)')
          write(*,'(" Discrete  ",a9,":",f12.5)') reactionstring(type),
     +      xsdisctot(type)
          write(*,'(" Continuum ",a9,":",f12.5)') reactionstring(type),
     +      xsexclcont(type)
          write(*,'(" Total     ",a9,":",f12.5)') reactionstring(type),
     +      xsexclusive(type)
        endif
  110 continue
c
c Write results to separate file
c
c filediscrete: flag for discrete level cross sections on separate file
c numinclow   : number of incident energies below Elow
c Atarget     : mass number of target nucleus
c Starget     : symbol of target nucleus
c Ztarget     : charge number of target nucleus
c Qres        : Q-value for residual nucleus
c Ethresh     : threshold incident energy for residual nucleus
c numinc      : number of incident energies
c eninc,Einc  : incident energy in MeV
c
      do 130 type=0,6
        if (parskip(type)) goto 130
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        do 140 nex=0,Nlast(Zix,Nix,0)
          if (type.eq.k0.and.nex.eq.Ltarget) goto 140
          if (filediscrete(nex)) then
            discfile='  .L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(type)
            write(discfile(5:6),'(i2.2)') nex
            if (nin.eq.numinclow+1) then
              open (unit=1,status='unknown',file=discfile)
              write(1,'("# ",a1," + ",i3,a2,": Discrete ",a9,
     +          " cross section - Level",i3)') parsym(k0),Atarget,
     +          Starget,reactionstring(type),nex
              write(1,'("# Q-value    =",1p,e12.5,0p," Spin=",f5.1,
     +          " Parity= ",a1)') Qres(Zix,Nix,nex),jdis(Zix,Nix,nex),
     +          cparity(parlev(Zix,Nix,nex))
              write(1,'("# E-threshold=",1p,e12.5)')
     +          Ethresh(Zix,Nix,nex)
              write(1,'("# # energies =",i3)') numinc
              write(1,'("#    E         xs        Direct    Compound")')
              do 150 nen=1,numinclow
                write(1,'(1p,e10.3,3e12.5)') eninc(nen),
     +            fxsdisc(nen,type,nex),fxsdirdisc(nen,type,nex),
     +            fxscompdisc(nen,type,nex)
  150         continue
            else
              open (unit=1,status='old',file=discfile)
              do 160 nen=1,nin+4
                read(1,*)
  160         continue
            endif
            write(1,'(1p,e10.3,3e12.5)') Einc,xsdisc(type,nex),
     +        xsdirdisc(type,nex),xscompdisc(type,nex)
            close (unit=1)
          endif
  140   continue
  130 continue
c
c Write continuum cross sections to separate file
c
c filechannels: flag for exclusive channel cross sections on
c               separate file
c
      if (filechannels) then
        do 210 type=0,6
          if (parskip(type)) goto 210
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          if (type.eq.k0) then
            NL=0
          else
            NL=Nlast(Zix,Nix,0)
          endif
          contfile='  .con'
          write(contfile(1:2),'(2a1)') parsym(k0),parsym(type)
          if (nin.eq.numinclow+1) then
            open (unit=1,status='unknown',file=contfile)
            write(1,'("# ",a1," + ",i3,a2,": Continuum ",a9,
     +        " cross section")') parsym(k0),Atarget,Starget,
     +        reactionstring(type)
            write(1,'("# Q-value    =",1p,e12.5)') Qres(Zix,Nix,NL)
            write(1,'("# E-threshold=",1p,e12.5)') Ethresh(Zix,Nix,NL)
            write(1,'("# # energies =",i3)') numinc
            write(1,'("#    E         xs      Continuum      (",a1,
     +        ",g",a1,")")') parsym(k0),parsym(type)
            do 220 nen=1,numinclow
              write(1,'(1p,e10.3,3e12.5)') eninc(nen),
     +          fxsexclcont(nen,type)+fxsngn(nen,type),
     +          fxsexclcont(nen,type),fxsngn(nen,type)
  220       continue
          else
            open (unit=1,status='old',file=contfile)
            do 230 nen=1,nin+4
              read(1,*)
  230       continue
          endif
          write(1,'(1p,e10.3,3e12.5)') Einc,
     +      xsexclcont(type)+xsngn(type),xsexclcont(type),xsngn(type)
          close (unit=1)
  210   continue
      endif
c
c Write cumulated binary cross sections to separate file
c
      if (filechannels) then
        do 310 type=0,6
          if (parskip(type)) goto 310
          totfile='  .tot'
          write(totfile(1:2),'(2a1)') parsym(k0),parsym(type)
          if (nin.eq.numinclow+1) then
            open (unit=1,status='unknown',file=totfile)
            write(1,'("# ",a1," + ",i3,a2,": Total exclusive ",a9,
     +        " cross section")') parsym(k0),Atarget,Starget,
     +        reactionstring(type)
            write(1,'("#         ")')
            write(1,'("#         ")')
            write(1,'("# # energies =",i3)') numinc
            write(1,'("#    E       Total     Discrete     ",
     +        "  Continuum  (",a1,",g",a1,")")') parsym(k0),parsym(type)
            do 320 nen=1,numinclow
              write(1,'(1p,e10.3,4e12.5)') eninc(nen),
     +          fxsexclusive(nen,type),fxsdisctot(nen,type),
     +          fxsexclcont(nen,type),fxsngn(nen,type)
  320       continue
          else
            open (unit=1,status='old',file=totfile)
            do 330 nen=1,nin+4
              read(1,*)
  330       continue
          endif
          write(1,'(1p,e10.3,4e12.5)') Einc,xsexclusive(type),
     +      xsdisctot(type),xsexclcont(type),xsngn(type)
          close (unit=1)
  310   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
