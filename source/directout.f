      subroutine directout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 7, 2006
c | Task  : Output of direct reaction cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,i,ilev,plev(numlev2),iang,nset,nrest,n1,n2,nen
      real    elev(numlev2),jlev(numlev2),dad(numlev2,0:numang)
c
c *********************** Inelastic cross sections *********************
c
c Zindex,Zix   : charge number index for residual nucleus
c Nindex,Nix   : neutron number index for residual nucleus
c k0           : index of incident particle
c numlev2      : maximum number of levels
c xsdirdisc    : direct cross section for discrete state
c edis,elev    : energy of level
c eoutdis      : outgoing energy of discrete state reaction
c jdis,jlev    : spin of level
c cparity      : parity of level (character)
c parlev,plev  : parity of level
c deftype      : deformation length (D) or parameter (B)
c deform       : deformation parameter
c xsdirdisctot : direct cross section summed over discrete states
c nlev         : number of levels for nucleus
c xscollconttot: total collective cross section in the continuum
c flagang      : flag for output of angular distributions
c nangle       : number of angles
c directad,dad : direct angular distribution
c nset,nrest   : help variables
c angle        : angle
c
      write(*,'(/" ++++++++++ DIRECT CROSS SECTIONS ++++++++++"/)')
      write(*,'(" Direct inelastic cross sections"/)')
      write(*,'(" Level  Energy   E-out      J/P   Cross section",
     +  "  Def. par."/)')
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      do 10 i=1,numlev2
        if (xsdirdisc(k0,i).ne.0.)
     +    write(*,'(1x,i3,2f10.5,f7.1,a1,f12.5,5x,a1,f9.5)')
     +    i,edis(Zix,Nix,i),eoutdis(k0,i),jdis(Zix,Nix,i),
     +    cparity(parlev(Zix,Nix,i)),xsdirdisc(k0,i),deftype(Zix,Nix),
     +    deform(Zix,Nix,i)
   10 continue
      write(*,'(/" Discrete direct inelastic cross section:",
     +  f12.5,"   Level 1-",i3)') xsdirdisctot(k0),nlev(Zix,Nix)
      write(*,'(" Collective cross section in continuum  :",f12.5)')
     +  xscollconttot
      if (flagang) then
        write(*,'(/" Direct inelastic angular distributions")')
        ilev=0
        do 20 i=1,numlev2
          if (xsdirdisc(k0,i).ne.0.) then
            ilev=ilev+1
            elev(ilev)=edis(Zix,Nix,i)
            jlev(ilev)=jdis(Zix,Nix,i)
            plev(ilev)=parlev(Zix,Nix,i)
            do 30 iang=0,nangle
              dad(ilev,iang)=directad(k0,i,iang)
   30       continue
          endif
   20   continue
        nset=ilev/10
        nrest=mod(ilev,10)
        do 40 n1=1,nset
          n2=10*(n1-1)
          write(*,'(/" Angle",10(" Ex=",f6.3,"  "))')
     +      (elev(i),i=n2+1,n2+10)
          write(*,'(4x,10("    JP=",f4.1,a1)/)')
     +      (jlev(i),cparity(plev(i)),i=n2+1,n2+10)
          do 50 iang=0,nangle
            write(*,'(1x,f5.1,1p,10e12.5)') angle(iang),
     +        (dad(i,iang),i=n2+1,n2+10)
   50     continue
   40   continue
        if (nrest.gt.0) then
          write(*,'(/" Angle  ",10("Ex=",f6.3,"   "))')
     +      (elev(i),i=10*nset+1,10*nset+nrest)
          write(*,'(4x,10("    JP=",f4.1,a1)/)') (jlev(i),
     +      cparity(plev(i)),i=10*nset+1,10*nset+nrest)
          do 60 iang=0,nangle
            write(*,'(1x,f5.1,1p,10e12.5)') angle(iang),
     +        (dad(i,iang),i=nset+1,nset+nrest)
   60     continue
        endif
      endif
c
c *********************** Giant resonances *****************************
c
c flaggiant : flag for collective contribution from giant resonances
c xsgrcoll  : giant resonance cross section
c Egrcoll   : energy of giant resonance
c eoutgr    : emission energy
c Ggrcoll   : width of giant resonance
c betagr    : deformation parameter for giant resonance
c xsgrtot   : total smoothed giant resonance cross section
c xscollcont: collective cross section in the continuum
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c anglecont : angle in degrees for continuum
c grcollad  : giant resonance angular distribution
c flagspec  : flag for output of spectra
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c egrid     : outgoing energy grid
c xsgrstate : smoothed giant resonance cross section per state
c xsgr      : smoothed giant resonance cross section
c
      if (.not.flaggiant) return
      write(*,'(/" ++++++++++ GIANT RESONANCES ++++++++++"/)')
      write(*,'("      Cross section   Exc. energy Emis. energy",
     +  "   Width    Deform. par."/)')
      write(*,'(" GMR  :",5f12.5)') xsgrcoll(k0,0,1),Egrcoll(0,1),
     +  eoutgr(k0,0,1),Ggrcoll(0,1),betagr(0,1)
      write(*,'(" GQR  :",5f12.5)') xsgrcoll(k0,2,1),Egrcoll(2,1),
     +  eoutgr(k0,2,1),Ggrcoll(2,1),betagr(2,1)
      write(*,'(" LEOR :",5f12.5)') xsgrcoll(k0,3,1),Egrcoll(3,1),
     +  eoutgr(k0,3,1),Ggrcoll(3,1),betagr(3,1)
      write(*,'(" HEOR :",5f12.5)') xsgrcoll(k0,3,2),Egrcoll(3,2),
     +  eoutgr(k0,3,2),Ggrcoll(3,2),betagr(3,2)
      write(*,'(/" Total:",f12.5/)') xsgrtot(k0)-xscollconttot
      if (flagddx) then
        write(*,'(" Average angular distributions",/)')
        write(*,'(" Angle    GMR         GQR         LEOR      HEOR"/)')
        do 110 iang=0,nanglecont
          write(*,'(1x,f5.1,1p,4e12.5)') anglecont(iang),
     +      grcollad(k0,0,1,iang),grcollad(k0,2,1,iang),
     +      grcollad(k0,3,1,iang),grcollad(k0,3,2,iang)
  110   continue
      endif
      if (flagspec) then
        write(*,'(/" Giant resonance spectra",/)')
        write(*,'("  Energy   Total       GMR        GQR       ",
     +    "LEOR       HEOR     Collective"/)')
        do 120 nen=ebegin(k0),eend(k0)
          write(*,'(1x,f7.3,1p,6e11.4)') egrid(nen),xsgr(k0,nen),
     +      xsgrstate(k0,0,1,nen),xsgrstate(k0,2,1,nen),
     +      xsgrstate(k0,3,1,nen),xsgrstate(k0,3,2,nen),xscollcont(nen)
  120   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
