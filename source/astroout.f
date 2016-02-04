      subroutine astroout
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely, Stephane Hilaire, Arjan Koning
c | Date  : June 10, 2012
c | Task  : Output of astrophysical reaction rates
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*13  astrofile
      character*1   targetparity,yesno
      integer       Acomp,Zcomp,Ncomp,i,Z,A,ires,iresprod,iwriterp,
     +              zrespro((numZ+1)*(numN+1)),
     +              arespro((numZ+1)*(numN+1))
      real          rateastrorp(numT,(numZ+1)*(numN+1)),rateastroeps
c
c ************************ Output of reaction rates ********************
c
c Acomp      : mass number index for compound nucleus
c maxA       : maximal number of nucleons away from initial compound
c              nucleus
c Zcomp      : proton number index for compound nucleus
c maxZ       : maximal number of protons away from initial compound
c              nucleus
c Ncomp      : neutron number index for compound nucleus
c maxN       : maximal number of neutrons away from initial compound
c              nucleus
c nTmax      : effective number of temperatures for Maxwellian
c rateastro  : thermonuclear reaction rate factor
c macsastro  : thermonuclear reaction cross section
c flagastrogs: flag for calculation of astrophysics reaction rate with
c              target in ground state only
c ZZ,Z       : charge number of residual nucleus
c AA,A       : mass number of residual nucleus
c nuc        : symbol of nucleus
c T9         : Temperature grid in 10**9 K
c partf      : integrated partition function
c
      write(*,'(/" 8. Thermonuclear reaction rates")')
      do 10 Acomp=0,maxA
        do 20 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 20
          do 30 i=1,nTmax
            if (rateastro(Zcomp,Ncomp,i).gt.0.) goto 40
   30     continue
          goto 20
   40     Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
         if (nTmax.eq.1) then
            write(*,'(/" Reaction rate for Z=",i3," A=",i3," (",i3,a2,
     +        ") at <E>=",f8.5," MeV (Excited States Contribution : ",
     +        a1,")"/)') Z,A,A,nuc(Z),astroE,yesno(.not.flagastrogs)
          else
            write(*,'(/" Reaction rate for Z=",
     +        i3," A=",i3," (",i3,a2,")"/)') Z,A,A,nuc(Z)
          endif
          write(*,'("    T        G(T)        Rate       MACS "/)')
          do 50 i=1,nTmax
            write(*,'(1x,f8.4,1p,3e12.5)') T9(i),partf(i),
     +        rateastro(Zcomp,Ncomp,i),macsastro(Zcomp,Ncomp,i)
   50     continue
          if (flagracap.and.Zcomp.eq.0.and.Acomp.eq.0) then
            write(*,'(/"    T      Rate(Eq)    Rate(DC)  ",
     +        "  MACS(Eq)    MACS(DC)  "/)')
            do 60 i=1,nTmax
              write(*,'(1x,f8.4,1p,4e12.5)') T9(i),
     +          rateastro(Zcomp,Ncomp,i)-rateastroracap(i),
     +          rateastroracap(i),
     +          macsastro(Zcomp,Ncomp,i)-macsastroracap(i),
     +          macsastroracap(i)
   60       continue
          endif
   20   continue
   10 continue
c
c Write results to separate files
c
c rateastroeps: cutoff value
c zresprod    : help variable
c Atarget     : mass number of target nucleus
c Ztarget     : charge number of target nucleus
c Starget     : symbol of target nucleus
c parsym      : symbol of particle
c k0          : index of incident particle
c flagfission : flag for fission
c
      astrofile='astrorate.g'
      open (unit=1,status='unknown',file=astrofile)
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",g) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs)
      else
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",g)")')
     +    Atarget,Starget,parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS")')
      do 110 i=1,nTmax
        write(1,'(f8.4,1p,2e12.5)') T9(i),rateastro(0,0,i),
     +    macsastro(0,0,i)
  110 continue
      close (unit=1)
      astrofile='astrorate.p'
      open (unit=1,status='unknown',file=astrofile)
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",p) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs)
      else
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",p)")')
     +    Atarget,Starget,parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS")')
      do 120 i=1,nTmax
        write(1,'(f8.4,1p,2e12.5)') T9(i),rateastro(1,0,i),
     +    macsastro(1,0,i)
  120 continue
      close (unit=1)
      astrofile='astrorate.a'
      open (unit=1,status='unknown',file=astrofile)
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",a) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs)
      else
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",a)")')
     +  Atarget,Starget,parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS")')
      do 130 i=1,nTmax
        write(1,'(f8.4,1p,2e12.5)') T9(i),rateastro(2,2,i),
     +    macsastro(2,2,i)
  130 continue
      close (unit=1)
      if (flagfission) then
        astrofile='astrorate.f'
        open (unit=1,status='unknown',file=astrofile)
        if (nTmax.eq.1) then
          write(1,'("# Reaction rate for ",i3,a2,"(",a1,",f) at <E>=",
     +      f8.5," MeV (Excited States Contribution : ",a1,")")')
     +      Atarget,Starget,parsym(k0),astroE,
     +      yesno(.not.flagastrogs)
        else
          write(1,'("# Reaction rate for ",i3,a2,"(",a1,",a)")')
     +    Atarget,Starget,parsym(k0)
        endif
        write(1,'("#    T       Rate       MACS")')
        do 140 i=1,nTmax
          write(1,'(f8.4,1p,2e12.5)') T9(i),rateastrofis(i),
     +      macsastrofis(i)
  140   continue
        close (unit=1)
      endif
      rateastroeps=1.e-10
      iresprod=0
      do 150 Acomp=0,maxA
        do 150 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 150
          if (Zcomp.eq.0.and.Ncomp.eq.0) goto 150
          if (Zcomp.eq.0.and.Ncomp.eq.1) goto 150
          if (Zcomp.eq.1.and.Ncomp.eq.0) goto 150
          if (Zcomp.eq.2.and.Ncomp.eq.2) goto 150
          iwriterp=0
          do 160 i=1,nTmax
            if (rateastro(Zcomp,Ncomp,i).gt.rateastroeps) iwriterp=1
  160     continue
          if (iwriterp.eq.1) then
            iresprod=iresprod+1
            do 170 i=1,nTmax
              rateastrorp(i,iresprod)=rateastro(Zcomp,Ncomp,i)
  170       continue
            zrespro(iresprod)=ZZ(Zcomp,Ncomp,0)
            arespro(iresprod)=AA(Zcomp,Ncomp,0)
          endif
  150 continue
      astrofile='astrorate.tot'
      targetparity='+'
      if (targetP.eq.-1) targetparity='-'
      open (unit=1,status='unknown',file=astrofile)
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +    " reactions  Jp(GS)=",f5.1,a1," at <E>=",f8.5,
     +    " MeV (Excited States Contribution : ",a1,")")')
     +    Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +    targetparity,astroE,yesno(.not.flagastrogs)
      else
        write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +    " reactions  Jp(GS)=",f5.1,a1)')
     +    Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +    targetparity
      endif
      write(1,'("     T9      G(T)    (",a1,",g)",i3,a2,"  (",a1,
     +  ",n)",i3,a2,"  (",a1,",p)",i3,a2,"  (",a1,",a)",i3,a2,1x,
     +  "   Fission  ",
     +  200(6x,i3,a2,1x))')
     +  parsym(k0),AA(0,0,0),nuc(ZZ(0,0,0)),
     +  parsym(k0),AA(0,1,0),nuc(ZZ(0,1,0)),
     +  parsym(k0),AA(1,0,0),nuc(ZZ(1,0,0)),
     +  parsym(k0),AA(2,2,0),nuc(ZZ(2,2,0)),
     +  (arespro(ires),nuc(zrespro(ires)),ires=1,iresprod)
      do 180 i=1,nTmax
        write(1,'(f8.4,1p,200e12.5)') T9(i),partf(i),rateastro(0,0,i),
     +    rateastro(0,1,i),rateastro(1,0,i),rateastro(2,2,i),
     +    rateastrofis(i),(rateastrorp(i,ires),ires=1,iresprod)
  180 continue
      close (unit=1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
