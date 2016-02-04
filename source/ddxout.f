      subroutine ddxout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Output of double-differential cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*26 ddxfile
      integer      type,nen,iang,i
      real         Eo(0:numen2),enf,fac,xsa,xsb,xs1,xs2,xs3,xs4,xs5,
     +             Eout,angf
c
c ******************* Double-differential cross sections ***************
c
c ddxmode      : mode for double-differential cross sections: 0: None,
c                1: Angular distributions, 2: Spectra per angle, 3: Both
c parskip      : logical to skip outgoing particle
c xsparticle   : total particle production cross section
c ebegin       : first energy point of energy grid
c eendout      : last energy point of energy grid
c espec,Eo     : outgoing energy grid
c xssumout     : cross section summed over mechanisms
c parname      : name of particle
c nanglecont   : number of angles for continuum
c anglecont    : angle in degrees for continuum
c xssumoutad   : angular distribution summed over mechanisms
c xsdiscoutad  : smoothed angular distribution for discrete state
c xspreeqoutad : preequilibrium angular distribution per particle type
c xsmpreeqoutad: multiple preequilibrium angular distribution
c xscompoutad  : compound emission angular distribution
c ddxecount    : counter for double-differential cross section files
c fileddxe     : designator for double-differential cross sections on
c                separate file: angular distribution
c Starget      : symbol of target nucleus
c natstring    : string extension for file names
c iso          : counter for isotope
c locate       : subroutine to find value in ordered table
c deltaE       : energy bin around outgoing energies
c parsym       : symbol of particle
c Einc         : incident energy in MeV
c xsa,....     : help variables
c
c 1. Angular distributions per outgoing energy
c
      if (ddxmode.eq.1.or.ddxmode.eq.3) then
        write(*,'(/" 9. Double-differential cross sections per",
     +    " outgoing energy")')
        do 10 type=1,6
          if (parskip(type)) goto 10
          if (xsparticle(type).eq.0.) goto 10
          do 20 nen=ebegin(type),eendout(type)
            Eo(nen)=espec(type,nen)
            if (xssumout(type,nen).eq.0.) goto 20
            write(*,'(/" DDX for outgoing ",a8," at ",f7.3," MeV"/)')
     +        parname(type),Eo(nen)
            write(*,'(" Angle   Total      Direct  ",
     +        "   Pre-equil.  Mult. preeq   Compound"/)')
            do 30 iang=0,nanglecont
              write(*,'(1x,f5.1,1p,5e12.5)') anglecont(iang),
     +          xssumoutad(type,nen,iang),xsdiscoutad(type,nen,iang),
     +          xspreeqoutad(type,nen,iang),
     +          xsmpreeqoutad(type,nen,iang),xscompoutad(type,nen,iang)
   30       continue
   20     continue
          do 40 i=1,ddxecount(type)
            enf=fileddxe(type,i)
            call locate(Eo,ebegin(type),eendout(type),enf,nen)
            fac=(enf-Eo(nen))/deltaE(nen)
            ddxfile=' ddxE000.000E000.0.MeV'//natstring(iso)
            write(ddxfile(1:1),'(a1)') parsym(type)
            write(ddxfile(6:12),'(f7.3)') Einc
            write(ddxfile(6:8),'(i3.3)') int(Einc)
            write(ddxfile(14:18),'(f5.1)') enf
            write(ddxfile(14:16),'(i3.3)') int(enf)
            open (unit=1,status='unknown',file=ddxfile)
            write(1,'("# ",a1," + ",i3,a2,": ",a8," DDX spectrum")')
     +        parsym(k0),Atarget,Starget,parname(type)
            write(1,'("# E-incident = ",f7.3)') Einc
            write(1,'("# E-emission = ",f7.3)') enf
            write(1,'("# # angles =",i3)') nanglecont+1
            write(1,'("# Angle    Total       Direct    Pre-equil.",
     +        "  Mult. preeq  Compound")')
            do 50 iang=0,nanglecont
              xsa=xssumoutad(type,nen,iang)
              xsb=xssumoutad(type,nen+1,iang)
              xs1=xsa+fac*(xsb-xsa)
              xsa=xsdiscoutad(type,nen,iang)
              xsb=xsdiscoutad(type,nen+1,iang)
              xs2=xsa+fac*(xsb-xsa)
              xsa=xspreeqoutad(type,nen,iang)
              xsb=xspreeqoutad(type,nen+1,iang)
              xs3=xsa+fac*(xsb-xsa)
              xsa=xsmpreeqoutad(type,nen,iang)
              xsb=xsmpreeqoutad(type,nen+1,iang)
              xs4=xsa+fac*(xsb-xsa)
              xsa=xscompoutad(type,nen,iang)
              xsb=xscompoutad(type,nen+1,iang)
              xs5=xsa+fac*(xsb-xsa)
              write(1,'(f5.1,1p,5e12.5)') anglecont(iang),xs1,xs2,xs3,
     +          xs4,xs5
   50       continue
            close (unit=1)
   40     continue
c
c Results in LAB frame
c
c flagrecoil   : flag for calculation of recoils
c flaglabddx   : flag for calculation of DDX in LAB system
c iejlab       : number of ejectile lab bins
c Eejlab       : center of ejectile lab bin
c ddxejlab     : LAB ddx array
c
          if (flagrecoil.and.flaglabddx) then
            do 60 nen=1,iejlab(type)
              Eo(nen)=Eejlab(type,nen)
              write(*,'(/" DDX for outgoing ",a8," at ",f7.3," MeV",
     +          " in LAB frame")') parname(type),Eo(nen)
              write(*,'(" Angle   Cross section"/)')
              do 70 iang=0,nanglecont
                write(*,'(1x,f5.1,1p,e12.5)') anglecont(iang),
     +            ddxejlab(type,nen,iang)
   70         continue
   60       continue
            do 80 i=1,ddxecount(type)
              enf=fileddxe(type,i)
              call locate(Eo,1,iejlab(type),enf,nen)
              fac=(enf-Eo(nen))/deltaE(nen)
              ddxfile=' ddxE000.000E000.0.lab'//natstring(iso)
              write(ddxfile(1:1),'(a1)') parsym(type)
              write(ddxfile(6:12),'(f7.3)') Einc
              write(ddxfile(6:8),'(i3.3)') int(Einc)
              write(ddxfile(14:18),'(f5.1)') enf
              write(ddxfile(14:16),'(i3.3)') int(enf)
              open (unit=1,status='unknown',file=ddxfile)
              write(1,'("# ",a1," + ",i3,a2,": ",a8," DDX spectrum",
     +          " in LAB system")') parsym(k0),Atarget,Starget,
     +          parname(type)
              write(1,'("# E-incident = ",f7.3)') Einc
              write(1,'("# E-emission = ",f7.3)') enf
              write(1,'("# # angles =",i3)') nanglecont+1
              write(1,'("# Angle    Total")')
              do 90 iang=0,nanglecont
                xsa=ddxejlab(type,nen,iang)
                xsb=ddxejlab(type,nen+1,iang)
                xs1=xsa+fac*(xsb-xsa)
                write(1,'(f5.1,1p,e12.5)') anglecont(iang),xs1
   90         continue
   80       continue
          endif
   10   continue
      endif
c
c 2. Emission spectra per outgoing angle
c
c anginc   : angle increment
c ddxacount: counter for double-differential cross section files
c fileddxa : designator for double-differential cross sections on
c            separate file: spectrum per angle
c
      if (ddxmode.eq.2.or.ddxmode.eq.3) then
        anginc=180./nanglecont
        write(*,'(/" 9. Double-differential cross sections per",
     +    " outgoing angle")')
        do 110 type=1,6
          if (parskip(type)) goto 110
          if (xsparticle(type).eq.0.) goto 110
          do 120 iang=0,nanglecont
            write(*,'(/" DDX for outgoing ",a8," at ",f7.3," degrees")')
     +        parname(type),anglecont(iang)
            write(*,'(/"   E-out    Total      Direct  ",
     +        "   Pre-equil. Mult. preeq   Compound"/)')
            do 130 nen=ebegin(type),eendout(type)
              Eout=espec(type,nen)
              write(*,'(1x,f7.3,1p,5e12.5)') Eout,
     +          xssumoutad(type,nen,iang),xsdiscoutad(type,nen,iang),
     +          xspreeqoutad(type,nen,iang),
     +          xsmpreeqoutad(type,nen,iang),xscompoutad(type,nen,iang)
  130       continue
  120     continue
          do 140 i=1,ddxacount(type)
            angf=fileddxa(type,i)
            call locate(anglecont,0,nanglecont,angf,iang)
            fac=(angf-anglecont(iang))/anginc
            ddxfile=' ddxE000.000A000.0.deg'//natstring(iso)
            write(ddxfile(1:1),'(a1)') parsym(type)
            write(ddxfile(6:12),'(f7.3)') Einc
            write(ddxfile(6:8),'(i3.3)') int(Einc)
            write(ddxfile(14:18),'(f5.1)') angf
            write(ddxfile(14:16),'(i3.3)') int(angf)
            open (unit=1,status='unknown',file=ddxfile)
            write(1,'("# ",a1," + ",i3,a2,": ",a8," DDX spectrum")')
     +        parsym(k0),Atarget,Starget,parname(type)
            write(1,'("# E-incident = ",f7.3)') Einc
            write(1,'("# Angle      = ",f7.3)') angf
            write(1,'("# # energies =",i3)')
     +        eendout(type)-ebegin(type)+1
            write(1,'("# E-out    Total       Direct    Pre-equil.",
     +        "  Mult. preeq  Compound")')
            do 150 nen=ebegin(type),eendout(type)
              Eout=espec(type,nen)
              xsa=xssumoutad(type,nen,iang)
              xsb=xssumoutad(type,nen,iang+1)
              xs1=xsa+fac*(xsb-xsa)
              xsa=xsdiscoutad(type,nen,iang)
              xsb=xsdiscoutad(type,nen,iang+1)
              xs2=xsa+fac*(xsb-xsa)
              xsa=xspreeqoutad(type,nen,iang)
              xsb=xspreeqoutad(type,nen,iang+1)
              xs3=xsa+fac*(xsb-xsa)
              xsa=xsmpreeqoutad(type,nen,iang)
              xsb=xsmpreeqoutad(type,nen,iang+1)
              xs4=xsa+fac*(xsb-xsa)
              xsa=xscompoutad(type,nen,iang)
              xsb=xscompoutad(type,nen,iang+1)
              xs5=xsa+fac*(xsb-xsa)
              write(1,'(f7.3,1p,5e12.5)') Eout,xs1,xs2,xs3,xs4,xs5
  150       continue
            close (unit=1)
  140     continue
c
c Results in LAB frame
c
          if (flagrecoil.and.flaglabddx) then
            do 160 iang=0,nanglecont
              write(*,'(/" DDX for outgoing ",a8," at ",f7.3,
     +          " degrees in LAB frame"/)') parname(type),
     +          anglecont(iang)
              write(*,'(" Energy  Cross section"/)')
              do 170 nen=1,iejlab(type)
                write(*,'(1x,f7.3,1p,e12.5)') Eejlab(type,nen),
     +            ddxejlab(type,nen,iang)
  170         continue
  160       continue
            do 180 i=1,ddxacount(type)
              angf=fileddxa(type,i)
              call locate(anglecont,0,nanglecont,angf,iang)
              fac=(angf-anglecont(iang))/anginc
              ddxfile=' ddxE000.000A000.0.lab'//natstring(iso)
              write(ddxfile(1:1),'(a1)') parsym(type)
              write(ddxfile(6:12),'(f7.3)') Einc
              write(ddxfile(6:8),'(i3.3)') int(Einc)
              write(ddxfile(14:18),'(f5.1)') angf
              write(ddxfile(14:16),'(i3.3)') int(angf)
              open (unit=1,status='unknown',file=ddxfile)
              write(1,'("# ",a1," + ",i3,a2,": ",a8," DDX spectrum",
     +          " in LAB system")') parsym(k0),Atarget,Starget,
     +          parname(type)
              write(1,'("# E-incident = ",f7.3)') Einc
              write(1,'("# Angle      = ",f7.3)') angf
              write(1,'("# # energies =",i3)') iejlab(type)
              write(1,'("# E-out    Total")')
              do 190 nen=1,iejlab(type)
                xsa=ddxejlab(type,nen,iang)
                xsb=ddxejlab(type,nen,iang+1)
                xs1=xsa+fac*(xsb-xsa)
                write(1,'(f7.3,1p,e12.5)') Eejlab(type,nen),xs1
  190         continue
  180       continue
            close (unit=1)
          endif
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
