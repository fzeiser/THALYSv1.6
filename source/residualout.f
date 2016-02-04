      subroutine residualout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 20, 2012
c | Task  : Output of residual production cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*16 rpfile,isofile
      integer      Acomp,Zcomp,Ncomp,nex,Z,A,Ares,nen
c
c *************** Residual production cross sections *******************
c
c Acomp   : mass number index for compound nucleus
c maxA    : maximal number of nucleons away from the initial compound
c           nucleus
c Zcomp   : charge number index for compound nucleus
c maxZ    : maximal number of protons away from the initial compound
c           nucleus
c Ncomp   : neutron number index for compound nucleus
c maxN    : maximal number of neutrons away from the initial compound
c           nucleus
c xspopnuc: population cross section per nucleus
c xseps   : limit for cross sections
c
      write(*,'(/" 4. Residual production cross sections"/)')
      write(*,'("   a. Per isotope"/)')
      write(*,'("   Z   A  nuclide    total     level   ",
     +  "isomeric    isomeric    lifetime")')
      write(*,'(17x,"cross section",7x,"cross section  ratio"/)')
      do 10 Acomp=0,maxA
        do 20 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 20
          if (xspopnuc(Zcomp,Ncomp).lt.xseps) goto 20
c
c A. Total and ground state
c
c ZZ,Z    : charge number of residual nucleus
c AA,A    : mass number of residual nucleus
c nuc     : symbol of nucleus
c xspopex : population cross section summed over spin and parity
c xsbranch: branching ratio for isomeric cross section
c
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          write(*,'(1x,2i4," (",i3,a2,")",1p,e12.5,"    0   ",
     +      1p,e12.5,0p,f9.5)') Z,A,A,nuc(Z),xspopnuc(Zcomp,Ncomp),
     +      xspopex(Zcomp,Ncomp,0),xsbranch(Zcomp,Ncomp,0)
c
c B. Per isomer
c
c Nlast     : last discrete level
c tau       : lifetime of state in seconds
c xsmassprod: residual production cross section per mass unit
c Ares      : mass number of residual nucleus
c
          do 30 nex=1,Nlast(Zcomp,Ncomp,0)
            if (tau(Zcomp,Ncomp,nex).ne.0.) then
              write(*,'(31x,i3,3x,1p,e12.5,0p,f9.5,2x,1p,e12.5,
     +          " sec. ")') nex,xspopex(Zcomp,Ncomp,nex),
     +          xsbranch(Zcomp,Ncomp,nex),tau(Zcomp,Ncomp,nex)
            endif
   30     continue
   20   continue
   10 continue
      write(*,'(/"   b. Per mass"/)')
      write(*,'("   A  cross section"/)')
      do 40 Acomp=0,maxA
        if (xsmassprod(Acomp).gt.xseps) then
          Ares=Ainit-Acomp
          write(*,'(1x,i4,1p,e12.5)') Ares,xsmassprod(Acomp)
        endif
   40 continue
c
c ************* Check of residual production cross section *************
c
c xsresprod  : total residual production (= reaction) cross section
c flagfission: flag for fission
c xsfistot   : total fission cross section
c xsnonel    : non-elastic cross section
c
      write(*,'(/" Total residual production cross section:",f12.5)')
     +  xsresprod
      if (flagfission) then
        write(*,'(" Total fission cross section            :",f12.5)')
     +    xsfistot
        write(*,'(" Fission + res. production cross section:",f12.5)')
     +    xsresprod+xsfistot
      endif
      write(*,'(" Non-elastic cross section              :",f12.5)')
     +  xsnonel
c
c Write results to separate file
c
c fileresidual  : flag for residual production cross sections  on
c                 separate file
c rpexist       : flag for existence of residual production cross
c                 section
c rpisoexist    : flag for existence of isomeric residual production
c                 cross section
c natstring     : string extension for file names
c iso           : counter for isotope
c parsym        : symbol of particle
c k0            : index of incident particle
c Atarget       : mass number of target nucleus
c Ztarget       : charge number of target nucleus
c Qres          : Q-value for residual nucleus
c Ethresh       : threshold incident energy for residual nucleus
c edis          : energy of level
c numinc        : number of incident energies
c numinclow     : number of incident energies below Elow
c nin           : counter for incident energy
c Starget       : symbol of target nucleus
c eninc,Einc    : incident energy in MeV
c flagcompo     : flag for output of cross section components
c xspopdir      : direct population cross section per nucleus
c xspoppreeq    : preequilibrium population cross section per nucleus
c xspopcomp     : compound population cross section per nucleus
c
      if (fileresidual) then
        do 110 Acomp=0,maxA
          do 120 Zcomp=0,maxZ
            Ncomp=Acomp-Zcomp
            if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 120
            if (xspopnuc(Zcomp,Ncomp).lt.xseps.and.
     +        .not.rpexist(Zcomp,Ncomp)) goto 120
c
c A. Total
c
            if (Zcomp.eq.parZ(k0).and.Ncomp.eq.parN(k0)) then
              nex=1
            else
              nex=0
            endif
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            rpfile='rp000000.tot'//natstring(iso)
            write(rpfile(3:8),'(2i3.3)') Z,A
            if (.not.rpexist(Zcomp,Ncomp)) then
              rpexist(Zcomp,Ncomp)=.true.
              open (unit=1,status='unknown',file=rpfile)
              write(1,'("# ",a1," + ",i3,a2,": Production of ",i3,a2,
     +          " - Total")') parsym(k0),Atarget,Starget,A,nuc(Z)
              write(1,'("# Q-value    =",1p,e12.5,0p," mass=",f11.6)')
     +          Qres(Zcomp,Ncomp,nex),nucmass(Zcomp,Ncomp)
              write(1,'("# E-threshold=",1p,e12.5)')
     +          Ethresh(Zcomp,Ncomp,nex)
              write(1,'("# # energies =",i3)') numinc
              if (flagcompo) then
                write(1,'("#    E         xs                    ",
     +            " Direct  Preequilibrium Compound")')
                do 130 nen=1,numinclow
                  write(1,'(1p,e10.3,e12.5,12x,3e12.5)') eninc(nen),
     +              fxspopnuc(nen,Zcomp,Ncomp),0.,0.,
     +              fxspopnuc(nen,Zcomp,Ncomp)
  130           continue
                do 140 nen=numinclow+1,nin-1
                  write(1,'(1p,e10.3,e12.5,12x,3e12.5)') eninc(nen),
     +              0.,0.,0.,0.
  140           continue
              else
                write(1,'("#    E         xs")')
                do 150 nen=1,numinclow
                  write(1,'(1p,e10.3,e12.5)') eninc(nen),
     +              fxspopnuc(nen,Zcomp,Ncomp)
  150           continue
                do 160 nen=numinclow+1,nin-1
                  write(1,'(1p,e10.3,e12.5)') eninc(nen),0.
  160           continue
              endif
            else
              open (unit=1,status='old',file=rpfile)
              do 170 nen=1,nin+4
                read(1,*,end=200,err=200)
  170         continue
            endif
            if (flagcompo) then
              write(1,'(1p,e10.3,e12.5,12x,3e12.5)') Einc,
     +          xspopnuc(Zcomp,Ncomp),xspopdir(Zcomp,Ncomp),
     +          xspoppreeq(Zcomp,Ncomp),xspopcomp(Zcomp,Ncomp)
            else
              write(1,'(1p,e10.3,e12.5)') Einc,xspopnuc(Zcomp,Ncomp)
            endif
  200       close (unit=1)
c
c B. Per ground state and isomer
c
            do 210 nex=1,Nlast(Zcomp,Ncomp,0)
              if (tau(Zcomp,Ncomp,nex).ne.0.) goto 220
  210       continue
            goto 120
  220       do 230 nex=0,Nlast(Zcomp,Ncomp,0)
              if (nex.eq.0.or.tau(Zcomp,Ncomp,nex).ne.0.) then
                isofile='rp000000.L00'//natstring(iso)
                write(isofile(3:8),'(2i3.3)') Z,A
                write(isofile(11:12),'(i2.2)') nex
                if (.not.rpisoexist(Zcomp,Ncomp,nex)) then
                  rpisoexist(Zcomp,Ncomp,nex)=.true.
                  open (unit=1,status='unknown',file=isofile)
                  if (nex.eq.0) then
                    write(1,'("# ",a1," + ",i3,a2,": Production of ",
     +                i3,a2," - Ground state")') parsym(k0),Atarget,
     +                Starget,A,nuc(Z)
                  else
                    write(1,'("# ",a1," + ",i3,a2,": Production of ",
     +                i3,a2," - Level",i3,f12.5," MeV")') parsym(k0),
     +                Atarget,Starget,A,nuc(Z),nex,edis(Zcomp,Ncomp,nex)
                  endif
                  write(1,'("# Q-value    =",1p,e12.5,0p," mass=",
     +              f11.6)') Qres(Zcomp,Ncomp,nex),nucmass(Zcomp,Ncomp)
                  write(1,'("# E-threshold=",1p,e12.5)')
     +              Ethresh(Zcomp,Ncomp,nex)
                  write(1,'("# # energies =",i3)') numinc
                  write(1,'("#    E         xs      Branching")')
                  do 240 nen=1,numinclow
                    write(1,'(1p,e10.3,e12.5,0p,f9.5)') eninc(nen),
     +                fxspopex(nen,Zcomp,Ncomp,nex),
     +                fxsbranch(nen,Zcomp,Ncomp,nex)
  240             continue
                  do 250 nen=numinclow+1,nin-1
                    write(1,'(1p,e10.3,e12.5,0p,f9.5)') eninc(nen),0.,
     +                fxsbranch(max(numinclow,1),Zcomp,Ncomp,nex)
  250             continue
                else
                  open (unit=1,status='old',file=isofile)
                  do 260 nen=1,nin+4
                    read(1,*,end=270,err=270)
  260             continue
                endif
                write(1,'(1p,e10.3,e12.5,0p,f9.5)') Einc,
     +            xspopex(Zcomp,Ncomp,nex),xsbranch(Zcomp,Ncomp,nex)
  270           close (unit=1)
              endif
  230       continue
  120     continue
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
