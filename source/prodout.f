      subroutine prodout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 24, 2013
c | Task  : Output of isotope production
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*3  rstr,ystr
      character*13 state
      character*15 Yfile
      integer      i,Zix,Nix,Z,N,A,is,k,it
c
c ************************* Main output ********************************
c
c Tirrad   : irradiation time per time unit
c Tcool    : cooling time per unit
c Ebeam    : incident energy in MeV for isotope production
c Eback    : lower end of energy range in MeV for isotope production
c Ibeam    : beam current in mA
c rhotarget: target density
c Area     : target area in cm^2
c targetdx : effective thickness of target
c Vtar     : active target volume
c Mtar     : active target mass
c Ntar0    : number of original target atoms
c projnum  : number of incident particles [s^-1]
c heat     : produced heat
c
      write(*,'(/" Summary of isotope production "/)')
      write(*,'(" Maximal irradiation time    : ",i3," years ",i3,
     +  " days ",i3," hours ",i3," minutes ",i3," seconds ")')
     +  (Tirrad(i),i=1,5)
      write(*,'(" Cooling time                : ",i3," years ",i3,
     +  " days ",i3," hours ",i3," minutes ",i3," seconds ")')
     +  (Tcool(i),i=1,5)
      write(*,'(" Energy range                : ",f7.3," --> ",f7.3,
     +  " MeV")') Ebeam,Eback
      write(*,'(" Beam current                : ",f7.3," mA")')
     +  Ibeam
      write(*,'(" Target material density     : ",f7.3," g/cm^3")')
     +  rhotarget
      write(*,'(" Target area                 : ",f7.3," cm^2")')
     +  Area
      write(*,'(" Effective target thickness  : ",f7.3," cm")')
     +  targetdx
      write(*,'(" Effective target volume     : ",f7.3," cm^3")')
     +  Vtar
      write(*,'(" Effective target mass       : ",f7.3," g   ")')
     +  Mtar
      write(*,'(" Number of target atoms      : ",1p,e12.5)')
     +  Ntar0
      write(*,'(" Number of incident particles: ",1p,e12.5," s^-1")')
     +  projnum
      write(*,'(" Produced heat in target     : ",f7.3," kW")')
     +  heat
c
c Output to files per residual product
c
c radiounit: unit for radioactivity: Bq, kBq, MBq, Gbq,
c            mCi, Ci or kCi
c yieldunit: unit for isotope yield: num (number),
c            mug (micro-gram), mg, g, or kg
c maxZ     : maximal number of protons away from the initial
c            compound nucleus
c Z        : charge number of nucleus
c Zinit    : charge number of initial compound nucleus
c maxN     : maximal number of neutrons away from the initial
c            compound nucleus
c N        : neutron number of nucleus
c Ninit    : neutron number of initial compound nucleus
c A        : mass number of nucleus
c Lisomer  : level number of isomer
c Nisomer  : number of isomers for this nuclide
c Yexist   : flag for existence of yield
c Yfile    : file with production yields
c state    : state of final nuclide
c parsym   : symbol of particle
c Atarget  : mass number of target nucleus
c nuc      : symbol of nucleus
c Ztarget  : charge number of target nucleus
c Tp       : irradiation time with maximal yield per time unit
c Td       : half life per time unit
c prate    : production rate per isotope
c lambda   : decay rate per isotope
c Ntime    : number of time points
c Tgrid    : time
c activity : activity of produced isotope in MBq
c yield    : yield of produced isotope in MBq/(mA.h)
c Niso     : number of isotopes produced after irradiation
c Nisorel  : fraction of number of produced isotopes per element
c
      rstr='MBq'
      ystr='   '
      if (radiounit.eq.'bq') rstr=' Bq'
      if (radiounit.eq.'kbq') rstr='KBq'
      if (radiounit.eq.'gbq') rstr='GBq'
      if (radiounit.eq.'ci') rstr=' Ci'
      if (radiounit.eq.'kci') rstr='KCi'
      if (radiounit.eq.'mci') rstr='mCi'
      if (yieldunit.eq.'g') ystr='  g'
      if (yieldunit.eq.'mug') ystr='mug'
      if (yieldunit.eq.'mg') ystr=' mg'
      if (yieldunit.eq.'kg') ystr=' kg'
      do 110 Zix=0,maxZ
        Z=Zinit-Zix
        do 120 Nix=0,maxN
          N=Ninit-Nix
          A=Z+N
          do 130 is=-1,Nisomer(Zix,Nix)
            if (.not.Yexist(Zix,Nix,is)) goto 130
            Yfile='Y000000.tot'//natstring(iso)
            write(Yfile(2:7),'(2i3.3)') Z,A
            if (is.ge.0) Yfile(9:11)='L00'
            if (is.ge.1) then
              write(Yfile(10:11),'(i2.2)') Lisomer(Zix,Nix,is)
              state=' Isomer= x   '
              write(state(10:10),'(i1)') is
            else
              state=' Ground state'
            endif
            open (unit=1,status='unknown',file=Yfile)
            write(1,'("# Reaction: ",a1," + ",i3,a2," Production of ",
     +        i3,a2,a13)') parsym(k0),Atarget,nuc(Ztarget),
     +        A,nuc(Z),state
            if (Thalf(Zix,Nix,is).gt.1.e17) then
              write(1,'("# Half life            : stable")')
            else
              write(1,'("# Half life            : ",i6," years ",i3,
     +          " days",i3," hours",i3," minutes",i3," seconds ")')
     +          (Td(Zix,Nix,is,k),k=1,5)
            endif
            if (Tmax(Zix,Nix,is).gt.1.e17) then
              write(1,'("# Maximum production at: infinity")')
            else
              write(1,'("# Maximum production at: ",i6," years ",i3,
     +          " days",i3," hours",i3," minutes",i3," seconds ")')
     +          (Tp(Zix,Nix,is,k),k=1,5)
            endif
            write(1,'("# Initial production rate: ",1p,e12.5,
     +        " [s^-1] Decay rate: ",e12.5," [s^-1]")')
     +        prate(Zix,Nix,is),lambda(Zix,Nix,is)
            write(1,'("# # time points =",i3)') numtime
            write(1,'("# Time [h] Activity [",a3,"] #isotopes ["a3,"]",
     +        "  Yield [",a3,"/mAh]  Isotopic frac.")') rstr,ystr,rstr
            do 140 it=1,numtime
              write(1,'(f8.1,1p,3e15.5,0p,f15.5)') Tgrid(it),
     +          activity(Zix,Nix,is,it),Niso(Zix,Nix,is,it),
     +          yield(Zix,Nix,is,it),Nisorel(Zix,Nix,is,it)
  140       continue
            close (unit=1)
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C) 2010  A.J. Koning
