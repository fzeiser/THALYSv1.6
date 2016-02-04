      subroutine msdout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 5, 2006
c | Task  : Output of multi-step direct cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ns,type,nen,iang
c
c ************** Angle-integrated multi-step cross sections ************
c
c parname   : name of particle
c maxmsd    : number of MSD steps
c msdstepint: n-step direct cross section integrated over energy
c msdsum    : multi-step direct cross section summed over steps and
c             integrated over energy
c msdall    : total multi-step direct cross section
c parskip   : logical to skip outgoing particle
c parsym    : symbol of particle
c k0        : index of incident particle
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c egrid     : outgoing energy grid
c msdtot    : multi-step direct cross section summed over steps
c msdstep   : continuum n-step direct cross section
c
      write(*,'(/" ++++++++++ MULTI-STEP DIRECT MODEL ++++++++++")')
      write(*,'(/" 1. Total multi-step direct cross sections"/)')
      write(*,'(" Step ",2(5x,a8)/)') (parname(type),type=1,2)
      do 10 ns=1,maxmsd
        write(*,'(1x,i3,1p,3x,2(1x,e12.5))') ns,(msdstepint(type,ns),
     +    type=1,2)
   10 continue
      write(*,'(/" Total ",1p,2(1x,e12.5))') (msdsum(type),type=1,2)
      write(*,'(/" Total MSD cross section:",f12.5/)') msdall
      write(*,'(" 2. Multi-step direct spectra")')
      do 20 type=1,2
        if (parskip(type)) goto 20
        write(*,'(/," Angle-integrated (",a1,",",a1,") MSD spectra"/)')
     +    parsym(k0),parsym(type)
        write(*,'("  Energy    Total",4(5x,i2,"-step")/)') (ns,ns=1,4)
        do 30 nen=ebegin(type),eend(type)
          write(*,'(1x,f7.3,1p,5e12.5)') egrid(nen),msdtot(type,nen),
     +      (msdstep(type,ns,nen),ns=1,4)
   30   continue
   20 continue
c
c ******************* Multi-step angular distributions *****************
c
c flagddx     : flag for output of double-differential cross sections
c nanglecont  : number of angles for continuum
c anglecont   : angle in degrees for continuum
c msdtotad    : multi-step direct angular distribution summed over steps
c msdstepad   : continuum n-step direct angular distribution
c msdtotintad : multi-step direct angular distribution summed over steps
c               and integrated over energy
c msdstepintad: n-step direct angular distribution integrated over
c               energy
c
      if (.not.flagddx) return
      write(*,'(/" 3. Multi-step direct angular distributions")')
      do 110 type=1,2
        if (parskip(type)) goto 110
        do 120 nen=ebegin(type),eend(type)
          write(*,'(/," (",a1,",",a1,") MSD angular distribution ",
     +      "for E-out= ",f7.3,/)') parsym(k0),parsym(type),egrid(nen)
          write(*,'(" Angle    Total",4(5x,i2,"-step")/)') (ns,ns=1,4)
          do 130 iang=0,nanglecont
            write(*,'(1x,f5.1,1p,5e12.5)') anglecont(iang),
     +        msdtotad(type,nen,iang),(msdstepad(type,ns,nen,iang),
     +        ns=1,4)
  130     continue
  120   continue
        write(*,'(/," (",a1,",",a1,") MSD angular distribution ",
     +    "integrated over energy"/)') parsym(k0),parsym(type)
        write(*,'(" Angle    Total",4(5x,i2,"-step")/)') (ns,ns=1,4)
        do 140 iang=0,nanglecont
          write(*,'(1x,f5.1,1p,5e12.5)') anglecont(iang),
     +      msdtotintad(type,iang),(msdstepintad(type,ns,iang),ns=1,4)
  140   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
