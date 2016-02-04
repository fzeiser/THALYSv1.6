      subroutine decaydata
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 20, 2012
c | Task  : Decay data
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*4  decaychar
      character*80 decayfile,string
      integer      Zix,Z,Nbegin,Nend,Abegin,Aend,N,Nix,ia,is,NC,iline,i
      real         TT,rrt
c
c ******************************** Time units **************************
c
c minutesec: number of seconds in a minute
c hoursec  : number of seconds in an hour
c daysec   : number of seconds in a day
c yearsec  : number of seconds in a year
c
      minutesec=60.
      hoursec=60.*minutesec
      daysec=24.*hoursec
      yearsec=365.25*daysec
c
c ************** Find half lives for all involved nuclides *************
c
c Zix      : charge number index for residual nucleus
c maxZ     : maximal number of protons away from initial compound
c            nucleus
c Z        : charge number of residual nucleus
c Zinit    : charge number of initial compound nucleus
c Nbegin   : first N to be included
c Ninit    : neutron number of initial compound nucleus
c maxN     : maximal number of neutrons away from initial compound
c            nucleus
c Nend     : last N to be included
c Abegin   : first A to be included
c Aend     : last A to be included
c decaychar: help variable
c decayfile: decay data file
c lenpath  : length of pathname
c lambda   : decay rate per isotope
c rtyp     : type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
c Thalf    : half life of nuclide in sec.
c Td       : half life per time unit
c A        : mass number of nucleus
c
c Read decay constants from JEFF-3.1.1 radioactive decay data library
c
      do 10 Zix=-1,maxZ
        Z=Zinit-Zix
        Nbegin=Ninit-maxN
        Nend=Ninit
        Abegin=Z+Nbegin
        Aend=Z+Nend
        decaychar='z   '
        write(decaychar(2:4),'(i3.3)') Z
        decayfile=path(1:lenpath)//'decay/'//decaychar
        open (unit=1,status='unknown',file=decayfile)
  100   read(1,'(a80)',end=200) string
        if (string(72:80).ne.'1451    5') goto 100
        read(string(8:10),'(i3)') ia
        if (ia.lt.Abegin) goto 100
        if (ia.gt.Aend) goto 200
        N=ia-Z
        Nix=Ninit-N
        is=-1
        if (string(11:11).eq.'M') is=1
  110   read(1,'(a80)',end=200) string
        if (string(72:80).ne.'8457    2') goto 110
        read(string(1:11),*) Thalf(Zix,Nix,is)
        read(string(45:55),'(i11)') NC
        iline=1+(NC-1)/6
        do 130 i=1,iline
          read(1,'(a80)',end=200) string
  130   continue
        read(1,'(a80)',end=200) string
        read(1,'(a80)',end=200) string
        read(string(1:11),*) rrt
        rtyp(Zix,Nix,is)=int(rrt)
        if (is.eq.-1) rtyp(Zix,Nix,0)=rtyp(Zix,Nix,-1)
        if (Thalf(Zix,Nix,is).eq.0.) Thalf(Zix,Nix,is)=1.e30
c
c Write half life in years, days, etc.
c
c minutesec: number of seconds in a minute
c hoursec  : number of seconds in an hour
c daysec   : number of seconds in a day
c yearsec  : number of seconds in a year
c
        TT=Thalf(Zix,Nix,is)
        Td(Zix,Nix,is,1)=int(TT/yearsec)
        TT=TT-Td(Zix,Nix,is,1)*yearsec
        Td(Zix,Nix,is,2)=int(TT/daysec)
        TT=TT-Td(Zix,Nix,is,2)*daysec
        Td(Zix,Nix,is,3)=int(TT/hoursec)
        TT=TT-Td(Zix,Nix,is,3)*hoursec
        Td(Zix,Nix,is,4)=int(TT/minutesec)
        TT=TT-Td(Zix,Nix,is,4)*minutesec
        Td(Zix,Nix,is,5)=int(TT)
c
c Calculate decay rates
c Nuclides with a lifetime longer 1.e17 sec are considered stable
c
        if (Thalf(Zix,Nix,is).gt.1.e17) then
          lambda(Zix,Nix,is)=0.
        else
          lambda(Zix,Nix,is)=log(2.)/Thalf(Zix,Nix,is)
        endif
        goto 100
  200   close (unit=1)
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
