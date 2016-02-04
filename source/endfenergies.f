      subroutine endfenergies
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 9, 2013
c | Task  : Energy grid for ENDF-6 file
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          nen,i,k,n,nensub,nenint,type,Zix,Nix,nex,idc,
     +                 nen2
      double precision Ein,degrid,Eeps,ratio,emax,ee,e6tmp
c
c ************************ Basic ENDF-6 energy grid ********************
c
c ompenergyfile: file with energies for OMP calculation (ENDF files
c                only)
c eninclow     : minimal incident energy for nuclear model calculations
c Ein,e6,Eeps  : energies of ENDF-6 energy grid in MeV
c eninclow     : minimal incident energy for nuclear model calculations
c degrid       : energy increment
c enincmax     : maximum incident energy
c Emaxtalys    : maximum acceptable energy for TALYS
c
c The basic ENDF-6 energy grid we use is:
c
c         0.001 -   0.01 MeV  : dE= 0.001 MeV
c         0.01  -   0.1 MeV   : dE= 0.01 MeV
c         0.1   -   1 MeV     : dE= 0.05 MeV
c         1     -   4 MeV     : dE= 0.1 MeV
c         4     -   7 MeV     : dE= 0.2 MeV
c         7     -  20 MeV     : dE= 0.5 MeV
c        20     - 100 MeV     : dE= 1.0 MeV
c       100     - 200 MeV     : dE= 2.0 MeV
c       200     - 300 MeV     : dE= 5.0 MeV
c           above 300 MeV     : dE=10.0 MeV
c
c This grid ensures that the total, elastic and reaction cross section
c are calculated on a sufficiently precise energy grid.
c However, this grid can be overwritten using the 'ompenergyfile'.
c
      e6(1)=eninclow
      nen=1
      if (ompenergyfile(1:1).ne.' ') then
        open (unit=2,status='unknown',file=ompenergyfile)
  10    read(2,*,end=20,err=510) Ein
        if (Ein.le.eninclow) goto 10
        nen=nen+1
        if (nen.gt.numen6) goto 520
        e6(nen)=Ein
        goto 10
  20    close (unit=2)
c
c Sort incident energies in ascending order and remove double points
c
        do 30 i=1,nen
          do 40 k=1,i
            if (e6(i).ge.e6(k)) goto 40
            e6tmp=e6(i)
            e6(i)=e6(k)
            e6(k)=e6tmp
   40     continue
   30   continue
        n=nen
        do 50 i=1,nen-1
          if (e6(i).eq.e6(i+1)) then
            do 60 k=i+1,nen
              e6(k)=e6(k+1)
   60       continue
            n=n-1
          endif
   50   continue
        nen=n
        goto 100
      else
        Ein=0.
        degrid=0.001
   70   Ein=Ein+degrid
        if (Ein.gt.e6(1)) then
          nen=nen+1
          e6(nen)=Ein
        endif
        Eeps=Ein+1.e-4
        if (Eeps.gt.enincmax) goto 100
        if (Eeps.gt.0.01) degrid=0.01
        if (Eeps.gt.0.1) degrid=0.05
        if (Eeps.gt.1.) degrid=0.1
        if (Eeps.gt.4.) degrid=0.2
        if (Eeps.gt.8.) degrid=0.5
        if (Eeps.gt.20.) degrid=1.
        if (Eeps.gt.100.) degrid=2.
        if (Eeps.gt.200.) degrid=5.
        if (Eeps.gt.300.) degrid=10.
        if (Eeps.gt.Emaxtalys) goto 100
        goto 70
      endif
c
c Add possible denser energy point from original energy grid
c
c nensub: intermediate number of incident energies
c numinc: number of incident energies
c
  100 nensub=nen
      if (k0.eq.1) then
        do 110 nin=1,numinc-1
          do 120 nen=1,nensub
            ratio=e6(nen)/eninc(nin)
            if (ratio.ge.0.999.and.ratio.le.1.001) goto 110
  120     continue
          do 130 nen=1,nensub
            if (eninc(nin).gt.e6(nen).and.eninc(nin).le.e6(nen+1)) then
              nenint=nen+1
              goto 140
            endif
  130     continue
          goto 110
  140     nensub=nensub+1
          if (nensub.gt.numen6) goto 520
          do 150 nen=nensub,nenint+1,-1
            e6(nen)=e6(nen-1)
  150     continue
          e6(nenint)=eninc(nin)
  110   continue
        nen=nensub
c
c *************** Add partial thresholds to energy grid ****************
c
c k0        : index of incident particle
c emax      : help variable
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c Nlast     : last discrete level
c Ltarget   : excited level of target
c Ethresh   : threshold incident energy for residual nucleus
c idnum     : counter for exclusive channel
c idchannel : identifier for exclusive channel
c Ethrexcl  : threshold incident energy for exclusive channel
c nen6      : total number of energies
c
        emax=min(enincmax,20.)
        do 210 type=1,6
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          do 220 nex=0,Nlast(Zix,Nix,0)
            if (type.eq.k0.and.nex.eq.Ltarget) goto 220
            ee=Ethresh(Zix,Nix,nex)
            if (ee.gt.eninclow.and.ee.le.emax) then
              nen=nen+1
              if (nen.gt.numen6) goto 520
              e6(nen)=ee
            endif
  220     continue
  210   continue
        do 230 idc=0,idnum
          if (idchannel(idc).eq.100000) goto 230
          if (idchannel(idc).eq.10000) goto 230
          if (idchannel(idc).eq.1000) goto 230
          if (idchannel(idc).eq.100) goto 230
          if (idchannel(idc).eq.10) goto 230
          if (idchannel(idc).eq.1) goto 230
          ee=Ethrexcl(idc,0)
          if (ee.gt.eninclow.and.ee.le.emax) then
            nen=nen+1
            if (nen.gt.numen6) goto 520
            e6(nen)=ee
          endif
  230   continue
      endif
      nen6=nen
c
c *************************** Sort energies ****************************
c
      do 310 nen=1,nen6
        do 310 nen2=nen,nen6
          if (e6(nen).le.e6(nen2)) goto 310
          e6tmp=e6(nen)
          e6(nen)=e6(nen2)
          e6(nen2)=e6tmp
  310 continue
      return
  510 write(*,'(" TALYS-error: Problem in file ",a73)') ompenergyfile
      write(*,'(" after E= ",e12.5)') Ein
      stop
  520 write(*,'(" TALYS-error: there are more than",i4,
     +  " incident energies in file ",a73)') numen6,ompenergyfile
      write(*,'(" numen6 in talys.cmb should be increased")')
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
