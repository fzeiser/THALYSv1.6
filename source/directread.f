      subroutine directread
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Eric Bauge
c | Date  : April 27, 2013
c | Task  : Read ECIS results for direct cross section
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical lexist
      integer          type,Zix,Nix,NL,i,iS,nS,k,nleg,l,iang,itype
      real             levelenergy,xsdwbatot
      double precision xs,ddl
c
c **************** Read direct, inelastic cross sections ***************
c
c k0              : index of incident particle
c parskip         : logical to skip outgoing particle
c Zindex,Zix      : charge number index for residual nucleus
c Nindex,Nix      : neutron number index for residual nucleus
c Nlast,NL        : last discrete level
c numlev2         : maximum number of levels
c Ltarget         : excited level of target
c deform          : deformation parameter
c eoutdis         : outgoing energy of discrete state reaction
c edis,levelenergy: energy of level
c eninccm         : center-of-mass incident energy in MeV
c parA            : mass number of particle
c xs              : help variable
c xsdirdisc       : direct cross section for discrete state
c dorigin         : origin of direct cross section (Direct or Preeq)
c
      inquire (file='ecis.dirang',exist=lexist)
      if (.not.lexist) return
      open (unit=8,status='unknown',file='ecis.dirang')
      open (unit=9,status='unknown',file='ecis.dirleg')
      open (unit=10,status='unknown',file='ecis.dirin')
      do 10 type=k0,k0
        if (parskip(type)) goto 10
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        NL=Nlast(Zix,Nix,0)
c
c 1. Direct collective states
c
        do 20 i=0,numlev2
          if (i.eq.Ltarget.and.type.eq.k0) goto 20
          if (deform(Zix,Nix,i).eq.0.) goto 20
          if (eoutdis(type,i).le.0.) goto 20
          levelenergy=edis(Zix,Nix,i)
          if (eninccm.le.levelenergy+0.1*parA(type)) goto 20
          read(10,'()')
          read(10,*) xs
          xsdirdisc(type,i)=real(xs)
          if (i.le.NL) dorigin(type,i)='Direct'
c
c ******************* Direct reaction Legendre coefficients ************
c
c We read the Legendre coefficients for the direct component of the
c reaction only, the compound nucleus coefficients are calculated by
c TALYS later on.
c
c nS      : number of states
c nleg    : number of Legendre coefficients
c i       : level
c l       : l-value
c dleg,ddl: direct reaction Legendre coefficient
c

          read(9,'(55x,i5)') nS
          do 110 iS=1,nS
            if (iS.eq.1) then
              read(9,'(5x,i5)')  nleg
              do 120 k=1,nleg
                read(9,'()')
  120         continue
            else
              read(9,'(5x,i5)')  nleg
              do 130 k=1,nleg
                read(9,'(5x,i5,e20.10)') l,ddl
                if (i.le.NL) dleg(type,i,l)=real(ddl)
  130         continue
            endif
  110     continue
c
c ************************ Angular distributions ***********************
c
c nangle  : number of angles
c itype,xs: help variables
c directad: direct angular distribution
c
c We first skip the elastic angular distribution
c
          read(8,'()')
          read(8,'(12x,i3)') nS
          do 210 iang=0,nangle
            do 210 k=1,nS
              read(8,'()')
  210     continue
c
c Read direct angular distribution
c
          read(8,'(12x,i3)') nS
          do 220 iang=0,nangle
            do 220 k=1,nS
              read(8,'(i3,12x,e12.5)',end=220,err=220) itype,xs
              if (itype.eq.0) directad(type,i,iang)=real(xs)
  220     continue
   20   continue
c
c 2. Giant resonance states
c
c Egrcoll: energy of giant resonance
c
        if (type.eq.k0) then
          do 310 l=0,3
            do 310 i=1,2
              if (betagr(l,i).eq.0.) goto 310
              levelenergy=Egrcoll(l,i)
              if (eninccm.le.levelenergy+0.1*parA(type)) goto 310
c
c Giant resonance cross section
c
c xsgrcoll: giant resonance cross section
c
              read(10,'()')
              read(10,*) xs
              xsgrcoll(k0,l,i)=real(xs)
c
c Giant resonance angular distribution
c
c nanglecont: number of angles for continuum
c grcollad  : giant resonance angular distribution
c
              read(8,'()')
              read(8,'(12x,i3)') nS
              do 330 iang=0,nanglecont
                do 330 k=1,nS
                  read(8,'()')
  330         continue
              read(8,'(12x,i3)') nS
              do 340 iang=0,nanglecont
                do 340 k=1,nS
                  read(8,'(i3,12x,e12.5)',err=340,end=340) itype,xs
                  if (itype.eq.0) grcollad(k0,l,i,iang)=real(xs)
  340         continue
  310     continue
        endif
c
c ************* Create total direct inelastic cross section ************
c
c xsdwbatot    : direct DWBA cross section summed over discrete states
c xsdirdisctot : direct cross section summed over discrete states
c xscollconttot: total collective cross section in the continuum
c xsdirdiscsum : total direct cross section
c ecisstatus   : status of ECIS file
c
        xsdwbatot=0.
        do 410 i=0,numlev2
          if (i.eq.0.and.type.eq.k0) goto 410
          if (deform(Zix,Nix,i).ne.0.) then
            if (i.le.NL) then
              xsdwbatot=xsdwbatot+xsdirdisc(type,i)
            else
              xscollconttot=xscollconttot+xsdirdisc(type,i)
            endif
          endif
  410   continue
        xsdirdisctot(type)=xsdirdisctot(type)+xsdwbatot
        xsdirdiscsum=xsdirdiscsum+xsdwbatot
   10 continue
      close (unit=8,status=ecisstatus)
      close (unit=9,status=ecisstatus)
      close (unit=10,status=ecisstatus)
      open (unit=3,status='unknown',file='ecis.dircs')
      close (unit=3,status=ecisstatus)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
