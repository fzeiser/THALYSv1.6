      subroutine prodres
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 1, 2012
c | Task  : Residual production cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist,flagpositive
      character*16 rpfile
      character*18 nonfile
      character*80 string
      integer      iE,Zix,Nix,Z,N,A,is,nen
      real         E,xs
c
c **************** Read residual production cross sections *************
c
c prodexist   : flag for existence of residual product
c maxZ        : maximal number of protons away from the initial
c               compound nucleus
c Z           : charge number of nucleus
c Zinit       : charge number of initial compound nucleus
c maxN        : maximal number of neutrons away from the initial
c               compound nucleus
c N           : neutron number of nucleus
c Ninit       : neutron number of initial compound nucleus
c numen       : maximum number of energies
c A           : mass number of nucleus
c is          : isotope counter: -1=total, 0=ground state 1=isomer
c Erp         : incident energy
c xsrp        : residual production cross section
c Nenrp       : number of incident energies for residual production
c               cross sections
c nonfile     : file with nonelastic cross sections
c natstring   : string extension for file names
c rpfile      : file with residual production cross sections
c flagpositive: flag for existence of non-zero cross sections
c Estart      : incident energy in MeV
c Eend        : lower end of energy range in MeV
c Ztarget     : charge number of target nucleus
c Atarget     : mass number of target nucleus
c isonuc      : isomeric state of nucleus
c
      do 110 Zix=-1,maxZ
        Z=Zinit-Zix
        do 120 Nix=-1,maxN
          do 130 is=-1,numisom
            prodexist(Zix,Nix,is)=.false.
            Nenrp(Zix,Nix,is)=0
            do 140 nen=1,numen
              Erp(Zix,Nix,is,nen)=0.
              xsrp(Zix,Nix,is,nen)=0.
  140       continue
  130     continue
c
c For convenience in later loops, we store the non-elastic cross
c section in the (-1,-1) element of xsrp. In the next loop, we subtract
c the inelastic cross section from this, i.e. xsrp will contain
c all non-elastic cross sections other than inelastic.
c
          if (Zix.eq.-1.and.Nix.eq.-1) then
            is=-1
            prodexist(Zix,Nix,is)=.true.
            nonfile='nonelastic.tot'//natstring(iso)
            open (unit=1,status='unknown',file=nonfile)
            iE=0
  150       read(1,'(a80)',end=160) string
            if (string(1:1).eq.'#') goto 150
            read(string,*) E,xs
            iE=iE+1
            if (iE.gt.numen) goto 160
            Erp(Zix,Nix,is,iE)=E
            xsrp(Zix,Nix,is,iE)=xs
            goto 150
  160       close (unit=1)
            Nenrp(Zix,Nix,is)=iE
            goto 120
          endif
          if (Nix.eq.-1) goto 120
c
c Residual production cross sections
c
          N=Ninit-Nix
          A=Z+N
          rpfile='rp000000.L00'//natstring(iso)
          write(rpfile(3:5),'(i3.3)') Z
          write(rpfile(6:8),'(i3.3)') A
c
c Check for existence of isomeric cross sections. The total cross
c section is considered for is=-1.
c
c Lisomer: level number of isomer
c Nisomer: number of isomers for this nuclide
c
          do 170 is=-1,Nisomer(Zix,Nix)
            if (is.eq.-1) rpfile(10:12)='tot'
            if (is.eq.0) rpfile(10:12)='L00'
            if (is.gt.0) write(rpfile(11:12),'(i2.2)')
     +        Lisomer(Zix,Nix,is)
            inquire (file=rpfile,exist=lexist)
            if (lexist) then
              prodexist(Zix,Nix,is)=.true.
              flagpositive=.false.
              iE=0
              open (unit=1,status='unknown',file=rpfile)
  180         read(1,'(a80)',end=190) string
              if (string(1:1).eq.'#') goto 180
              read(string,*,err=200) E,xs
              if (E.ge.Eback.and.E.le.Ebeam.and.xs.gt.0.)
     +          flagpositive=.true.
              iE=iE+1
              if (iE.gt.numen) goto 190
              Erp(Zix,Nix,is,iE)=E
              xsrp(Zix,Nix,is,iE)=xs
c
c Subtract inelastic cross section from non-elastic cross section
c
              if (Z.eq.Ztarget.and.A.eq.Atarget.and.is.le.0) then
                xsrp(-1,-1,-1,iE)=xsrp(-1,-1,-1,iE)-xsrp(Zix,Nix,is,iE)
              endif
              goto 180
  190         close (unit=1)
              Nenrp(Zix,Nix,is)=iE
              if (.not.flagpositive) prodexist(Zix,Nix,is)=.false.
            endif
  170     continue
  120   continue
  110 continue
      return
  200 write(*,'(" TALYS-error: Problem in cross section file ",a80)')
     +  rpfile
      stop
      end
Copyright (C) 2010  A.J. Koning
