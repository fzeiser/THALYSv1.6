      subroutine dwbaecis
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 27, 2013
c | Task  : ECIS calculations of DWBA for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*13 outfile
      integer      ii,itype,type,Zix,Nix,nen1end,nen1,nen2
      real         QQ
c
c ********************** Set ECIS input parameters *********************
c
c title      : title of ECIS input file
c ecis1,ecis2: 100 input flags ('T' or 'F') for ECIS
c flagrel    : flag for relativistic kinematics
c ncoll      : number of nuclear states
c maxJmsd    : maximal spin for MSD calculation
c iterm      : number of iterations
c npp        : number of optical potentials
c hint       : integration step size h
c rmatch     : matching radius
c
      open (unit=9,status='unknown',file='ecisdwba.inp')
      title='DWBA cross sections for MSD                       '
      ecis1='FFFFFFFFFFFTFFFFFFFFFFFFFFFTFTFFFFFFFFFFFFFFFFFFFF'
      ecis2='FFFFFFFFTFFFFTFFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
      if (flagrel) ecis1(8:8)='T'
      ncoll=maxJmsd+2
      iterm=1
      npp=3
      hint=0.
      rmatch=15.
c
c We use a simple formula to estimate the required number of j-values:
c    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
c and we always take a minimum of njmax=20.
c
c njmax     : maximal number of j-values in ECIS
c Atarget   : mass number of target nucleus
c onethird  : 1/3
c projmass  : mass of projectile
c Einc      : incident energy in MeV
c numl      : maximum l-value (set in talys.cmb)
c betamsd   : deformation parameter
c k0        : index of incident particle
c angbeg    : first angle
c anginc    : angle increment
c nanglecont: number of angles for continuum
c angend    : last angle
c
      njmax=int(2.4*1.25*(real(Atarget)**onethird)*0.22*
     +  sqrt(projmass*Einc))
      njmax=max(njmax,20)
      njmax=min(njmax,numl)
      betamsd=0.02
      if (k0.eq.1) then
        angbeg=0.
      else
        angbeg=0.00001
      endif
      anginc=180./nanglecont
      angend=180.
c
c ***************** Loop over possible one-step reactions **************
c
c flagoutdwba  : flag for output of DWBA cross sections for MSD
c parskip      : logical to skip outgoing particle
c Zindex,Zix   : charge number index for residual nucleus
c Nindex,Nix   : neutron number index for residual nucleus
c QQ           : Q-value
c S            : separation energy per particle
c flagonestep  : flag for continuum one-step direct only
c nen1end      : help variable
c msdbins2     : number of energy points for MSD calculation
c Emsdin       : incident MSD energy
c specmass     : specific mass for target nucleus
c Emsd         : MSD energy grid
c Emsdout      : outgoing MSD energy
c Exmsd        : excitation energy for MSD energy grid
c flagecisdwba : flag for new ECIS calculation for DWBA for MSD
c ecisdwbamac  : subroutine to create ECIS input file for macroscopic
c                DWBA calculation for MSD
c dwbaread     : subroutine to read ECIS results for DWBA for MSD
c dwbaout      : subroutine for output of ECIS results for DWBA for MSD
c dwbaint      : subroutine to interpolate DWBA cross sections for MSD
c onecontinuumA: subroutine for unnormalized one-step direct cross
c                sections for MSD
c onestepA     : subroutine for unnormalized one-step direct cross
c                sections for outgoing energy grid
c
      if (flagoutdwba) write(*,
     +  '(/" ++++++++++ DWBA CROSS SECTIONS FOR MSD ++++++++++")')
c
c *************************** Macroscopic MSD **************************
c
c 2. First exchange one-step reaction for multi-step
c
      do 10 ii=1,2
        if (ii.eq.2) then
          inquire (file='ecis.msdin',exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: The first calculation of a run",
     +        " should always be done with ecissave y and ecisdwba y")')
            stop
          endif
          open (unit=8,status='unknown',file='ecis.msdang')
          open (unit=10,status='unknown',file='ecis.msdin')
        endif
        itype=k0
        do 20 type=1,2
          if (parskip(type)) goto 20
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          QQ=S(0,0,itype)-S(0,0,type)
          if (flagonestep) then
            nen1end=0
          else
            nen1end=msdbins2
          endif
          do 30 nen1=0,nen1end,2
            Emsdin=real(specmass(Zix,Nix,itype)*Emsd(nen1))
            do 40 nen2=nen1+2,msdbins2,2
              Emsdout=Emsd(nen2)
              Exmsd=Emsdin-Emsdout+QQ
              if (Exmsd.lt.0..or.(Exmsd+0.1).ge.Emsdin) goto 40
              if (ii.eq.1.and.flagecisdwba) call ecisdwbamac(itype,type)
              if (ii.eq.2) then
                call dwbaread(nen1,nen2)
                if (flagoutdwba) call dwbaout(itype,type,nen1,nen2)
              endif
   40       continue
   30     continue
          if (ii.eq.2) then
            call dwbaint
            if (.not.flagonestep) call onecontinuumA(itype,type)
            call onestepA(type)
          endif
   20   continue
        if (flagonestep) goto 300
c
c 3. Inelastic one-step reaction for multi-step
c
        do 110 type=1,2
          if (parskip(type)) goto 110
          if (type.eq.k0) goto 110
          itype=type
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          do 120 nen1=0,msdbins2,2
            Emsdin=real(specmass(Zix,Nix,itype)*Emsd(nen1))
            do 130 nen2=nen1+2,msdbins2,2
              Emsdout=Emsd(nen2)
              Exmsd=Emsdin-Emsdout
              if (Exmsd.lt.0..or.(Exmsd+0.1).ge.Emsdin) goto 130
              if (ii.eq.1.and.flagecisdwba) call ecisdwbamac(itype,type)
              if (ii.eq.2) then
                call dwbaread(nen1,nen2)
                if (flagoutdwba) call dwbaout(itype,type,nen1,nen2)
              endif
  130       continue
  120     continue
          if (ii.eq.2) then
            call dwbaint
            call onecontinuumA(itype,type)
          endif
  110   continue
c
c 4. Second exchange one-step reaction for multi-step
c
        do 210 itype=1,2
          if (parskip(itype)) goto 210
          if (itype.eq.k0) goto 210
          type=k0
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          QQ=S(0,0,itype)-S(0,0,type)
          do 220 nen1=0,msdbins2,2
            Emsdin=real(specmass(Zix,Nix,itype)*Emsd(nen1))
            do 230 nen2=nen1+2,msdbins2,2
              Emsdout=Emsd(nen2)
              Exmsd=Emsdin-Emsdout+QQ
              if (Exmsd.lt.0..or.(Exmsd+0.1).ge.Emsdin) goto 230
              if (ii.eq.1.and.flagecisdwba) call ecisdwbamac(itype,type)
              if (ii.eq.2) then
                call dwbaread(nen1,nen2)
                if (flagoutdwba) call dwbaout(itype,type,nen1,nen2)
              endif
  230       continue
  220     continue
          if (ii.eq.2) then
            call dwbaint
            call onecontinuumA(itype,type)
          endif
  210   continue
  300   if (ii.eq.1.and.flagecisdwba) then
          write(9,'("fin")')
          close (unit=9)
c
c **************** ECIS calculation for DWBA for MSD *******************
c
c flagoutecis: flag for output of ECIS results
c outfile    : output file
c nulldev    : null device
c ecist      : subroutine ecis, adapted for TALYS
c ecisstatus : status of ECIS file
c
          if (flagoutecis) then
            outfile='ecisdwba.out '
          else
            outfile=nulldev
          endif
          call ecist('ecisdwba.inp ',outfile,
     +      'ecis.msdcs   ','ecis.msdin   ','null         ',
     +      'ecis.msdang  ','null         ')
          open (unit=9,status='unknown',file='ecisdwba.inp')
          close (unit=9,status=ecisstatus)
        endif
        if (ii.eq.2) then
          close (unit=8,status=ecisstatus)
          close (unit=10,status=ecisstatus)
        endif
  10  continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
