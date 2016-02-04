      subroutine grid
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : April 21, 2012
c | Task  : Energy and angle grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical lexist
      integer nen,nen0,type,iang,mt,is
      real    Eout,degrid,Eeps,coulfactor,dang,angval,angmin,angmax,val
c
c ************************ Basic outgoing energy grid ******************
c
c egrid,Eout,Eeps: energies of basic energy grid in MeV
c degrid         : energy increment
c segment        : number of segments to divide emission energy grid
c maxen          : total number of energies
c
c The basic outgoing energy grid we use is:
c
c         0.001, 0.002, 0.005    MeV
c         0.01,  0.02,  0.05     MeV
c         0.1-   2 MeV : dE= 0.1 MeV
c         2  -   4 MeV : dE= 0.2 MeV
c         4  -  20 MeV : dE= 0.5 MeV
c        20  -  40 MeV : dE= 1.0 MeV
c        40  - 200 MeV : dE= 2.0 MeV
c        200 - 300 MeV : dE= 5.0 MeV
c        above 300 MeV : dE=10.0 MeV
c
c This grid ensures that the calculation for outgoing energies of a few
c MeV (around the evaporation peak) is sufficiently precise, whereas at
c higher energies a somewhat coarser energy grid can be used. For the
c same reason, this energy grid is used for the calculation of
c transmission coefficients.
c The grid can be further subdivided with the segment keyword.
c
      Eout=0.
      degrid=0.001
      nen=0
   10 Eout=Eout+degrid
      Eeps=Eout+1.e-4
      if (Eeps.gt.enincmax+12.) goto 20
      if (nen.eq.numen) goto 20
      nen=nen+1
      egrid(nen)=Eout
      if (Eeps.gt.0.002) degrid=0.003
      if (Eeps.gt.0.005) degrid=0.005
      if (Eeps.gt.0.01) degrid=0.01
      if (Eeps.gt.0.02) degrid=0.03
      if (Eeps.gt.0.05) degrid=0.05
      if (Eeps.gt.0.1) degrid=0.1
      if (Eeps.gt.2.) degrid=0.2
      if (Eeps.gt.4.) degrid=0.5/segment
      if (Eeps.gt.20.) degrid=1./segment
      if (Eeps.gt.40.) degrid=2./segment
      if (Eeps.gt.200.) degrid=5./segment
      if (Eeps.gt.300.) degrid=10./segment
      goto 10
   20 maxen=nen
      if (maxen.gt.numen) then
        write(*,'(" TALYS-error: number of energies too large - ",
     +    "increase numen in talys.cmb or lower Emaxtalys")')
        stop
      endif
c
c The widths of the bins around the emission energies are set.
c
c deltaE : energy bin around outgoing energies
c Etop   : top of outgoing energy bin
c Ebottom: bottom of outgoing energy bin
c
      do 30 nen=2,maxen-1
        deltaE(nen)=0.5*(egrid(nen+1)-egrid(nen-1))
        Etop(nen)=0.5*(egrid(nen)+egrid(nen+1))
        Ebottom(nen)=0.5*(egrid(nen)+egrid(nen-1))
   30 continue
      deltaE(1)=egrid(1)+0.5*(egrid(2)-egrid(1))
      Etop(1)=deltaE(1)
      Ebottom(1)=0.
      deltaE(0)=0.
      deltaE(maxen)=0.5*(egrid(maxen)-egrid(maxen-1))
      Etop(maxen)=egrid(maxen)
      Ebottom(maxen)=0.5*(egrid(maxen)+egrid(maxen-1))
c
c ******************** Set lower limit for energy grid *****************
c
c ebegin    : first energy point of energy grid
c Atarget   : mass number of target nucleus
c coulfactor: constant for Coulomb barrier
c coullimit : energy limit for charged particle OMP calculation
c parskip   : logical to skip outgoing particle
c coulbar   : Coulomb barrier
c
c For charged particles it is not necessary, or even numerically
c possible, to calculate transmission coefficients and cross sections
c for very low energies. Therefore, we relate their energy grids to the
c corresponding Coulomb barrier. The first outgoing energy is a factor
c of coulfactor lower than the Coulomb barrier. The last outgoing energy
c of the grid depends on the incident energy and is initialized later
c in the energies subroutine.
c
      ebegin(0)=1
      ebegin(1)=1
      do 110 type=2,6
        ebegin(type)=0
        if (parskip(type)) goto 110
        coulfactor=0.01*(1.+Atarget/200.)
        coullimit(type)=coulfactor*coulbar(type)
        do 120 nen=1,maxen
          if (egrid(nen).gt.coullimit(type)) then
            ebegin(type)=nen
            goto 110
          endif
  120   continue
  110 continue
c
c ********* Set upper limit for energy grid and other energies *********
c
c Einc    : incident energy in MeV
c enincmax: maximum incident energy
c energies: subroutine for energies
c eendmax : last energy point of energy grid for maximum incident energy
c eend    : last energy point of energy grid
c
      Einc=enincmax
      call energies
      do 210 type=0,6
        if (parskip(type)) goto 210
        eendmax(type)=eend(type)
  210 continue
c
c *************** Determine number of low incident energies ************
c
c TALYS performs full nuclear reaction calculations for incident
c energies above Elow only. We keep track of the number of lower
c energies, for which simple empirical cross section estimates are made.
c
c eninclow : minimal incident energy for nuclear model calculations
c D0       : s-wave resonance spacing in keV
c D0theo   : theoretical s-wave resonance spacing
c E1v      : energy at end of 1/v region
c eninc    : incident energy in MeV
c locate   : subroutine to find value in ordered table
c numinc   : number of incident energies
c numinclow: number of incident energies below Elow
c
c For excitation functions that extend to very low incident energies,
c the energy at the end of the 1/v region and the energy at the end
c of the resonance region are also inserted.
c
      if (eninclow.eq.0.) then
        if (D0(0,0).eq.0.) then
          eninclow=min(D0theo(0,0)*1.e-6,1.)
        else
          eninclow=min(D0(0,0)*1.e-6,1.)
        endif
      endif
      if (numinc.ge.numenlow-2) eninclow=min(eninclow,eninc(numenlow-2))
      E1v=0.2*eninclow
      if (eninc(1).lt.E1v) then
        call locate(eninc,1,numinc,E1v,nen0)
        if (eninc(nen0)/E1v.lt.0.99) then
          numinc=numinc+1
          do 310 nen=numinc,nen0+2,-1
            eninc(nen)=eninc(nen-1)
  310     continue
          eninc(nen0+1)=E1v
        endif
        call locate(eninc,1,numinc,eninclow,nen0)
        if (eninc(nen0)/eninclow.lt.0.99) then
          numinc=numinc+1
          do 320 nen=numinc,nen0+2,-1
            eninc(nen)=eninc(nen-1)
  320     continue
          eninc(nen0+1)=eninclow
        endif
      endif
      numinclow=0
      do 330 nen=1,numinc
        if (eninc(nen).lt.eninclow) numinclow=numinclow+1
  330 continue
c
c ************** Set limit for transmission coefficients ***************
c
c translimit: limit for transmission coefficient
c transpower: power for transmission coefficient limit
c
  400 translimit=1./(10**transpower)
c
c **************************** Basic angle grid ************************
c
c 1. Discrete angular distributions
c
c dang  : delta angle
c nangle: number of angles
c angle : angle in degrees
c
      dang=180./nangle
      do 410 iang=0,nangle
        angle(iang)=iang*dang
  410 continue
c
c Angular information necessary for recoil calculation
c
c flagrecoil: flag for the calculation of recoils
c angval    : angle
c angmin    : minimum of angular bin
c angmax    : maximum of angular bin
c deg2rad   : conversion factor for degrees to radians
c cosangmin : cosine of minimum of angular bin
c cosangmax : cosine of maximum of angular bin
c sinangmin : sine of minimum of angular bin
c sinangmax : sine of maximum of angular bin
c dcosang   : width of cosine bin
c
      if (flagrecoil) then
        do 420 iang=0,nangle
          angval=angle(iang)
          angmin=max(0.,angval-0.5*dang)
          angmax=min(180.,angval+0.5*dang)
          cosangmin(iang)=cos(angmin*deg2rad)
          cosangmax(iang)=cos(angmax*deg2rad)
          sinangmin(iang)=sin(angmin*deg2rad)
          sinangmax(iang)=sin(angmax*deg2rad)
          dcosang(iang)=abs(cosangmin(iang)-cosangmax(iang))
  420   continue
        sinangmin(0)=1.e-6
        sinangmax(nangle)=1.e-6
      endif
c
c 2. Continuum angular distributions
c
c nanglecont: number of angles for continuum
c anglecont : angle in degrees for continuum
c
      dang=180./nanglecont
      do 430 iang=0,nanglecont
        anglecont(iang)=iang*dang
  430 continue
c
c Angular information necessary for recoil calculation
c
c angcontmin   : minimum of angular bin
c angcontmax   : maximum of angular bin
c cosangcontmin: cosine of minimum of angular bin
c cosangcontmax: cosine of maximum of angular bin
c sinangcontmin: sine of minimum of angular bin
c sinangcontmax: sine of maximum of angular bin
c dcosangcont  : width of cosine bin
c
      if (flagrecoil) then
        do 440 iang=0,nanglecont
          angval=anglecont(iang)
          angcontmin(iang)=max(0.,angval-0.5*dang)
          angcontmax(iang)=min(180.,angval+0.5*dang)
          cosangcontmin(iang)=cos(angcontmin(iang)*deg2rad)
          cosangcontmax(iang)=cos(angcontmax(iang)*deg2rad)
          sinangcontmin(iang)=sin(angcontmin(iang)*deg2rad)
          sinangcontmax(iang)=sin(angcontmax(iang)*deg2rad)
          dcosangcont(iang)=abs(cosangcontmin(iang)-cosangcontmax(iang))
  440   continue
        do 450 iang=nanglecont+1,2*nanglecont+1
          cosangcontmax(iang)=-cosangcontmax(iang-nanglecont-1)
          cosangcontmin(iang)=-cosangcontmin(iang-nanglecont-1)
          sinangcontmin(iang)=-sinangcontmin(iang-nanglecont-1)
          sinangcontmax(iang)=-sinangcontmax(iang-nanglecont-1)
          dcosangcont(iang)=dcosangcont(2*nanglecont+1-iang)
  450   continue
        sinangcontmin(0)=1.e-6
        sinangcontmax(nanglecont)=1.e-6
        sinangcontmax(2*nanglecont+1)=-1.e-6
        sinangcontmin(nanglecont+1)=-1.e-6
      endif
c
c *** Open files with basic reaction information for incident channel **
c
c flagecissave: flag for saving ECIS input and output files
c ecisstatus  : status of ECIS file
c flaginccalc : flag for new ECIS calculation for incident channel
c flagendf    : flag for information for ENDF-6 file
c flagendfecis: flag for new ECIS calculation for ENDF-6 files
c
      if (flagecissave) then
        ecisstatus='keep'
      else
        ecisstatus='delete'
      endif
      if (.not.flaginccalc) then
        inquire (file='incident.cs',exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: The first calculation of a run",
     +      " should always be done with ecissave y and inccalc y")')
          stop
        endif
        open (unit=13,status='old',file='incident.cs')
        open (unit=17,status='old',file='incident.tr')
        open (unit=18,status='old',file='incident.ang')
        open (unit=19,status='old',file='incident.leg')
        open (unit=20,status='old',file='incident.in')
      endif
      if (flagendf.and..not.flagendfecis) then
        inquire (file='endf.cs',exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: The first calculation of a run",
     +      " should always be done with ecissave y and endfecis y")')
          stop
        endif
        open (unit=23,status='old',file='endf.cs')
      endif
c
c To fit data very precisely, a "normal" TALYS calculation may not
c suffice. Therefore, as a final "rescue" an option is included to read
c incident energy dependent adjustment factors for the most important
c cross sections. The MT number from the ENDF-6 format is used as index.
c
c flagrescue: flag for final rescue: normalization to data
c Crescue   : adjustment factor for this incident energy
c Nrescue   : number of energies for adjustment factors
c Erescue   : energy grid for adjustment factors
c frescue   : adjustment factor
c grescue   : global multiplication factor for incident energy dependent
c             adjustment factors
c rescuefile: file with incident energy dependent adjustment factors
c
      if (flagrescue) then
        do 510 mt=1,nummt
          do 520 is=-1,1
            Crescue(mt,is)=1.
            do 530 nen=1,numen6
              Erescue(mt,is,nen)=0.
              frescue(mt,is,nen)=grescue(mt,is)
  530       continue
            if (rescuefile(mt,is)(1:1).ne.' ') then
              open (unit=2,status='old',file=rescuefile(mt,is))
              nen=1
  540         read(2,*,end=550) Erescue(mt,is,nen),val
              frescue(mt,is,nen)=grescue(mt,is)*val
              if (nen.gt.1.and.(Erescue(mt,is,nen).le.
     +          Erescue(mt,is,nen-1))) then
                write(*,'(" TALYS-error: energies in rescuefile must",
     +            " be given in ascending order")')
                stop
              endif
              nen=nen+1
              if (nen.gt.numen6+1) then
                write(*,'(" TALYS-error: number of lines in",
     +            " rescuefile >",i6)') numen6
                stop
              endif
              goto 540
  550         Nrescue(mt,is)=nen-1
              close (unit=2)
            endif
  520     continue
  510   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
