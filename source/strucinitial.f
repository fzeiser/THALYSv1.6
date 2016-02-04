      subroutine strucinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Initialization of arrays for various structure parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*2  phstring1(14)
      character*4  phstring2(72)
      character*90 denfile
      integer      Zix,Nix,Z,N,A,type,i,k,l,irad,nen,nex,ipar,j,in,ip,
     +             id,it,ih,ia,is
c
c *********** Initialization of nuclear structure arrays ***************
c
c Masses and separation energies
c
c Nix         : neutron number index for residual nucleus
c numN        : maximal number of neutrons away from the initial
c               compound nucleus
c Zix         : charge number index for residual nucleus
c numZ        : maximal number of protons away from the initial
c               compound nucleus
c nucmass     : mass of nucleus
c expmass     : experimental mass
c thmass      : theoretical mass
c expmexc     : experimental mass excess
c thmexc      : theoretical mass excess
c dumexc      : theoretical mass excess from Duflo-Zuker formula
c beta4       : deformation parameters
c gsspin      : ground state spin
c gsparity    : ground state parity
c ldparexist  : flag for existence of tabulated level density parameters
c numpar      : number of particles
c specmass    : specific mass for target nucleus
c redumass    : reduced mass
c S           : separation energy per particle
c numen       : maximum number of outgoing energies
c egrid       : energies of basic energy grid in MeV
c Q           : Q value
c coullimit   : energy limit for charged particle OMP calculation
c numelem     : number of elements
c nummass     : number of masses
c lmaxinc     : maximal l-value for transmission coefficients for
c               incident channel
c flagreaction: flag for calculation of nuclear reactions
c
      do 10 Nix=0,numN+4
        do 10 Zix=0,numZ+4
          nucmass(Zix,Nix)=0.
          expmass(Zix,Nix)=0.
          thmass(Zix,Nix)=0.
          expmexc(Zix,Nix)=0.
          thmexc(Zix,Nix)=0.
          dumexc(Zix,Nix)=0.
          beta4(Zix,Nix)=0.
          Z=Zinit-Zix
          N=Ninit-Nix
          A=Z+N
          if (mod(A,2).eq.0) then
            gsspin(Zix,Nix)=0.
          else
            gsspin(Zix,Nix)=0.5
          endif
          gsparity(Zix,Nix)=1
   10 continue
      do 20 Nix=0,numN
        do 20 Zix=0,numZ
          ldparexist(Zix,Nix)=.false.
          do 20 type=0,numpar
            specmass(Zix,Nix,type)=0.
            redumass(Zix,Nix,type)=0.
            S(Zix,Nix,type)=0.
   20 continue
      do 30 nen=0,numen
        egrid(nen)=0.
   30 continue
      do 40 type=0,numpar
        Q(type)=0.
        coullimit(type)=0.
   40 continue
      lmaxinc=0
c
c If no nuclear reaction calculation is requested, we skip a large
c part of this subroutine to speed up the calculation.
c
      if (.not.flagreaction) goto 600
c
c Level and deformation parameters
c
c flaginvecis: logical for calculating inverse channel OMP
c numlev2    : maximum number of levels
c jassign    : flag for assignment of spin
c passign    : flag for assignment of parity
c parlev     : parity of level
c edis       : energy of level
c jdis       : spin of level
c leveltype  : type of level (rotational (R) or vibrational (V))
c indexlevel : level index
c indexcc    : level index for coupled channel
c vibband    : band number of level
c lband      : angular momentum
c Kmag,Kband : magnetic quantum number
c iph,iphonon: phonon (1 or 2)
c deform     : deformation parameter
c defpar     : deformation parameter
c numlev     : maximum number of included discrete levels
c tau        : lifetime of state in seconds
c bassign    : flag for assignment of branching ratio
c conv       : conversion coefficient
c colltype   : type of collectivity (D, V or R)
c deftype    : deformation length (D) or parameter (B)
c nlevmax2   : maximum number of levels
c ndef       : number of collective levels
c nrot       : number of deformation parameters
c Lisomer    : level number of isomer
c Nisomer    : number of isomers for this nuclide
c deltaEx    : excitation energy bin for population arrays
c sfactor    : spin factor
c numrotcc   : number of rotational deformation parameters
c rotpar     : deformation parameters for rotational nucleus
c
      flaginvecis=.true.
      do 110 i=0,numlev2
        do 110 Nix=0,numN
          do 110 Zix=0,numZ
            jassign(Zix,Nix,i)=' '
            passign(Zix,Nix,i)=' '
            parlev(Zix,Nix,i)=1
            edis(Zix,Nix,i)=0.
            jdis(Zix,Nix,i)=0.
            leveltype(Zix,Nix,i)='D'
            indexlevel(Zix,Nix,i)=0
            indexcc(Zix,Nix,i)=0
            vibband(Zix,Nix,i)=0
            lband(Zix,Nix,i)=0
            Kband(Zix,Nix,i)=0
            iphonon(Zix,Nix,i)=0
            deform(Zix,Nix,i)=0.
            defpar(Zix,Nix,i)=0.
  110 continue
      do 120 i=0,numlev
        do 120 Nix=0,numN
          do 120 Zix=0,numZ
            tau(Zix,Nix,i)=0.
  120 continue
      do 130 k=0,numlev
        do 130 i=0,numlev
          do 130 Nix=0,numN
            do 130 Zix=0,numZ
              bassign(Zix,Nix,i,k)=' '
              conv(Zix,Nix,i,k)=0.
  130 continue
      do 140 Nix=0,numN
        do 140 Zix=0,numZ
          colltype(Zix,Nix)=' '
          deftype(Zix,Nix)='B'
          nlevmax2(Zix,Nix)=0
          ndef(Zix,Nix)=0
          nrot(Zix,Nix)=0.
          do 140 i=0,numex
            deltaEx(Zix,Nix,i)=0.
            do 140 ipar=-1,1,2
              do 140 j=0,numJ
                sfactor(Zix,Nix,i,j,ipar)=0.
  140 continue
      do 145 Nix=-1,numN
        do 145 Zix=-1,numZ
          Nisomer(Zix,Nix)=0
          do 147 i=0,numisom
            Lisomer(Zix,Nix,i)=0
  147     continue
  145 continue
      do 150 i=0,numrotcc
        do 150 Nix=0,numN
          do 150 Zix=0,numZ
            rotpar(Zix,Nix,i)=0.
  150 continue
c
c Decay data parameters
c
c lambda: decay rate per isotope
c prate : production rate per isotope
c rtyp  : type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
c Thalf : half life of nuclide in sec.
c Td    : half life per time unit
c
      do 160 Nix=0,numN
        do 160 Zix=-1,numZ
          do 170 is=-1,numisom
            lambda(Zix,Nix,is)=0.
            prate(Zix,Nix,is)=0.
            rtyp(Zix,Nix,is)=0
            Thalf(Zix,Nix,is)=1.e30
            do 180 it=1,5
              Td(Zix,Nix,is,it)=0.
  180       continue
  170     continue
  160 continue
c
c Resonance parameters
c
c gamgamth: theoretical total radiative width
c swaveth : theoretical strength function for s-wave
c dD0     : uncertainty in D0
c dgamgam : uncertainty in gamgam
c
      do 210 Nix=0,numN
        do 210 Zix=0,numZ
          gamgamth(Zix,Nix)=0.
          swaveth(Zix,Nix)=0.
          dD0(Zix,Nix)=0.
          dgamgam(Zix,Nix)=0.
  210 continue
c
c Gamma parameters
c
c numgam    : maximum number of l-values for gamma multipolarity
c ngr       : number of GR
c kgr       : constant for gamma-ray strength function
c qrpaexist : flag for existence of tabulated QRPA strength functions
c numgamqrpa: number of energies for QRPA strength function
c eqrpa     : energy grid for QRPA strength function
c fe1qrpa   : tabulated QRPA strength function
c xsracap   : direct-semidirect radiative capture cross section
c xsracapEM : direct-semidirect radiative capture cross section as
c             function of type
c
      do 310 l=1,numgam
        do 310 irad=0,1
          do 310 Nix=0,numN
            do 310 Zix=0,numZ
              ngr(Zix,Nix,irad,l)=1
              kgr(Zix,Nix,irad,l)=0.
  310 continue
      do 320 Nix=0,numN
        do 320 Zix=0,numZ
          qrpaexist(Zix,Nix)=.false.
  320 continue
      do 330 nen=0,numgamqrpa
        do 340 Nix=0,numN
          do 340 Zix=0,numZ
            eqrpa(Zix,Nix,nen)=0.
            do 340 it=1,numTqrpa
              fe1qrpa(Zix,Nix,nen,it)=0.
  340   continue
  330 continue
      do 350 nen=1,numenin
        xsracap(nen)=0.
        do 360 l=1,numgam
          do 360 irad=0,1
            xsracapEM(nen,irad,l)=0.
  360   continue
  350 continue
c
c Optical model parameters
c
c ompglobal  : flag for use of global optical model
c rc0,rv0,...: optical model parameters
c disp       : flag for dispersive optical model
c jlmexist   : flag for existence of tabulated radial matter density
c normjlm    : JLM potential normalization factors
c numomp     : number of energies on optical model file
c eomp       : energies on optical model file
c vomp       : optical model parameters from file
c numNph     : maximal number of neutrons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c numZph     : maximal number of protons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c wvol       : absorption part of the optical potential averaged over
c              the volume
c flagjlm    : flag for using semi-microscopic JLM OMP
c flagracap  : flag for radiative capture model
c alphaomp   : alpha optical model (1=normal, 2= McFadden-Satchler,
c              3-5= folding potential)
c rhojlmn    : density for neutrons
c rhojlmp    : density for protons
c potjlm     : JLM potential depth values
c radjlm     : radial points for JLM potential
c xstotadjust: total cross section adjustment
c xseladjust : elastic cross section adjustment
c xsnonadjust: nonelastic cross section adjustment
c threshnorm : normalization factor at trheshold
c Rprime     : potential scattering radius
c
      do 410 k=1,numpar
        do 410 Nix=0,numN
          do 410 Zix=0,numZ
            ompglobal(Zix,Nix,k)=.false.
            rc0(Zix,Nix,k)=0.
            rv0(Zix,Nix,k)=0.
            av0(Zix,Nix,k)=0.
            v1(Zix,Nix,k)=0.
            v2(Zix,Nix,k)=0.
            v3(Zix,Nix,k)=0.
            w1(Zix,Nix,k)=0.
            w2(Zix,Nix,k)=0.
            w3(Zix,Nix,k)=0.
            w4(Zix,Nix,k)=0.
            rvd0(Zix,Nix,k)=0.
            avd0(Zix,Nix,k)=0.
            d1(Zix,Nix,k)=0.
            d2(Zix,Nix,k)=0.
            d3(Zix,Nix,k)=0.
            rvso0(Zix,Nix,k)=0.
            avso0(Zix,Nix,k)=0.
            vso1(Zix,Nix,k)=0.
            vso2(Zix,Nix,k)=0.
            wso1(Zix,Nix,k)=0.
            wso2(Zix,Nix,k)=0.
            disp(Zix,Nix,k)=.false.
            jlmexist(Zix,Nix,k)=.false.
            normjlm(Zix,Nix,k)=1.
  410 continue
      do 420 i=1,19
        do 420 nen=0,numomp
          do 420 k=1,numpar
            do 420 Nix=0,numNph
              do 420 Zix=0,numZph
                eomp(Zix,Nix,k,nen)=0.
                vomp(Zix,Nix,k,nen,i)=0.
  420 continue
      do 430 type=1,2
        V0(type)=0.
        Vjoin(type)=0.
        Wjoin(type)=0.
        do 430 nen=-200,10*numen
          wvol(type,nen)=0.
  430 continue
      if (flagjlm.or.flagracap.or.alphaomp.ge.3) then
        do 440 k=1,6
          do 440 nen=1,numjlm
            do 440 Nix=0,numN
              do 440 Zix=0,numZ
                rhojlmp(Zix,Nix,nen,k)=0.
                rhojlmn(Zix,Nix,nen,k)=0.
                potjlm(Zix,Nix,nen,k)=0.
  440   continue
        do 450 nen=1,numjlm
          do 450 Nix=0,numN
            do 450 Zix=0,numZ
              radjlm(Zix,Nix,nen)=0.
  450   continue
      endif
      do 460 nen=1,numen6
        xstotadjust(nen)=0.
        xseladjust(nen)=0.
        xsnonadjust(nen)=0.
  460 continue
      do 470 type=0,numpar
        threshnorm(type)=1.
  470 continue
      Rprime=0.
c
c Fission parameters
c
c flagfission: flag for fission
c nfisbar    : number of fission barrier parameters
c nclass2    : number of sets of class2 states
c numbar     : number of fission barriers
c nfistrhb   : number of head band transition states for barrier
c nfisc2hb   : number of class2 states for barrier
c minertia   : moment of inertia of fission barrier deformation
c fecont     : start of continuum energy
c minertc2   : moment of inertia for class2 states
c nfistrrot  : number of rotational transition states for barrier
c nfisc2rot  : number of rotational class2 states per set
c Emaxclass2 : maximum energy for class2 states
c pfistrhb   : parity of head band transition states
c pfisc2hb   : parity of class2 states
c efistrhb   : energy of head band transition states
c jfistrhb   : spin of head band transition states
c efisc2hb   : energy of class2 states
c jfisc2hb   : spin of class2 states
c pfistrrot  : parity of rotational transition states
c efistrrot  : energy of rotational transition states
c jfistrrot  : spin of rotational transition states
c pfisc2rot  : parity of rotational class2 states
c efisc2rot  : energy of rotational class2 states
c jfisc2rot  : spin of rotational class2 states
c nubarexist : flag for existence of nubar file
c
      do 510 Nix=0,numN
        do 510 Zix=0,numZ
          nfisbar(Zix,Nix)=0
          nclass2(Zix,Nix)=0
  510 continue
      if (flagfission) then
        do 520 i=1,numbar
          do 520 Nix=0,numN
            do 520 Zix=0,numZ
              nfistrhb(Zix,Nix,i)=0
              nfisc2hb(Zix,Nix,i)=0
              minertia(Zix,Nix,i)=0.
              fecont(Zix,Nix,i)=0.
              minertc2(Zix,Nix,i)=0.
              nfistrrot(Zix,Nix,i)=0.
              nfisc2rot(Zix,Nix,i)=0.
              Emaxclass2(Zix,Nix,i)=0.
  520   continue
        do 530 k=0,numlev
          do 530 i=1,numbar
            do 530 Nix=0,numN
              do 530 Zix=0,numZ
                pfistrhb(Zix,Nix,i,k)=1
                pfisc2hb(Zix,Nix,i,k)=1
                efistrhb(Zix,Nix,i,k)=0.
                jfistrhb(Zix,Nix,i,k)=0.
                efisc2hb(Zix,Nix,i,k)=0.
                jfisc2hb(Zix,Nix,i,k)=0.
  530   continue
        do 540 k=0,numrot
          do 540 i=1,numbar
            do 540 Nix=0,numN
              do 540 Zix=0,numZ
                pfistrrot(Zix,Nix,i,k)=1
                efistrrot(Zix,Nix,i,k)=0.
                jfistrrot(Zix,Nix,i,k)=0.
                pfisc2rot(Zix,Nix,i,k)=1
                efisc2rot(Zix,Nix,i,k)=0.
                jfisc2rot(Zix,Nix,i,k)=0.
  540   continue
        do 550 k=1,numbar
          do 560 i=1,numbinfis
            eintfis(i,k)=0.
            do 570 j=0,numJ
              do 570 l=-1,1
                rhofis(i,j,l,k)=0.
  570       continue
  560     continue
  550   continue
      endif
      do 580 k=1,numbeta
        betafis(k)=0.
        vfis(k)=0.
  580 continue
      do 590 k=1,2*numbar
        Vpos(k)=0.
        Vheight(k)=0.
        Vwidth(k)=0.
  590 continue
      nubarexist=.false.
c
c Level density parameters
c
c Nlast      : last discrete level
c Ediscrete  : energy of middle of discrete level region
c scutoffdisc: spin cutoff factor for discrete level region
c delta      : energy shift
c ldexist    : flag for existence of level density table
c edens      : energy grid for tabulated level densities
c ldmodel    : level density model
c nendens    : number of energies for level density grid
c nenphdens  : number of energies for particle-hole state density grid
c Edensmax   : maximum energy on level density table
c Ephdensmax : maximum energy on particle-hole state density table
c ENSDF      : string from original ENSDF discrete level file
c D0theo     : theoretical s-wave resonance spacing
c Econd      : condensation energy
c Ucrit      : critical U
c Scrit      : critical entropy
c Dcrit      : critical determinant
c aldcrit    : critical level density parameter
c ldtable    : level density from table
c ldtottableP: total level density per parity from table
c ldtottable : total level density from table
c rhogrid    : integrated level density
c
  600 do 610 i=0,numbar
        do 610 Nix=0,numN
          do 610 Zix=0,numZ
            Nlast(Zix,Nix,i)=0.
            Ediscrete(Zix,Nix,i)=0.
            scutoffdisc(Zix,Nix,i)=1.
            delta(Zix,Nix,i)=0.
            ldexist(Zix,Nix,i)=.false.
  610 continue
c
c Set energy grid for tabulated level densities
c
      edens(0)=0.
      do 620 nex=1,20
        edens(nex)=0.25*nex
  620 continue
      do 630 nex=21,30
        edens(nex)=5.+0.5*(nex-20)
  630 continue
      do 640 nex=31,40
        edens(nex)=10.+nex-30
  640 continue
      edens(41)=22.5
      edens(42)=25.
      edens(43)=30.
      do 650 nex=44,60
        edens(nex)=30.+10.*(nex-43)
  650 continue
      nenphdens=60
      Ephdensmax=200.
      do 660 i=0,numlev
        do 660 Nix=0,numN
          do 660 Zix=0,numZ
            ENSDF(Zix,Nix,i)='                  '
  660 continue
      do 670 Nix=0,numN
        do 670 Zix=0,numZ
          if (ldmodel(Zix,Nix).le.3) nendens(Zix,Nix)=60
          if (ldmodel(Zix,Nix).eq.4) then
            nendens(Zix,Nix)=55
            Edensmax(Zix,Nix)=150.
          endif
          if (ldmodel(Zix,Nix).ge.5) then
            nendens(Zix,Nix)=60
            Edensmax(Zix,Nix)=200.
          endif
  670 continue
      if (.not.flagreaction) return
      do 710 Nix=0,numN
        do 710 Zix=0,numZ
          D0theo(Zix,Nix)=0.
          if (ldmodel(Zix,Nix).eq.3) then
            Ucrit(Zix,Nix)=0.
            Tcrit(Zix,Nix)=0.
            Econd(Zix,Nix)=0.
            Scrit(Zix,Nix)=0.
            aldcrit(Zix,Nix)=0.
          endif
  710 continue
      do 720 Nix=0,numN
        do 720 Zix=0,numZ
          if (ldmodel(Zix,Nix).ge.4) then
            do 725 i=0,numbar
              do 725 ipar=-1,1,2
                do 725 j=0,29
                  do 725 nex=0,numdens
                    ldtottable(Zix,Nix,nex,i)=0.
                    ldtottableP(Zix,Nix,nex,ipar,i)=0.
                    ldtable(Zix,Nix,nex,j,ipar,i)=0.
  725       continue
          endif
  720 continue
      do 730 ipar=-1,1,2
        do 730 j=0,numJ
          do 730 nex=0,numdens
            do 730 Nix=0,numN
              do 730 Zix=0,numZ
                rhogrid(Zix,Nix,nex,j,ipar)=0.
  730 continue
c
c Weak coupling parameters
c
c jcore: spin of level of core nucleus
c pcore: parity of level of core nucleus
c
      do 740 i=0,numlev2
        do 740 Nix=0,numN
          do 740 Zix=0,numZ
            jcore(Zix,Nix,i)=0.
            pcore(Zix,Nix,i)=1
  740 continue
c
c Particle-hole state densities
c
c phmodel : particle-hole state density model
c phexist2: flag for existence of particle-hole state density table
c phtable2: particle-hole state density from table
c
      if (phmodel.eq.2) then
        do 750 l=0,numexc
          do 750 k=0,numexc
            do 750 j=0,numexc
              do 750 i=0,numexc
                do 755 Nix=0,numN
                  do 755 Zix=0,numZ
                    phexist2(Zix,Nix,i,j,k,l)=.false.
                    phexist1(Zix,Nix,i,j)=.false.
  755           continue
                do 760 Nix=0,numNph
                  do 760 Zix=0,numZph
                    do 765 nex=0,numdens
                      phtable2(Zix,Nix,i,j,k,l,nex)=0.
                      phtable1(Zix,Nix,i,j,nex)=0.
  765               continue
  760           continue
  750   continue

c Configurations for microscopic particle-hole state densities
c
c Nphconf2 : number of 2-component particle-hole configurations
c Nphconf1 : number of 1-component particle-hole configurations
c phstring2: help variable
c ppitable : proton particle number from table
c hpitable : proton hole number from table
c pnutable : neutron particle number from table
c hnutable : neutron hole number from table
c pptable  : particle number from table
c hhtable  : hole number from table
c
        denfile=path(1:lenpath)//'density/ph/z026'
        inquire (file=denfile,exist=lexist)
        if (lexist) then
          Nphconf2=72
          Nphconf1=14
          open (unit=2,status='old',file=denfile)
          read(2,'(/////,10x,72(a4,5x),1x,14(a2,7x))')
     +      (phstring2(i),i=1,72),(phstring1(k),k=1,14)
          do 770 i=1,Nphconf2
            read(phstring2(i),'(4i1)') ppitable(i),hpitable(i),
     +        pnutable(i),hnutable(i)
  770     continue
          do 780 i=1,Nphconf1
            read(phstring1(i),'(2i1)') pptable(i),hhtable(i)
  780     continue
        else
          Nphconf2=0
          Nphconf1=0
        endif
      endif
c
c Giant resonance sum rules
c
c betagr : deformation parameter for giant resonance
c Egrcoll: energy of giant resonance
c Ggrcoll: width of giant resonance
c
      do 810 l=0,3
        do 810 i=1,2
          betagr(l,i)=0.
          Egrcoll(l,i)=0.
          Ggrcoll(l,i)=0.
  810 continue
c
c Q-values and threshold energies
c
c Qres   : Q-value for residual nucleus
c Ethresh: threshold incident energy for residual nucleus
c
      do 910 i=0,numlev
        do 910 Nix=0,numN
          do 910 Zix=0,numZ
            Qres(Zix,Nix,i)=0.
            Ethresh(Zix,Nix,i)=0.
  910 continue
c
c Reaction flags
c
c flagwidth  : flag for width fluctuation calculation
c flagpreeq  : flag for pre-equilibrium calculation
c flagcompang: flag for compound angular distribution calculation
c flaggiant  : flag for collective contribution from giant resonances
c flagmulpre : flag for multiple pre-equilibrium calculation
c
c These flags will be reset later in subroutine energies.f, depending
c on the incident energy.
c
      flagwidth=.false.
      flagpreeq=.false.
      flagcompang=.false.
      flaggiant=.false.
      flagmulpre=.false.
c
c Flags for existence of files
c
c rpexist     : flag for existence of residual production cross section
c fisexist    : flag for existence of fission cross section
c rpisoexist  : flag for existence of isomeric residual production cross
c               section
c gamexist    : flag for existence of gamma production cross section
c numin,....  : maximal number of ejectile in channel description
c chanexist   : flag for existence of exclusive cross section
c gamchanexist: flag for existence of exclusive discrete gamma-rays
c chanfisexist: flag for existence of exclusive fission cross section
c chanisoexist: flag for existence of exclusive isomeric cross section
c fpexist     : flag for existence of fission product
c Nrescue     : number of energies for adjustment factors
c urrexist    : flag for existence of URR
c lminU,lmaxU : minimal and maximal orbital angular momentum
c JminU,JmaxU : minimal and maximal total angular momentum
c flagurrendf : flag for URR info to ENDF
c
      do 1010 Nix=0,numN
        do 1010 Zix=0,numZ
          rpexist(Zix,Nix)=.false.
          fisexist(Zix,Nix)=.false.
          do 1020 nex=0,numlev
            rpisoexist(Zix,Nix,nex)=.false.
 1020     continue
          do 1030 i=0,numlev
            do 1030 j=0,i
              gamexist(Zix,Nix,i,j)=.false.
 1030     continue
 1010 continue
      do 1040 in=0,numin
        do 1040 ip=0,numip
          do 1040 id=0,numid
            do 1040 it=0,numit
              do 1040 ih=0,numih
                do 1040 ia=0,numia
                  chanexist(in,ip,id,it,ih,ia)=.false.
                  gamchanexist(in,ip,id,it,ih,ia)=.false.
                  chanfisexist(in,ip,id,it,ih,ia)=.false.
                  do 1050 nex=0,numlev
                    chanisoexist(in,ip,id,it,ih,ia,nex)=.false.
 1050             continue
 1040 continue
      do 1060 i=1,numelem
        do 1060 j=1,nummass
          fpexist(i,j)=.false.
 1060 continue
      lminU=numl
      lmaxU=0
      flagurrendf=.false.
      do 1070 l=0,numl
        JminU(l)=numJ
        JmaxU(l)=0
        do 1080 type=-1,11
          urrexist(type,l)=.false.
 1080   continue
 1070 continue
      do 1090 i=1,nummt
        do 1100 is=-1,1
          Nrescue(i,is)=0
 1100   continue
 1090 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
