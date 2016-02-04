      subroutine reacinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Initialization of arrays for various cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,ispin,LL,i,type,iang,Nix,Zix,nen,n1,n2,parity,J,p,nex,
     +        n3,n4,nen1,nen2,nexout,in,ip,id,it,ih,ia
c
c *************** Initialization of primary reaction *******************
c
c primary: flag to designate primary (binary) reaction
c
      primary=.true.
c
c *************** Initialize incident reaction arrays ******************
c
c lmaxinc     : maximal l-value for transmission coefficients for
c               incident channel
c xstotinc    : total cross section (neutrons only) for incident channel
c xsreacinc   : reaction cross section for incident channel
c xsoptinc    : optical model reaction cross section for incident
c               channel
c xsracape    : direct radiative capture cross section
c xselasinc   : total elastic cross section (neutrons only) for
c               incident channel
c xsinitpop   : initial population cross section
c numl        : maximal number of l-values
c Sstrength   : s,p,d,etc-wave strength function
c strengthl   : l neutron strength function
c strengthlj  : (l,j) neutron strength function
c sigurrs     : (l,j) scattering cross section for URR
c sigurrc     : (l,j) capture cross section for URR
c sigurrf     : (l,j) fission cross section for URR
c Tjlinc      : transmission coefficients as a function of spin
c               and l for the incident channel
c Turr        : (l,j) transmission coefficient for URR calculation
c Purrlj      : (l,j) parity for URR calculation
c xsbinarylj  : (l,j) cross section for URR calculation
c nulj        : (l,j) number of degrees of freedom for URR calculation
c urrwidth    : channel width in URR
c Tlinc       : transmission coefficients as a function of l for
c               the incident channel, averaged over spin
c numlev2     : maximum number of included discrete levels
c dleg        : direct reaction Legendre coefficient
c directad    : direct angular distribution
c ruth        : elastic/Rutherford ratio
c dorigin     : origin of direct cross section (Direct or Preeq)
c xsdirdisc   : direct cross section for discrete state
c xsdirdisctot: direct cross section summed over discrete states
c xsdirdiscsum: total direct cross section
c
      lmaxinc=0
      xstotinc=0.
      xsreacinc=0.
      xsracape=0.
      xsracapedisc=0.
      xsracapecont=0.
      xsoptinc=0.
      xselasinc=0.
      xsinitpop=0.
      do 10 l=0,numl
        Sstrength(l)=0.
        strengthl(l)=0.
        do 20 ispin=-1,1
          Tjlinc(ispin,l)=0.
   20   continue
        do 30 J=0,numJ
          strengthlj(l,J)=0.
          sigurrs(l,j)=0.
          sigurrc(l,j)=0.
          sigurrf(l,j)=0.
          Turrljinc(l,J)=0.
          Purrlj(l,J)=1
          do 30 type=-1,numpar
            Turrlj(type,l,J)=0.
            xsbinarylj(type,l,J)=0.
            nulj(type,l,J)=0
            urrwidth(type,l,J)=0.
   30   continue
        Tlinc(l)=0.
        spot(l)=0.
   10 continue
      do 40 LL=0,3*numl
        do 40 i=0,numlev2
          do 40 type=0,numpar
            dleg(type,i,LL)=0.
   40 continue
      do 50 iang=0,numang
        do 60 i=0,numlev2
          do 60 type=0,numpar
            directad(type,i,iang)=0.
   60   continue
        ruth(iang)=0.
   50 continue
      do 70 i=0,numlev2
        do 70 type=0,numpar
          dorigin(type,i)='      '
          xsdirdisc(type,i)=0.
   70 continue
      do 80 type=0,numpar
        xsdirdisctot(type)=0.
   80 continue
      xsdirdiscsum=0.
c
c ********************** Initialization of energies ********************
c
c Nix    : neutron number index for residual nucleus
c numN   : maximal number of neutrons away from the initial
c          compound nucleus
c Zix    : charge number index for residual nucleus
c numZ   : maximal number of protons away from the initial
c          compound nucleus
c Exmax0 : maximum excitation energy for residual nucleus (including
c          negative energies)
c Exmax  : maximum excitation energy for residual nucleus
c maxex  : maximum excitation energy bin for compound nucleus
c numex  : maximal number of excitation energies
c Ex     : excitation energy
c Etotal : total energy of compound system (target + projectile)
c Exinc  : excitation energy of entrance bin
c deltaEx: excitation energy bin for population arrays
c
      do 90 Nix=0,numN
        do 90 Zix=0,numZ
          Exmax0(Zix,Nix)=0.
          Exmax(Zix,Nix)=0.
          maxex(Zix,Nix)=0.
          do 100 i=0,numex
            Ex(Zix,Nix,i)=0.
  100     continue
   90 continue
      Exmax0(0,0)=Etotal
      Exmax(0,0)=Etotal
      Exinc=Etotal
      do 105 i=0,numex
        deltaEx(0,0,i)=0.
  105 continue
c
c ****************** Initialize giant resonance arrays *****************
c
c xscollconttot: total collective cross section in the continuum
c xsgrcoll     : giant resonance cross section
c eoutgr       : outgoing energy
c flagddx      : flag for output of double-differential cross sections
c numangcont   : maximum number of angles for continuum
c grcollad     : giant resonance angular distribution
c numen        : maximum number of outgoing energies
c xsgrstate    : smoothed giant resonance cross section per state
c xsgr         : smoothed giant resonance cross section
c xsgrad       : smoothed giant resonance angular distribution
c xsgrtot      : total smoothed giant resonance cross section
c collcontad   : collective angular distribution in the continuum
c xsgrsum      : sum over giant resonance cross sections
c xscollcont   : collective cross section in the continuum
c
      xscollconttot=0.
      do 110 i=1,2
        do 110 l=0,3
          do 110 type=0,numpar
            xsgrcoll(type,l,i)=0.
            eoutgr(type,l,i)=0.
  110 continue
      if (flagddx) then
        do 120 iang=0,numangcont
          do 120 i=1,2
            do 120 l=0,3
              do 120 type=0,numpar
                grcollad(type,l,i,iang)=0.
  120   continue
      endif
      do 130 nen=0,numen
        do 130 i=1,2
          do 130 l=0,3
            do 130 type=0,numpar
              xsgrstate(type,l,i,nen)=0.
  130 continue
      do 140 nen=0,numen
        do 140 type=0,numpar
          xsgr(type,nen)=0.
  140 continue
      do 150 iang=0,numangcont
        do 150 nen=0,numen
          do 150 type=0,numpar
            xsgrad(type,nen,iang)=0.
  150 continue
      do 160 type=0,numpar
        xsgrtot(type)=0.
  160 continue
      do 170 iang=0,numangcont
        do 170 nen=0,numen
          collcontad(nen,iang)=0.
  170 continue
      xsgrsum=0.
      do 180 nen=0,numen
        xscollcont(nen)=0.
  180 continue
c
c *************** Initialize pre-equilibrium arrays ********************
c
c numparx       : maximum number of particles
c xsstep        : preequilibrium cross section per particle type, stage
c                 and outgoing energy
c xsstep2       : two-component preequilibrium cross section
c wemission     : emission rate
c wemission2    : two-component emission rate
c xspreeq       : preequilibrium cross section per particle type and
c                 outgoing energy
c xsmpreeq      : multiple pre-equilibrium emission spectrum
c xspreeqps     : preequilibrium cross section per particle type and
c                 outgoing energy for pickup and stripping
c xspreeqki     : preequilibrium cross section per particle type and
c                 outgoing energy for knockout and inelastic
c xspreeqbu     : preequilibrium cross section per particle type and
c                 outgoing energy for breakup
c parity        : parity
c J             : total angular momentum
c numJph        : maximal spin for particle-hole states
c xspreeqJP     :   preequilibrium cross section per particle type,
c                 outgoing energy, spin and parity
c xspreeqtot    : preequilibrium cross section per particle type
c xspreeqtotps  : preequilibrium cross section per particle type for
c                 pickup and stripping
c xspreeqtotki  : preequilibrium cross section per particle type for
c                 knockout and inelastic
c xspreeqtotbu  : preequilibrium cross section per particle type for
c                 breakup
c xssteptot     : preequilibrium cross section per particle type and
c                 stage
c numJ          : maximal J-value
c xscompad      : compound emission angular distribution
c xsbinemisad   : angular distribution for emission from first compound
c                 nucleus
c flagang       : flag for output of angular distributions
c xspreeqad     : preequilibrium angular distribution per particle type
c xsmpreeqad    : multiple preequilibrium angular distribution
c xspreeqdisc   : preequilibrium cross section for discrete state
c xspreeqdisctot: preequilibrium cross section summed over discrete
c                 states
c xspreeqdiscsum: total preequilibrium cross section for discrete states
c flagmulpre    : flag for multiple pre-equilibrium calculation
c numNph        : maximal number of neutrons away from the initial
c                 compound nucleus for multiple pre-equilibrium emission
c numZph        : maximal number of protons away from the initial
c                 compound nucleus for multiple pre-equilibrium emission
c xspopph       : population cross section per particle-hole
c                 configuration
c PP2           : total strength
c Spre          : time-integrated strength of two-component exciton
c                 state
c xspopph2      : population cross section per two-component
c                 particle-hole configuration
c xspreeqsum    : total preequilibrium cross section summed over
c                 particles
c Esurf         : well depth for surface interaction
c
      do 210 nen=0,numen
        do 210 n1=1,numparx
          do 210 type=0,numpar
          xsstep(type,n1,nen)=0.
  210 continue
      do 220 nen=0,numen
        do 220 n2=0,numparx
          do 220 n1=0,numparx
            do 220 type=0,numpar
              xsstep2(type,n1,n2,nen)=0.
              wemission(type,n1,n2,nen)=0.
              do 220 n4=0,numparx
                do 220 n3=0,numparx
                  wemission2(type,n1,n2,n3,n4,nen)=0.
  220 continue
      do 230 nen=0,numen
        do 230 type=0,numpar
          xspreeq(type,nen)=0.
          xsmpreeq(type,nen)=0.
          xspreeqps(type,nen)=0.
          xspreeqki(type,nen)=0.
          xspreeqbu(type,nen)=0.
  230 continue
      do 240 parity=-1,1,2
        do 240 J=0,numJph
          do 240 nen=0,numen
            do 240 type=0,numpar
              xspreeqJP(type,nen,J,parity)=0.
  240 continue
      do 250 type=0,numpar
        xspreeqtot(type)=0.
        xspreeqtotps(type)=0.
        xspreeqtotki(type)=0.
        xspreeqtotbu(type)=0.
  250 continue
      do 260 p=1,numparx
        do 260 type=0,numpar
          xssteptot(type,p)=0.
  260 continue
      if (flagang.or.flagddx) then
        do 270 iang=0,numangcont
          do 270 nen=0,numen
            do 270 type=0,numpar
              xscompad(type,nen,iang)=0.
              xsbinemisad(type,nen,iang)=0.
              xspreeqad(type,nen,iang)=0.
              xsmpreeqad(type,nen,iang)=0.
  270   continue
      endif
      do 300 i=0,numlev2
        do 300 type=0,numpar
          xspreeqdisc(type,i)=0.
  300 continue
      do 310 type=0,numpar
        xspreeqdisctot(type)=0.
  310 continue
      xspreeqdiscsum=0.
      if (flagmulpre) then
        do 320 n2=0,numparx
          do 320 n1=0,numparx
            do 320 nex=0,numex
              do 320 Nix=0,numNph
                do 320 Zix=0,numZph
                  xspopph(Zix,Nix,nex,n1,n2)=0.
  320   continue
        do 330 n4=0,numparx
          do 330 n3=0,numparx
            do 330 n2=0,numparx
              do 330 n1=0,numparx
                PP2(n1,n2,n3,n4)=0.
                Spre(n1,n2,n3,n4)=0.
                do 330 nex=0,numex
                  do 330 Nix=0,numNph
                    do 330 Zix=0,numZph
                      xspopph2(Zix,Nix,nex,n1,n2,n3,n4)=0.
  330   continue
      endif
      xspreeqsum=0.
      Esurf=0.
c
c Multi-step direct arrays
c
c preeqmode : designator for pre-equilibrium model
c numenmsd  : maximum number of energy points for DWBA calculation for
c             MSD
c numJmsd   : maximal spin for MSD
c xsdwin    : DWBA cross section as a function of incident energy,
c             outgoing energy and angular momentum
c xsdw      : DWBA angular distribution as a function of incident
c             energy, outgoing energy, angular momentum and angle
c Emsd      : MSD energy grid
c xscont    : continuum one-step direct cross section for MSD
c xscont1   : continuum one-step direct cross section for MSD
c             (unnormalized)
c xscontad  : continuum one-step direct angular distribution for MSD
c xscontad1 : continuum one-step direct angular distribution for MSD
c             (unnormalized)
c msdstep1  : continuum one-step direct cross section (unnormalized)
c msdstepad1: continuum one-step direct angular distribution
c             (unnormalized)
c msdstep0  : n-step cross section for MSD
c msdstep   : continuum n-step direct cross section
c msdstepad0: n-step angular distribution for MSD
c msdstepad : continuum n-step direct angular distribution
c
      if (preeqmode.eq.4) then
        do 410 i=0,2
          do 410 J=0,numJmsd
            do 410 nen2=0,numenmsd
              do 410 nen1=0,numenmsd
                xsdwin(nen1,nen2,J,i)=0.
  410   continue
        do 420 i=0,2
          do 420 iang=0,numangcont
            do 420 J=0,numJmsd
              do 420 nen2=0,numenmsd
                do 420 nen1=0,numenmsd
                  xsdw(nen1,nen2,J,iang,i)=0.
  420   continue
        do 430 nen2=0,numenmsd
          Emsd(nen2)=0.
          do 430 nen1=0,numenmsd
            do 430 J=0,numpar
              do 430 i=0,numpar
                xscont(i,J,nen1,nen2)=0.
                xscont1(i,J,nen1,nen2)=0.
                do 440 iang=0,numangcont
                  xscontad(i,J,nen1,nen2,iang)=0.
                  xscontad1(i,J,nen1,nen2,iang)=0.
  440           continue
  430   continue
        do 450 nen=0,numen
          do 450 i=0,numpar
            msdstep1(i,nen)=0.
            do 460 iang=0,numangcont
              msdstepad1(i,nen,iang)=0.
  460       continue
  450   continue
        do 470 nen=0,numen
          do 470 J=1,nummsd
            do 470 i=0,numpar
              msdstep(i,J,nen)=0.
              do 480 iang=0,numangcont
                msdstepad(i,J,nen,iang)=0.
  480         continue
  470   continue
        do 490 nen=0,numenmsd
          do 490 J=1,nummsd
            do 490 i=0,numpar
              msdstep0(i,J,nen)=0.
              do 500 iang=0,numangcont
                msdstepad0(i,J,nen,iang)=0.
  500         continue
  490   continue
      endif
c
c ********* Initialization of primary compound and binary arrays *******
c
c J2beg      : 2 * start of J summation
c J2end      : 2 * end of J summation
c xspopnuc   : population cross section per nucleus
c xspopdir   : direct population cross section per nucleus
c xspoppreeq : preequilibrium population cross section per nucleus
c xspopcomp  : compound population cross section per nucleus
c Fdir       : direct population fraction per nucleus
c Fpreeq     : preequilibrium population fraction per nucleus
c Fcomp      : compound population fraction per nucleus
c xspopex    : population cross section summed over spin and parity
c maxJ       : maximal J-value
c xspop      : population cross section
c xscompel   : compound elastic cross section
c xselastot  : total elastic cross section (shape + compound)
c xsnonel    : non-elastic cross section
c xscompnonel: total compound non-elastic cross section
c xsbinary   : cross section from initial compound to residual nucleus
c xsdisctot  : total cross section summed over discrete states
c xsdircont  : direct cross section for continuum
c xsdirect   : total direct cross section
c xsconttot  : total cross section for continuum
c xscompound : total compound cross section
c xscompcont : compound cross section for continuum
c lmaxhf     : maximal l-value for transmission coefficients
c contrib    : contribution to emission spectrum
c feedbinary : feeding from first compound nucleus
c flagspec   : flag for output of spectra
c binemis    : emission spectra from initial compound nucleus
c xscomp     : compound emission spectrum
c xsbinemis  : cross section for emission from first compound nucleus
c xsemis     : emission spectrum from compound nucleus
c xsngnspec  : total (projectile,gamma-ejectile) spectrum
c numlev     : maximum number of included discrete levels
c xsracappop : population cross section for radiative capture
c xspopex0   : binary population cross section for discrete states
c preeqpopex : pre-equilibrium population cross section summed over
c              spin and parity
c preeqpop   : pre-equilibrium population cross section
c xsdisc     : total cross section for discrete state
c xscompdisc : compound cross section for discrete state
c numang     : maximum number of angles
c compad     : compound angular distribution
c discad     : discrete state angular distribution
c cleg       : compound nucleus Legendre coefficient
c tleg       : total Legendre coefficient
c tlegnor    : total Legendre coefficient normalized to 1
c cleg0      : Legendre coefficient normalized to the first one
c transjl    : array for width fluctuation calculation
c
      J2beg=0
      J2end=0
      do 510 Nix=0,numN
        do 510 Zix=0,numZ
          xspopnuc(Zix,Nix)=0.
          xspopdir(Zix,Nix)=0.
          xspoppreeq(Zix,Nix)=0.
          xspopcomp(Zix,Nix)=0.
          Fdir(Zix,Nix)=0.
          Fpreeq(Zix,Nix)=0.
          Fcomp(Zix,Nix)=0.
  510 continue
      do 520 nex=0,numex
        xsracappopex(nex)=0.
        do 520 Nix=0,numN
          do 520 Zix=0,numZ
            xspopex(Zix,Nix,nex)=0.
            preeqpopex(Zix,Nix,nex)=0.
            maxJ(Zix,Nix,nex)=numJ
  520 continue
      do 530 parity=-1,1,2
        do 530 J=0,numJ
          do 530 nex=0,numex
            xsracappop(nex,J,parity)=0.
            do 530 Nix=0,numN
              do 530 Zix=0,numZ
                xspop(Zix,Nix,nex,J,parity)=0.
                preeqpop(Zix,Nix,nex,J,parity)=0.
  530 continue
      xscompel=0.
      xselastot=0.
      xsnonel=0.
      xscompnonel=0.
      xsbinary(-1)=0.
      do 540 type=0,numpar
        xsbinary(type)=0.
        xsdisctot(type)=0.
        xsdircont(type)=0.
        xsdirect(type)=0.
        xsconttot(type)=0.
        xscompound(type)=0.
        xscompcont(type)=0.
  540 continue
      do 550 nex=0,numex
        do 550 type=0,numpar
          lmaxhf(type,nex)=0.
          contrib(type,nex)=0.
          feedbinary(type,nex)=0.
  550 continue
      if (flagspec) then
        do 560 nen=0,numen
          do 560 nex=0,numex
            do 560 type=0,numpar
              binemis(type,nex,nen)=0.
  560   continue
        do 570 nen=0,numen
          do 570 type=0,numpar
            xscomp(type,nen)=0.
            xsbinemis(type,nen)=0.
            xsemis(type,nen)=0.
            xsngnspec(type,nen)=0.
  570   continue
      endif
      do 580 nex=0,numlev
        do 580 type=0,numpar
          xspopex0(type,nex)=0.
          xsdisc(type,nex)=0.
          xscompdisc(type,nex)=0.
  580 continue
      if (flagang.or.flagddx) then
        do 590 iang=0,numang
          do 590 nex=0,numlev
            do 590 type=0,numpar
              compad(type,nex,iang)=0.
              discad(type,nex,iang)=0.
  590   continue
      endif
      do 600 LL=0,3*numl
        do 600 nex=0,numlev
          do 600 type=0,numpar
            cleg(type,nex,LL)=0.
            tleg(type,nex,LL)=0.
            tlegnor(type,nex,LL)=0.
            cleg0(type,nex,LL)=0.
  600 continue
      do 610 l=1,numtrans
        do 610 i=0,5
          transjl(i,l)=0.
  610 continue
c
c ************* Initialization of multiple mission arrays **************
c
c xsngnsum    : sum over total (projectile,gamma-ejectile) cross
c               sections
c idnum       : counter for exclusive channel
c flagchannels: flag for exclusive channels calculation
c channelsum  : sum over exclusive channel cross sections
c numNchan    : maximal number of outgoing neutron units in individual
c               channel description
c numZchan    : maximal number of outgoing proton units in individual
c               channel description
c feedexcl    : feeding terms from compound excitation energy bin to
c               residual excitation energy bin
c popexcl     : population cross section of bin just before decay
c flagfission : flag for fission
c fisfeedex   : fission contribution from excitation energy bin
c tfis        : fission transmission coefficients
c xsfeed      : cross section from compound to residual nucleus
c xsngn       : total (projectile,gamma-ejectile) cross section
c xsgamdis    : discrete gamma-ray cross section
c xsgamdistot : total discrete gamma-ray cross section
c
      xsngnsum=0.
      idnum=-1
      if (flagchannels) then
        channelsum=0.
        do 710 nexout=0,numex+1
          do 710 nex=0,numex+1
            do 710 type=0,numpar
              do 710 Nix=0,numNchan
                do 710 Zix=0,numZchan
                  feedexcl(Zix,Nix,type,nex,nexout)=0.
  710   continue
        do 720 nex=0,numex
          do 720 Nix=0,numN
            do 720 Zix=0,numZ
              popexcl(Zix,Nix,nex)=0.
  720   continue
      endif
      if (flagfission) then
        do 730 nex=0,numex+1
          do 730 Nix=0,numN
            do 730 Zix=0,numZ
              fisfeedex(Zix,Nix,nex)=0.
              do 735 parity=-1,1,2
                do 735 J=0,numJ
                  fisfeedJP(Zix,Nix,nex,J,parity)=0.
  735         continue
  730   continue
      endif
      do 740 parity=-1,1,2
        do 740 J=0,numJ
          tfis(J,parity)=0.
  740 continue
      do 750 type=-1,6
        do 750 Nix=0,numN-2
          do 750 Zix=0,numZ-2
            xsfeed(Zix,Nix,type)=0.
  750 continue
      do 760 type=-1,6
        xsngn(type)=0.
  760 continue
      do 770 n2=0,numlev
        do 770 n1=0,numlev
          do 770 Nix=0,numN
            do 770 Zix=0,numZ
              xsgamdis(Zix,Nix,n1,n2)=0.
  770 continue
      do 780 Nix=0,numN
        do 780 Zix=0,numZ
          xsgamdistot(Zix,Nix)=0.
  780 continue
c
c ************* Initialization of channels and total arrays ************
c
c nin          : counter for incident energy
c idnumfull    : flag to designate maximum number of exclusive channels
c opennum      : total number of open channels
c numin,....   : maximal number of ejectile in channel description
c chanopen     : flag to open channel with first non-zero cross section
c numen2       : maximum number of outgoing energies
c xssumoutad   : angular distribution summed over mechanisms
c xsdiscoutad  : smoothed angular distribution for discrete state
c xscompoutad  : compound emission angular distribution
c xspreeqoutad : preequilibrium angular distribution per particle type
c xsmpreeqoutad: multiple preequilibrium angular distribution
c xsbranch     : branching ratio for isomeric cross section
c xsparcheck   : total particle production cross section
c xsspeccheck  : total particle production spectra
c
      if (nin.eq.1) then
        idnumfull=.false.
        opennum=-1
        do 810 in=0,numin
          do 810 ip=0,numip
            do 810 id=0,numid
              do 810 it=0,numit
                do 810 ih=0,numih
                  do 810 ia=0,numia
                    chanopen(in,ip,id,it,ih,ia)=.false.
  810   continue
      endif
      if (flagang.or.flagddx) then
        do 820 iang=0,numangcont
          do 840 nen=0,numen2
            do 840 type=0,numpar
              xssumoutad(type,nen,iang)=0.
              xsdiscoutad(type,nen,iang)=0.
              xscompoutad(type,nen,iang)=0.
              xspreeqoutad(type,nen,iang)=0.
              xsmpreeqoutad(type,nen,iang)=0.
  840     continue
  820   continue
      endif
      do 850 n1=0,numlev
        do 850 Nix=0,numN
          do 850 Zix=0,numZ
            xsbranch(Zix,Nix,n1)=0.
  850 continue
      do 860 type=0,numpar
        xsparcheck(type)=0.
        do 870 nen=0,numen
          xsspeccheck(type,nen)=0.
  870   continue
  860 continue
c
c ************* Initialization of total cross section arrays ***********
c
c nexmax       : maximum excitation energy bin for residual nucleus
c xsexclusive  : exclusive single channel cross section
c xsexclcont   : exclusive single channel cross section for continuum
c multiplicity : particle multiplicity
c xsparticle   : total particle production cross section
c xsfistot     : total fission cross section
c maxA         : maximal number of nucleons away from the initial
c                compound nucleus
c maxZ         : maximal number of protons away from the initial
c                compound nucleus
c maxN         : maximal number of protons away from the initial
c                compound nucleus
c xsresprod    : total residual production (= reaction) cross section
c xsmassprod   : residual production cross section per mass unit
c xstot6       : total cross section (neutrons only) for ENDF-6 file
c xsreac6      : reaction cross section for ENDF-6 file
c xsnon6       : non-elastic cross section for ENDF-6 file
c xsopt6       : optical model reaction cross section for ENDF-6 file
c xselassh6    : shape elastic cross section (neutrons only) for ENDF-6
c                file
c xselas6      : total elastic cross section (neutrons only) for ENDF-6
c                file
c xscompel6    : compound elastic cross section
c e6           : energies of ENDF-6 energy grid in MeV
c fxsbinary    : cross section from initial compound to residual nucleus
c fxsexclusive : exclusive single channel cross section
c fxsdisctot   : total cross section summed over discrete states
c fxsexclcont  : exclusive single channel cross section for continuum
c fxsngn       : total (projectile,gamma-ejectile) cross section
c fxsdisc      : total cross section for discrete state
c fxsdirdisc   : direct cross section for discrete state
c fxscompdisc  : compound cross section for discrete state
c fxspopnuc    : population cross section per nucleus
c fxspopex     : population cross section summed over spin and parity
c fxsbranch    : branching ratio for isomeric cross section
c idchannel    : identifier for exclusive channel
c fxschannel   : channel cross section
c fxsgamchannel: gamma channel cross section
c fxsratio     : ratio of exclusive cross section over residual
c                production cross section (for exclusive gamma ray
c                intensities)
c fxsgamdischan: discrete gamma channel cross section
c fxschaniso   : channel cross section per isomer
c fexclyield   : exclusive channel yield per isomer
c fxsgamdischan: discrete gamma channel cross section
c fxsnonel     : non-elastic cross section
c fxselastot   : total elastic cross section (shape + compound)
c fxstotinc    : total cross section (neutrons only) for incident
c                channel
c fxscompel    : compound elastic cross section
c fxselasinc   : total elastic cross section (neutrons only) for
c                incident channel
c fxsreacinc   : reaction cross section for incident channel
c fxscompnonel : total compound non-elastic cross section
c fxsdirdiscsum: total direct cross section
c fxspreeqsum  : total preequilibrium cross section summed over
c                particles
c fisstring    : string for exclusive fission reaction channel
c
      do 910 type=0,numpar
        nexmax(type)=-1
        xsexclusive(type)=0.
        xsexclcont(type)=0.
        multiplicity(type)=0.
        xsparticle(type)=0.
  910 continue
      xsfistot=0.
      maxA=maxZ+maxN
      xsresprod=0.
      do 920 ia=0,numA
        xsmassprod(ia)=0.
  920 continue
      if (nin.eq.1) then
        do 930 nen=1,numen6
          xstot6(nen)=0.
          xsreac6(nen)=0.
          xsnon6(nen)=0.
          xsnonel6(nen)=0.
          xsopt6(nen)=0.
          xselassh6(nen)=0.
          xselas6(nen)=0.
          xscompel6(nen)=0.
          e6(nen)=0.
  930   continue
        do 940 nen=1,numenlow
          do 950 type=0,numpar
            fxsbinary(nen,type)=0.
            fxsexclusive(nen,type)=0.
            fxsdisctot(nen,type)=0.
            fxsexclcont(nen,type)=0.
            fxsngn(nen,type)=0.
            do 960 n1=0,numlev
              fxsdisc(nen,type,n1)=0.
              fxsdirdisc(nen,type,n1)=0.
              fxscompdisc(nen,type,n1)=0.
  960       continue
  950     continue
          fxsngn(nen,-1)=0.
          fxsreacinc(nen)=0.
          fxselastot(nen)=0.
          fxsnonel(nen)=0.
          fxscompel(nen)=0.
          fxselasinc(nen)=0.
          fxstotinc(nen)=0.
          fxscompnonel(nen)=0.
          fxsdirdiscsum(nen)=0.
          fxspreeqsum(nen)=0.
          do 970 Nix=0,numN
            do 970 Zix=0,numZ
              fxspopnuc(nen,Zix,Nix)=0.
              do 980 n1=0,numlev
                fxspopex(nen,Zix,Nix,n1)=0.
                fxsbranch(nen,Zix,Nix,n1)=0.
  980         continue
  970     continue
          do 990 i=0,numchantot
            idchannel(i)=-1
            fxschannel(nen,i)=0.
            fxsgamchannel(nen,i)=0.
            fxsratio(nen,i)=0.
  990     continue
  940   continue
      endif
      do 1010 n2=0,numlev
        do 1020 n1=0,numlev
          do 1030 i=0,numchantot
            do 1040 nen=1,numenlow
              fxsgamdischan(nen,i,n1,n2)=0.
 1040       continue
 1030     continue
 1020   continue
 1010 continue
      do 1050 n1=0,numlev
        do 1060 i=0,numchantot
          do 1070 nen=1,numenlow
            fxschaniso(nen,i,n1)=0.
            fexclyield(nen,i,n1)=0.
 1070     continue
 1060   continue
 1050 continue
      do 1080 i=0,numchantot
        fisstring(i)='                  '
 1080 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
