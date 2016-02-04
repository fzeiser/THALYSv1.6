      subroutine comptarget
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : October 23, 2013
c | Task  : Compound reaction for initial compound nucleus
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          elastic
      integer          parspin2i,pspin2i,Zcomp,Ncomp,updown,l,parity,J2,
     +                 J,jj2beg,jj2end,l2beg,l2end,jj2primebeg,
     +                 jj2primeend,l2primebeg,l2primeend,jj2,l2,ihill,
     +                 type,parspin2o,pspin2o,Zix,Nix,NL,nexout,l2maxhf,
     +                 Pprimebeg,Pprimeend,Irspin2beg,Irspin2end,J2res,
     +                 Pprime,pardif2,Irspin2,Ir,jj2prime,l2prime,
     +                 lprime,modl,irad,updown2,ielas,LLmax,iphase,LL,
     +                 nex
      real             Tinc,fluxsum,ratio,Tout,Rspin,angfac,rJ,rlin,
     +                 rjin,rlout,rjout,phase,rLL,clin,clebsch,ra1in,
     +                 racah,ra2in,clout,ra1out,ra2out,Ablatt
      double precision factor1,fisterm,sumIPE,sumIP,sumIPas,rho,sumjl,
     +                 compterm
c
c *** Check if initial compound nucleus calculation needs to be done ***
c
c flaginitpop: flag for initial population distribution
c xsflux     : cross section flux
c xseps      : limit for cross sections
c
c This may occur when direct + pre-equilibrium reactions completely
c exhaust the reaction cross section for the binary channel, and
c when the initial system is a excitation energy population.
c
c
      if (flaginitpop) return
      if (xsflux.le.xseps) return
c
c *************************** Initialization ***************************
c
c Wab          : width fluctuation factor
c parspin2i    : 2 * particle spin for incident channel
c parspin      : spin of particle
c k0           : index of incident particle
c pspin2i,spin2: 2 * spin of particle (usually) for incident channel
c
c Initially, Wab is set to 1 (no width fluctuations).
c
      Wab=1.
      parspin2i=int(2.*parspin(k0))
      pspin2i=spin2(k0)
c
c The level densities and transmission coefficients can be prepared
c before the nested loops over all quantum numbers.
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c dExinc     : excitation energy bin for mother nucleus
c densprepare: subroutine to prepare energy grid, level density and
c              transmission coefficient information for compound nucleus
c
      Zcomp=0
      Ncomp=0
      dExinc=0.
      call densprepare(Zcomp,Ncomp,1)
c
c *** Output of flux conservation check of transmission coefficients ***
c
c flagcheck: flag for output of numerical checks
c flagwidth: flag for width fluctuation calculation
c wmode    : designator for width fluctuation model
c
      if (flagcheck.and.flagwidth) then
        write(*,'(/" ++++++++++ CHECK OF FLUX CONSERVATION",
     +    " OF TRANSMISSION COEFFICIENTS ++++++++++")')
        if (wmode.eq.0) write(*,'(" Hauser-Feshbach model"/)')
        if (wmode.eq.1) write(*,'(" Moldauer model")')
        if (wmode.eq.2) write(*,'(" HRTW model"/)')
        if (wmode.eq.3) write(*,'(" GOE model"/)')
      endif
c
c ************** Initialization of transmission coefficients ***********
c
c The transmission coefficients Tjlnex have values interpolated from
c the Tjl coefficients of the emission energy grid. For the incident
c channel, the Tjl are exactly calculated. Hence, we use exact values
c for the incident channel in the transmission coefficient array.
c
c updown     : spin index for transmission coefficient
c lmaxinc    : maximal l-value for transmission coefficients for
c              incident channel
c Tjlnex     : transmission coefficients as a function of
c              particle type, energy, spin and l-value
c Ltarget    : excited level of target
c Tjlinc,Tinc: transmission coefficients as a function of j and l
c              for the incident channel
c
      do 10 updown=-1,1
        do 10 l=0,lmaxinc
          Tjlnex(k0,Ltarget,updown,l)=Tjlinc(updown,l)
   10 continue
      nex=maxex(Zcomp,Ncomp)
c
c Initialisation of total astrophysical transmission coefficient
c
c flagastro  : flag for calculation of astrophysics reaction rate
c Tastrotot  : total transmission coefficient for astrophysical case
c rhoastrotot: total level density for astrophysical case
c
      if (flagastro) then
        Tastrotot=0.
        rhoastrotot=0.
      endif
c
c ************** Loop over incoming and outgoing channels **************
c
c 110: Sum over compound nucleus parity
c
c parity : parity
c pardif : difference between target and compound nucleus parity
c targetP: parity of target
c
c The variable pardif is used as an indicator of parity conservation
c for the incident channel.
c
      do 110 parity=-1,1,2
        pardif=abs(targetP-parity)/2
c
c 120: Sum over total angular momentum J (J2) of compound nucleus
c
c J2          : 2 * J
c flagfission : flag for fission
c nfisbar     : number of fission barrier parameters
c J2beg       : 2 * start of J summation
c J2end       : 2 * end of J summation
c tfission    : subroutine for fission transmission coefficients
c compprepare : subroutine to prepare information for initial
c               compound nucleus
c denomhf     : denominator for compound nucleus formula
c widthprepare: subroutine for preparation of width fluctuation
c               corrections
c tnumi       : counter for width fluctuation calculation
c
c There are two possible types of calculation for the initial compound
c nucleus. If either width fluctuation corrections or compound nucleus
c angular distributions are wanted, we need to sum explicitly over all
c possible quantum numbers before we calculate the width fluctuation
c or angular factor. If not, the sum over j,l of the transmission
c coefficients can be lumped into one factor, which decreases the
c calculation time. In the latter case, the partial decay widths
c enumhf are calculated in subroutine compprepare.
c
c In order to get do-loops running over integer values, certain quantum
c numbers are multiplied by 2, which can be seen from a 2 present in the
c corresponding variable names. For each loop, the begin and end point
c is determined from the triangular rule.
c
c For every J and P (parity), first the denominator (total width)
c denomhf for the Hauser-Feshbach formula is constructed in subroutine
c compprepare. Also width fluctuation variables that only depend on J
c and P and not on the other angular momentum quantum numbers are
c calculated in subroutine widthprepare. For the width fluctuation
c calculation, all transmission coefficients need to be placed in one
c sequential array. Therefore, a counter tnum needs to be followed
c to keep track of the proper index for the transmission coefficients.
c
        do 120 J2=J2beg,J2end,2
          J=J2/2
          if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0)
     +      call tfission(Zcomp,Ncomp,nex,J2,parity)
          call compprepare(Zcomp,Ncomp,J2,parity)
          if (denomhf.eq.0.) goto 120
          if (flagwidth) call widthprepare
          tnumi=0
c
c Initially, we assume that no width fluctuation corrections or
c angular distributions are calculated. This means the various loops
c over l and j for the incident and outgoing channel do not need to
c be performed. In this case, the terms needed for the Hauser-Feshbach
c calculation only depend on J and P, and the
c corresponding widths enumhf and denomhf have already been calculated
c in subroutine compprepare, where the loops over l and j were done.
c
c jj2beg     : 2 * start of j summation
c jj2end     : 2 * end of j summation
c l2beg      : 2 * start of l summation
c l2end      : 2 * end of l summation
c jj2primebeg: 2 * start of j' summation
c jj2primeend: 2 * end of j' summation
c l2primebeg : 2 * start of l summation
c l2primeend : 2 * end of l summation
c
          jj2beg=1
          jj2end=1
          l2beg=1
          l2end=1
          jj2primebeg=1
          jj2primeend=1
          l2primebeg=1
          l2primeend=1
c
c 130: Sum over j (jj2) of incident channel
c
c flagcompang  : flag for compound angular distribution calculation
c flagurr      : flag for output of unresolved resonance parameters
c targetspin2  : 2 * spin of target
c jj2          : 2 * j
c
c On-set of loop over j in the case of width fluctuations or
c angular distributions.
c
          if (flagwidth.or.flagcompang.or.flagurr) then
            jj2beg=abs(J2-targetspin2)
            jj2end=J2+targetspin2
          endif
          do 130 jj2=jj2beg,jj2end,2
c
c 140: Sum over l of incident channel
c
c l2: 2 * l
c
c On-set of loop over l in the case of width fluctuations or
c angular distributions.
c
            if (flagwidth.or.flagcompang.or.flagurr) then
              l2beg=abs(jj2-parspin2i)
              l2end=jj2+parspin2i
              l2end=min(l2end,2*lmaxinc)
            endif
            do 140 l2=l2beg,l2end,2
              l=l2/2
c
c Check parity conservation and make index for transmission coefficient
c for width fluctuation calculation.
c
c
c If the parity of the target nucleus is equal (unequal) to the parity
c of compound nucleus, i.e. pardif=0(1), the l-value must be even (odd).
c
              if (flagwidth.or.flagcompang.or.flagurr) then
                if (mod(l,2).ne.pardif) goto 140
                updown=(jj2-l2)/pspin2i
                Tinc=Tjlinc(updown,l)
                tnumi=tnumi+1
              endif
c
c 160: Sum over outgoing channels
c
c fluxsum      : check for conservation of flux per P,J,j,l
c tnumo,tnuminc: counters for width fluctuation calculation
c tnum         : total number of transmission coefficients
c numhill      : maximum number of Hill-Wheeler points
c ihill        : counter for Hill-Wheeler magnitude
c ratio,factor1: help variables
c tfisA        : transmission coefficient for Hill-Wheeler magnitude
c widthfluc    : subroutine for width fluctuation correction
c tfis         : fission transmission coefficients
c xsbinary     : cross section from initial compound to residual nucleus
c CNfactor     : factor for compound nucleus cross section:
c                pi/[ k**2 (2s+1)(2I+1) ]
c fisfeedJP    : fission contribution from excitation energy bin per J,P
c
              fluxsum=0.
c
c 1. Fission channel
c
c The fission contribution is calculated and added to the binary fission
c cross section. Also the transmission coefficient index for width
c fluctuations is increased.
c
              if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0) then
                if (flagwidth.or.flagcompang.or.flagurr) then
                  tnumo=tnum
                  if (flagwidth.and.wmode.ge.1) then
                    tnumo=tnum-numhill
                    factor1=0.
                    if (tfisA(J,parity,0).gt.0.) then
                      do 150 ihill=1,numhill
                        tnumo=tnumo+1
                        ratio=tfisA(J,parity,ihill)/tfisA(J,parity,0)
                        if (ratio.eq.0.) goto 150
                        call widthfluc(0)
                        factor1=factor1+
     +                    real(Tinc*tfis(J,parity)/denomhf*Wab)*ratio
  150                 continue
                      fluxsum=fluxsum+factor1
                    endif
                  else
                    factor1=real(Tinc*tfis(J,parity)/denomhf)
                  endif
                else
                  factor1=real(feed*tfis(J,parity)/denomhf)
                endif
                fisterm=CNfactor*(J2+1.)*factor1
                xsbinary(-1)=xsbinary(-1)+fisterm
                fisfeedJP(0,0,maxex(0,0)+1,J,parity)=fisterm
c
c Extract (L,J) dependent parameters for URR (Gilles Noguere)
c
c lurr      : maximal orbital angular momentum for URR
c Turr      : (l,j) transmission coefficient for URR calculation
c xsbinarylj: (l,j) cross section for URR calculation
c nulj      : (l,j) number of degrees of freedom for URR calculation
c
                if (flagurr.and.l.le.lurr) then
                  Turrlj(-1,l,J)=factor1*denomhf/(Tinc*Wab)
                  xsbinarylj(-1,l,J)=xsbinarylj(-1,l,J)+
     +              CNfactor*(J2+1.)*real(factor1)
                  nulj(-1,l,J)=1
                endif
              endif
c
c 2. Gamma and particle channels
c
c parskip   : logical to skip outgoing particle
c parspin2o : 2 * particle spin for outgoing channel
c pspin2o   : 2 * spin of particle (usually) for outgoing channel
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c Nlast,NL  : last discrete level
c
              tnumo=tnum+1
              do 160 type=0,6
                if (type.eq.1) tnumo=tnuminc
                if (parskip(type)) goto 160
                parspin2o=int(2.*parspin(type))
                pspin2o=spin2(type)
                Zix=Zindex(Zcomp,Ncomp,type)
                Nix=Nindex(Zcomp,Ncomp,type)
                NL=Nlast(Zix,Nix,0)
c
c 170: Sum over outgoing excitation energies
c
c sumIPE : compound contribution summed over residual spin and parity
c          and energy
c maxex  : maximum excitation energy bin for residual nucleus
c l2maxhf: 2 * lmaxhf
c lmaxhf : maximal l-value for transmission coefficients
c elastic: designator for elastic channel
c
c This loop is over all discrete levels and continuum bins of the
c final nucleus.
c
                sumIPE=0.
                do 170 nexout=0,maxex(Zix,Nix)
                  l2maxhf=2*lmaxhf(type,nexout)
                  elastic=(type.eq.k0.and.nexout.eq.Ltarget)
c
c Initialization of summations
c
c Pprimebeg : start of residual parity summation
c parlev    : parity of level
c Pprimeend : end of residual parity summation
c Irspin2beg: 2 * start of residual spin summation
c jdis      : spin of level
c Irspin2end: 2 * end of residual spin summation
c J2res     : help variable
c maxJ      : maximal J-value
c sumIP     : compound contribution summed over residual spin and parity
c sumIPas   : sumIP for astrophysics
c
c For discrete states, the begin and end points of the residual
c spin/parity summation are both set equal to the residual discrete
c level spin/parity.
c
                  if (nexout.le.NL) then
                    Pprimebeg=parlev(Zix,Nix,nexout)
                    Pprimeend=Pprimebeg
                    Irspin2beg=int(2.*jdis(Zix,Nix,nexout))
                    Irspin2end=Irspin2beg
                  else
c
c For the continuum, the begin and end points of the residual
c spin/parity summation are set to the maximally accessible values.
c
                    Pprimebeg=-1
                    Pprimeend=1
                    J2res=J2+parspin2o
                    Irspin2beg=mod(J2res,2)
                    Irspin2end=J2res+l2maxhf
                    Irspin2end=min(Irspin2end,2*maxJ(Zix,Nix,nexout))
                  endif
                  sumIP=0.
                  if (flagastro) sumIPas=0.
c
c 180: Sum over residual parity
c
c pardif2: difference between residual and compound nucleus parity
c
c The variable pardif2 is used as an indicator of parity conservation
c for the outgoing channel.
c
                  do 180 Pprime=Pprimebeg,Pprimeend,2
                    pardif2=abs(parity-Pprime)/2
c
c 190: Sum over residual spin
c
c Irspin2 : 2 * residual spin
c Ir      : residual spin
c rho,rho0: integrated level density
c
                    do 190 Irspin2=Irspin2beg,Irspin2end,2
                      Ir=Irspin2/2
                      rho=rho0(type,nexout,Ir,Pprime)
                      if (rho.lt.1.e-20) goto 190
c
c 200: Sum over j (jj2) of outgoing channel
c
c sumjl   : compound contribution summed over residual j and l
c jj2prime: 2 * j'
c
c On-set of loop over jprime in the case of width fluctuations or
c angular distributions.
c
                      sumjl=0.
                      if (flagwidth.or.flagcompang.or.flagurr) then
                        jj2primebeg=abs(J2-Irspin2)
                        jj2primeend=J2+Irspin2
                      endif
                      do 200 jj2prime=jj2primebeg,jj2primeend,2
c
c 210: Sum over l of outgoing channel
c
c l2prime: 2 * l'
c modl   : help variable
c
c On-set of loop over lprime in the case of width fluctuations or
c angular distributions.
c
                        if (flagwidth.or.flagcompang.or.flagurr) then
                          l2primebeg=abs(jj2prime-parspin2o)
                          l2primeend=jj2prime+parspin2o
                          l2primeend=min(l2primeend,l2maxhf)
                          if (type.eq.0) l2primebeg=max(l2primebeg,2)
                        endif
                        do 210 l2prime=l2primebeg,l2primeend,2
                          lprime=l2prime/2
                          modl=mod(lprime,2)
c
c We include photons as a special case, with the multipole radiation
c selection rules (irad=0: M-transition, irad=1: E-transition)
c
                          if (flagwidth.or.flagcompang.or.flagurr) then
c
c 1. Photons
c
c irad: variable to indicate M(=0) or E(=1) radiation
c Tout: transmission coefficient
c Tgam: gamma transmission coefficient
c
                            if (type.eq.0) then
                              if (pardif2.eq.modl) then
                                irad=1
                              else
                                irad=0
                              endif
                              Tout=Tgam(nexout,lprime,irad)
                            else
c
c 2. Particles
c
c If the parity of the residual nucleus is equal (unequal) to the parity
c of compound nucleus,i.e. pardif2=0(1), the l-value must be even (odd).
c
c updown2: spin index for transmission coefficient
c
                              if (modl.ne.pardif2) goto 210
                              updown2=(jj2prime-l2prime)/pspin2o
                              Tout=Tjlnex(type,nexout,updown2,lprime)
                            endif
c
c ** Populate the outgoing channels using the compound nucleus formula *
c
c ielas: designator for elastic channel
c
c We determine the index for the width fluctuation calculation and call
c the subroutine that calculates the correction factor. The contribution
c of this particular outgoing channel is added to the sum for the
c incident channel (to check for flux conservation).
c
                            if (type.ge.1) tnumo=tnumo+1
                            if (flagwidth.and.wmode.ge.1) then
                              ielas=0
                              if (elastic.and.jj2.eq.jj2prime.and.
     +                          l2.eq.l2prime) ielas=1
                              call widthfluc(ielas)
                            endif
                            factor1=real(Tinc*rho*Tout/denomhf*Wab)
                            fluxsum=fluxsum+factor1
                          else
c
c If NO width fluctuation corrections or angular distributions are
c required, the partial and total decay widths were already determined
c in subroutine compprepare. This means we are now in the short loop
c with l=j=lprime=jprime=1.
c
c enumhf     : enumerator for compound nucleus formula
c flagastrogs: flag for calculation of astrophysics reaction rate with
c              target in ground state only
c
                            factor1=feed*
     +                        enumhf(type,nexout,Ir,Pprime)/denomhf
                            if (flagastro.and..not.flagastrogs.and.
     +                        type.eq.k0) then
                              sumIPas=sumIPas+CNfactor*(J2+1)*
     +                          enumhf(type,nexout,Ir,Pprime)**2/denomhf
                            endif
                          endif
c
c The contribution is added to the total width for the particular
c incident l and j.
c
c compterm: partial contribution to compound nucleus term
c
                          compterm=CNfactor*(J2+1.)*factor1
                          sumjl=sumjl+compterm
c
c ******************* Compound angular distributions *******************
c
c Rspin       : residual spin
c angfac      : help variable
c fourpi      : 4.*pi
c LLmax       : maximal Legendre order
c rJ,rLL      : help variables
c rlin,rjin   : help variables
c rlout,rjout : help variables
c iphase      : help variable
c targetspin  : spin of target
c phase       : phase
c sgn         : +1 for even argument, -1 for odd argument
c clin,clout  : Clebsch-Gordan coefficients
c clebsch     : function for Clebsch-Gordan coefficients
c logfact     : factorial logarithm
c numfact     : number of terms for factorial logarithm
c ra1in,ra1out: Racah coefficients
c ra2in,ra2out: Racah coefficients
c racah       : function for racah coefficients
c Ablatt      : Blatt-Biedenharn A-factor
c cleg        : compound nucleus Legendre coefficient
c
c Compound angular distributions are calculated for discrete states
c only. It can be easily generalized to the continuum, but we postpone
c that until we find a reason to do so (the deviation from isotropy is
c assumed to be negligible). Note that we are still inside the most
c inner loop 210, i.e. as indicated in the manual, all quantum numbers
c are required for a proper calculation of the compound angular
c distribution. The Legendre coefficients are divide by (2LL+1) to
c bring them on the same level as the direct reaction Legendre
c coefficients that come out of ECIS.
c
                          if (flagcompang.and.nexout.le.NL) then
                            Rspin=0.5*Irspin2
                            angfac=(J2+1)*(jj2prime+1)*(jj2+1)*(l2+1)*
     +                        (l2prime+1)/fourpi
                            LLmax=min(l2,l2prime)
                            LLmax=min(LLmax,J2)
                            rJ=0.5*J2
                            rlin=real(l)
                            rjin=0.5*jj2
                            rlout=real(lprime)
                            rjout=0.5*jj2prime
                            iphase=int(abs(Rspin-parspin(type)-
     +                        targetspin+parspin(k0))+0.1)
                            phase=sgn(iphase)
                            do 220 LL=0,LLmax,2
                              rLL=real(LL)
                              clin=clebsch(rlin,rlin,rLL,0.,0.,0.,
     +                          logfact,numfact)
                              ra1in=racah(rJ,rjin,rJ,rjin,targetspin,
     +                          rLL,logfact,numfact)
                              ra2in=racah(rjin,rjin,rlin,rlin,rLL,
     +                          parspin(k0),logfact,numfact)
                              clout=clebsch(rlout,rlout,rLL,0.,0.,0.,
     +                          logfact,numfact)
                              ra1out=racah(rJ,rjout,rJ,rjout,Rspin,rLL,
     +                          logfact,numfact)
                              ra2out=racah(rjout,rjout,rlout,rlout,rLL,
     +                          parspin(type),logfact,numfact)
                              Ablatt=phase*angfac*clin*ra1in*ra2in*
     +                          clout*ra1out*ra2out
                              cleg(type,nexout,LL)=cleg(type,nexout,LL)+
     +                          Ablatt*compterm/(2*LL+1)
  220                       continue
                          endif
c
c ************************* End of all loops ***************************
c
c xspop     : population cross section
c xspopex   : population cross section summed over spin and parity
c xscompcont: compound cross section for continuum
c contrib   : contribution to emission spectrum
c xspopnuc  : population cross section per nucleus
c
c Increment of the the population arrays of the residual nuclei, the
c binary reaction cross sections and the contrib array, which will be
c used for the interpolation of the compound emission spectra.
c
  210                   continue
  200                 continue
                      xspop(Zix,Nix,nexout,Ir,Pprime)=
     +                  xspop(Zix,Nix,nexout,Ir,Pprime)+sumjl
                      sumIP=sumIP+sumjl
  190               continue
  180             continue
                  xspopex(Zix,Nix,nexout)=xspopex(Zix,Nix,nexout)+sumIP
                  if (nexout.gt.NL) then
                    xscompcont(type)=xscompcont(type)+sumIP
                    contrib(type,nexout)=contrib(type,nexout)+sumIP
                  endif
c
c Compound elastic scattering is excluded from the residual production
c cross sections.
c
                  if (.not.elastic) sumIPE=sumIPE+sumIP
                  if (flagastro.and..not.flagwidth.and.type.eq.k0.and.
     +              sumIP-sumIPas.gt.xseps)  sumIPE=sumIPE-sumIPas
  170           continue
                xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)+sumIPE
                xsbinary(type)=xsbinary(type)+sumIPE
c
c Extract (L,J) dependent parameters for URR (advice of Gilles Noguere)
c
                if (flagurr.and.l.le.lurr) then
                  Turrlj(type,l,J)=
     +              sumIPE*denomhf/(CNfactor*(J2+1.)*Tinc*Wab)
                  xsbinarylj(type,l,J)=xsbinarylj(type,l,J)+
     +              real(sumIPE)
                  if (type.eq.0) then
                    nulj(type,l,J)=nulj(type,l,J)+1
                  else
                    nulj(type,l,J)=1
                  endif
                endif
  160         continue
c
c
c Determine angular momentum range for URR calculation
c
c Turrljinc   : incident channel (l,j) transmission coefficient for
c               URR calculation
c Purrlj      : (l,j) parity for URR calculation
c lminU,lmaxU : minimal and maximal orbital angular momentum
c JminU,JmaxU : minimal and maximal total angular momentum
c
              if (flagurr.and.l.le.lurr) then
                Turrljinc(l,J)=Turrljinc(l,J)+Tinc
                Purrlj(l,J)=parity
                if (l.le.lminU) lminU=l
                if (l.ge.lmaxU) lmaxU=l
                if (J.le.JminU(l)) JminU(l)=J
                if (J.ge.JmaxU(l)) JmaxU(l)=J
              endif
c
c ****** Check of flux conservation of transmission coefficients *******
c
c cparity: parity (character)
c
c This check is included to test the stability of the width fluctuation
c calculation.
c
              if (flagcheck.and.flagwidth) then
                if (fluxsum.eq.0.) fluxsum=Tinc
                write(*,'(" Parity=",a1,"  J=",f4.1,"  j=",f4.1,
     +            "  l=",i2,"  T(j,l)=",1p,e12.5,
     +            "  Sum over outgoing channels=",1p,e12.5,0p,
     +            "  Ratio=",f8.5)') cparity(parity),0.5*J2,0.5*jj2,l,
     +            Tinc,fluxsum,Tinc/fluxsum
              endif
  140       continue
  130     continue
  120   continue
  110 continue
c
c ************************** Astrophysics ******************************
c
c transeps: absolute limit for transmission coefficient
c parZ    : charge number of particle
c k0      : index of incident particle
c parN    : neutron number of particle
c
      if (flagastro.and..not.flagastrogs) then
        if (Tastrotot.ge.transeps.and.rhoastrotot.ge.0..and.flagwidth)
     +    call astrotarget
        if (flagwidth) then
          if (maxex(parZ(k0),parN(k0))+1.gt.3) ewfc=Einc
        endif
      endif
c
c *********** Output of fission transmission coefficients **************
c
c flagfisout : flag for output of fission information
c tfissionout: subroutine for output of fission transmission
c              coefficients
c
      if (flagfisout) call tfissionout(Zcomp,Ncomp,nex)
c
c **** ECIS calculation of compound cross sections (reference only) ****
c
c In addition to a calculation by TALYS, a compound nucleus run by ECIS
c can be requested. The results will however not be used for TALYS but
c are just for comparison.
c
c flageciscomp: flag for compound nucleus calculation by ECIS
c raynalcomp  : subroutine for ECIS calculation of compound cross
c               sections (reference only)
c
      if (flageciscomp) call raynalcomp
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
