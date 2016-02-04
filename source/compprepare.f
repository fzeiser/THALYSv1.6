      subroutine compprepare(Zcomp,Ncomp,J2,parity)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : August 31, 2007
c | Task  : Prepare information for initial compound nucleus
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,J2,parity,J,parspin2i,pspin2i,jj2beg,
     +                 jj2end,jj2,l2beg,l2end,l2,l,updown,i,parspin2o,
     +                 pspin2o,type,Zix,Nix,nexout,l2maxhf,Pprimebeg,
     +                 Pprimeend,Irspin2beg,Irspin2end,J2res,Pprime,
     +                 pardif2,Irspin2,Ir,J2maxastro,jj2primebeg,
     +                 jj2primeend,jj2prime,l2primebeg,l2primeend,
     +                 l2prime,lprime,modl,irad,updown2,ihill
      real             Tinc,Tout
      double precision gamwidth,rho,factor,ratio,tfishill
c
c The transmission coefficients are put in arrays for possible width
c fluctuation calculations. Also the total width denomhf appearing in
c the denominator of the compound nucleus formula is created. Note that
c the complete subroutine is performed inside the loop over J and P in
c subroutine comptarget.
c
c *************************** Initialization ***************************
c
c J2           : 2 * J
c denomhf      : denominator for compound nucleus formula
c tnum         : counter for width fluctuation calculation
c feed         : feeding term for compound nucleus
c parspin2i    : 2 * particle spin for incident channel
c parspin      : spin of particle
c k0           : index of incident particle
c pspin2i,spin2: 2 * spin of particle (usually) for incident channel
c
      J=J2/2
      denomhf=0.
      tnum=0
      feed=0.
      parspin2i=int(2.*parspin(k0))
      pspin2i=spin2(k0)
c
c ********************* Loop over incident channels ********************
c
c In order to get do-loops running over integer values, certain quantum
c numbers are multiplied by 2, which can be seen from a 2 present in the
c corresponding variable names. For each loop, the begin and end point
c is determined from the triangular rule.
c
c We have to distinguish between incident particles and incident
c photons.
c
c A. Incident particles
c
      if (k0.gt.0) then
c
c 10: Sum over j (jj2) of incident channel
c
c jj2beg     : 2 * start of j summation
c targetspin2: 2 * spin of target
c jj2end     : 2 * end of j summation
c jj2        : 2 * j
c
        jj2beg=abs(J2-targetspin2)
        jj2end=J2+targetspin2
        do 10 jj2=jj2beg,jj2end,2
c
c 20: Sum over l of incident channel
c
c l2beg   : 2 * start of l summation
c l2end   : 2 * end of l summation
c lmaxinc : maximal l-value for transmission coefficients for
c           incident channel
c
          l2beg=abs(jj2-parspin2i)
          l2end=jj2+parspin2i
          l2end=min(l2end,2*lmaxinc)
          do 20 l2=l2beg,l2end,2
            l=l2/2
c
c Check parity conservation and make index for transmission
c coefficient.
c
c pardif     : difference between target and compound nucleus parity
c updown     : spin index for transmission coefficient
c Tjlinc,Tinc: transmission coefficients as a function of j and l
c              for the incident channel
c
c If the parity of the target nucleus is equal (unequal) to the parity
c of compound nucleus, i.e. pardif=0(1), the l-value must be even (odd).
c
            if (mod(l,2).ne.pardif) goto 20
            updown=(jj2-l2)/pspin2i
            Tinc=Tjlinc(updown,l)
c
c Information needed for width fluctuation calculation.
c
c flagwidth  : flag for width fluctuation calculation
c transjl    : array for width fluctuation calculation
c wpower     : power used for rho*(t**wpower)
c flagcompang: flag for compound angular distribution calculation
c
c For the width fluctuation calculation, all transmission coefficients
c need to be placed in one sequential array. Therefore, a counter tnum
c needs to be followed to keep track of the proper index for the
c transmission coefficients. The order inside the transjl array is:
c
c 1. Incident channel
c 2. Outgoing particle channels
c 3. Gamma channel
c 4. Fission channel (if present)
c
            if (flagwidth) then
              tnum=tnum+1
              transjl(0,tnum)=1.
              do 30 i=1,wpower
                transjl(i,tnum)=0.
                if (Tinc.gt.1.e-30**(1./i)) transjl(i,tnum)=Tinc**i
   30         continue
            else
              if (.not.flagcompang) feed=feed+Tinc
            endif
   20     continue
   10   continue
      else
c
c B. Incident photons
c
c gammax    : number of l-values for gamma multipolarity
c targetspin: spin of target
c irad      : variable to indicate M(=0) or E(=1) radiation
c targetP   : parity of target
c
c Note that for photons, the first index of Tjlinc represents
c radiation type (M or E)
c
        Tinc=0.
        if (J2.ne.targetspin2.or.J2.ne.0) then
          do 40 l=1,gammax
            if (0.5*J2.lt.targetspin-l.or.0.5*J2.gt.targetspin+l)
     +        goto 40
            irad=1
            if (parity.eq.targetP.and.mod(l,2).eq.1) irad=0
            if (parity.ne.targetP.and.mod(l,2).eq.0) irad=0
            Tinc=Tinc+Tjlinc(irad,l)
   40     continue
          feed=feed+Tinc
        endif
      endif
c
c There are two possible types of calculation for the initial compound
c nucleus. If either width fluctuation corrections or compound nucleus
c angular distributions are wanted, we need to sum explicitly over all
c possible quantum numbers before we calculate the width fluctuation
c or angular factor. If not, the sum over j,l of the transmission
c coefficients can be lumped into one factor, which decreases the
c calculation time. In the latter case, the partial decay widths are
c stored in enumhf.
c
c enumhf    : enumerator for compound nucleus formula
c Ltarget   : excited level of target
c tnuminc   : counter for width fluctuation calculation
c
      enumhf(k0,Ltarget,int(targetspin),targetP)=feed
      tnuminc=tnum
c
c ********************* Loop over outgoing channels ********************
c
c 110: Sum over outgoing particles (skip incident channel)
c
c 1. Fission
c
c flagfission: flag for fission
c nfisbar    : number of fission barrier parameters
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c fiswidth   : fission width
c tfis       : fission transmission coefficients
c
c The fission contribution is calculated and added to the total decay
c width.
c
      if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0) then
        fiswidth=tfis(J,parity)
        denomhf=denomhf+fiswidth
      endif
c
c 2. Gamma and particle channels
c
c gamwidth : sum over all gamma transmission coefficients
c parskip  : logical to skip outgoing particle
c parspin2o: 2 * particle spin for outgoing channel
c pspin2o  : 2 * spin of particle (usually) for incident channel
c Zix      : charge number index for residual nucleus
c Nix      : neutron number index for residual nucleus
c
      gamwidth=0.
      do 110 type=0,6
        if (parskip(type)) goto 110
        parspin2o=int(2.*parspin(type))
        pspin2o=spin2(type)
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
c
c 120: Sum over outgoing excitation energies
c
c maxex  : maximum excitation energy bin for residual nucleus
c l2maxhf: 2 * lmaxhf
c lmaxhf : maximal l-value for transmission coefficients
c
c This loop is over all discrete levels and continuum bins of the
c final nucleus.
c
        do 120 nexout=0,maxex(Zix,Nix)
          l2maxhf=2*lmaxhf(type,nexout)
c
c Initialization of summations
c
c Nlast     : last discrete level
c Pprimebeg : start of residual parity summation
c parlev    : parity of level
c Pprimeend : end of residual parity summation
c Irspin2beg: 2 * start of residual spin summation
c jdis      : spin of level
c Irspin2end: 2 * end of residual spin summation
c J2res     : help variable
c maxJ      : maximal J-value
c
c For discrete states, the begin and end points of the residual
c spin/parity summation are both set equal to the residual discrete
c level spin/parity.
c
          if (nexout.le.Nlast(Zix,Nix,0)) then
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
c
c 130: Sum over residual parity
c
c pardif2: difference between residual and compound nucleus parity
c
c The variable pardif2 is used as an indicator of parity conservation
c for the outgoing channel.
c
          do 130 Pprime=Pprimebeg,Pprimeend,2
            pardif2=abs(parity-Pprime)/2
c
c 140: Sum over residual spin
c
c Irspin2    : 2 * residual spin
c Ir         : residual spin
c rho,rho0   : integrated level density
c flagastro  : flag for calculation of astrophysics reaction rate
c flagastrogs: flag for calculation of astrophysics reaction rate
c              with target in ground state only
c
            do 140 Irspin2=Irspin2beg,Irspin2end,2
              Ir=Irspin2/2
              enumhf(type,nexout,Ir,Pprime)=0.
              rho=rho0(type,nexout,Ir,Pprime)
              if (rho.lt.1.e-20) goto 140
              if (flagastro.and..not.flagastrogs) then
                J2maxastro=int(2*(lmaxhf(k0,nexout)+parspin(k0)+Ir))
                J2maxastro=min(J2maxastro,numJ)
              endif
c
c 150: Sum over j (jj2) of outgoing channel
c
c jj2primebeg: 2 * start of j' summation
c jj2primeend: 2 * end of j' summation
c jj2prime   : 2 * j'
c
              jj2primebeg=abs(J2-Irspin2)
              jj2primeend=J2+Irspin2
              do 150 jj2prime=jj2primebeg,jj2primeend,2
c
c 160: Sum over l of outgoing channel
c
c l2primebeg: 2 * start of l summation
c l2primeend: 2 * end of l summation
c l2prime   : 2 * l'
c lprime    : 2 * l
c modl      : help variable
c
                l2primebeg=abs(jj2prime-parspin2o)
                if (type.eq.0) l2primebeg=max(l2primebeg,2)
                l2primeend=jj2prime+parspin2o
                l2primeend=min(l2primeend,l2maxhf)
                do 160 l2prime=l2primebeg,l2primeend,2
                  lprime=l2prime/2
                  modl=mod(lprime,2)
c
c 1. Photons
c
c We include photons as a special case, with the multipole radiation
c selection rules (irad=0: M-transition, irad=1: E-transition)
c
c Tout: transmission coefficients
c Tgam: gamma transmission coefficients
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
c updown2 : spin index for transmission coefficient
c Tjlnex  : transmission coefficients as a function of particle type,
c           energy, spin and l-value
c factor  : help variable
c numtrans: number of transmission coefficients
c numhill : maximum number of Hill-Wheeler points
c
                    if (modl.ne.pardif2) goto 160
                    updown2=(jj2prime-l2prime)/pspin2o
                    Tout=Tjlnex(type,nexout,updown2,lprime)
                  endif
c
c The contribution is added to the total width.
c
                  factor=rho*Tout
                  denomhf=denomhf+factor
c
c Information needed for width fluctuation calculation. Values for
c for rho*T**i where (i=0,5) are stored. The photon contribution is
c stored in a single gamma width.
c
                  if (flagwidth) then
                    if (type.eq.0) then
                      gamwidth=gamwidth+factor
                    else
                      if (tnum.lt.numtrans-numhill-2) then
                        tnum=tnum+1
                        transjl(0,tnum)=rho
                        do 170 i=1,wpower
                          transjl(i,tnum)=0.
                          if (Tout.gt.1.e-30**(1./i)) transjl(i,tnum)=
     +                      rho*(Tout**i)
  170                   continue
                      endif
                    endif
                  else
c
c If NO width fluctuation corrections or angular distributions are
c required, the partial decay width is used in subroutine comptarget.
c
                    if (.not.flagcompang) enumhf(type,nexout,Ir,Pprime)=
     +                enumhf(type,nexout,Ir,Pprime)+factor
                  endif
c
c Astrophysical case
c
c Tastrotot  : total transmission coefficient for astrophysical case
c rhoastrotot: total level density for astrophysical case
c
                  if (flagastro.and..not.flagastrogs.and.type.eq.k0)
     +              then
                    if (Tout.gt.0.) then
                      rhoastrotot=rhoastrotot+rho
                      if (nexout.ne.Ltarget) Tastrotot=Tastrotot+factor
                    endif
                    if (type.eq.k0.and.nexout.ne.Ltarget) then
                      if (J2.le.J2maxastro) feed=feed+factor
                    endif
                  endif
  160           continue
  150         continue
  140       continue
  130     continue
  120   continue
  110 continue
c
c ******** Add fission and gamma transmission coefficients to
c          transmission coefficient array for width fluctuations *******
c
c 1. Fission.
c
c ihill   : counter for Hill-Wheeler magnitude
c tfisA   : transmission coefficient for Hill-Wheeler magnitude
c tfishill: help variable
c rhofisA : integrated level density corresponding to tfisA
c rhofis  : integrated level density
c
      if (flagwidth) then
        if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0) then
          do 190 ihill=1,numhill
            ratio=0.
            if (tfisA(J,parity,0).gt.0.)
     +        ratio=tfisA(J,parity,ihill)/tfisA(J,parity,0)
            tfishill=ratio*fiswidth
            tnum=tnum+1
            transjl(0,tnum)=max(rhofisA(J,parity,ihill),1.d0)
            do 200 i=1,wpower
              transjl(i,tnum)=0.
              if (tfishill.gt.1.e-30**(1./i)) transjl(i,tnum)=
     +          transjl(0,tnum)*(tfishill**i)/(transjl(0,tnum)**i)
  200       continue
  190     continue
        endif
c
c 2. Photons.
c
        transjl(0,tnum+1)=1.
        do 210 i=1,wpower
          transjl(i,tnum+1)=0.
          if (gamwidth.gt.1.e-30**(1./i).and.
     +      gamwidth.lt.1.e30**(1./i)) transjl(i,tnum+1)=gamwidth**i
  210   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
