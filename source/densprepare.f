      subroutine densprepare(Zcomp,Ncomp,idfis)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 4, 2012
c | Task  : Prepare energy grid, level density and transmission
c |         coefficient information for compound nucleus
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80     key
      integer          Zcomp,Ncomp,idfis,type,Zix,Nix,A,NL,odd,ldmod,
     +                 nexout,Pprime,Ir,l,irad,updown,nen,na,
     +                 nb,nc,ibar,ibk,ibin,J,parity
      real             Ex0plus,Ex0min,Efs,SS,Rodd,dEx,dExhalf,Exout,
     +                 Rboundary,Ex1min,Ex1plus,Eout,emax,emin,Exm,
     +                 Rspin,factl(0:numl),Egamma,factor,gn0,fstrength,
     +                 Ea,Eb,Ec,ta,tb,tc,tint,exfis,dExmin,elowest,elow,
     +                 emid
      double precision density
c
c ************ Determine energetically allowed transitions *************
c
c  10: Loop over particle types
c
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c idfis     : fission identifier
c Ex0plus   : upper boundary of entrance bin
c Exinc     : excitation energy of entrance bin
c dExinc    : excitation energy bin for mother nucleus
c Ex0min    : lower boundary of entrance bin
c Efs       : fast particle energy for gamma ray strength function
c SS,S      : separation energy per particle
c parskip   : logical to skip outgoing particle
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c AA,A      : mass number of residual nucleus
c Nlast,NL  : last discrete level
c odd       : odd (1) or even (0) nucleus
c Rodd      : term to determine integer or half-integer spins
c
c The mother excitation energy bin is characterized by the middle Exinc,
c the top Ex0plus and the bottom Ex0min.
c
      Ex0plus=Exinc+0.5*dExinc
      Ex0min=Exinc-0.5*dExinc
      Efs=Exinc-S(Zcomp,Ncomp,1)
      do 10 type=0,6
        if (parskip(type)) goto 10
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        A=AA(Zcomp,Ncomp,type)
        NL=Nlast(Zix,Nix,0)
        SS=S(Zcomp,Ncomp,type)
        odd=mod(A,2)
        Rodd=0.5*odd
c
c  20: Loop over outgoing excitation energies
c
c deltaEx,dEx: excitation energy bin for population arrays
c dExhalf    : half of excitation energy bin for population arrays
c nexmax     : maximum excitation energy bin for residual nucleus
c Exout,Ex   : excitation energy
c Rboundary  : factor taking into count first accessible mother bin for
c              discrete state
c
c There are 4 types of decay:
c 1. From discrete state to discrete state: This happens for
c    the primary compound nucleus, which is formed at an energy
c    Etotal and can decay to a discrete state of a residual nucleus
c    (e.g. compound elastic scattering).
c 2. From discrete state to continuum: This happens for
c    the primary compound nucleus, which is formed at an energy
c    Etotal and can decay to the continuum of a residual nucleus.
c    (e.g. continuum inelastic scattering).
c 3. From continuum to discrete state: This happens for multiple
c    emission, where a residual nucleus can be populated in a continuum
c    excitation energy bin which can decay to a discrete state of
c    another residual nucleus.
c 4. From continuum to continuum: This happens for multiple emission,
c    where a residual nucleus can be populated in a continuum
c    excitation energy bin which can decay to a continuum bin of
c    another residual nucleus.
c
c Types 3 and 4 are subject to boundary effects, i.e. they can represent
c cases where not the entire mother bin can have decayed to the residual
c bin or level, because of the particle separation energy. The end
c points need to be taken care of by a proper normalization.
c
        do 20 nexout=0,nexmax(type)
          dEx=deltaEx(Zix,Nix,nexout)
          dExhalf=0.5*dEx
          Exout=Ex(Zix,Nix,nexout)
          Rboundary=1.
c
c Types 1 and 2. Decay from the primary compound nucleus.
c No special care needs to be taken. The residual bin/level is
c characterized by an excitation energy Exout. For transitions to the
c continuum, the bins is further characterized by a top Ex1plus and a
c bottom Ex1min.
c
c primary: flag to designate primary (binary) reaction
c Ex1min : lower boundary of residual bin
c Ex1plus: upper boundary of residual bin
c Eout   : emission energy
c
          if (primary) then
            if (nexout.gt.NL) then
              Ex1min=Exout-dExhalf
              Ex1plus=Exout+dExhalf
            endif
            Eout=Exinc-Exout-SS
          else
c
c Type 3. Decay from continuum to continuum. For most residual continuum
c bins, no special care needs to be taken and the emission energy Eout
c that characterizes the transition is simply the average between the
c highest energetic transition that is possible (emax, from the top of
c the mother bin to the bottom of the residual bin) and the lowest
c (emin). However, the highest residual bin (nexout=nexmax) is
c characterized by different energies (Ex1plus is the maximal residual
c excitation energy and Exout is shifted from its original position).
c If some transitions from mother to residual bin are forbidden, the
c factor Rboundary takes care of the correction.
c
c emax: maximal emission energy within bin decay
c emin: minimal emission energy within bin decay
c
            if (nexout.gt.NL) then
              Ex1min=Exout-dExhalf
              if (nexout.eq.nexmax(type).and.type.ge.1) then
                Ex1plus=Ex0plus-SS
                Exout=0.5*(Ex1plus+Ex1min)
              else
                Ex1plus=Exout+dExhalf
              endif
              emax=Ex0plus-SS-Ex1min
              emin=Ex0min-SS-Ex1plus
              Eout=0.5*(emin+emax)
              if (emin.lt.0.) then
                if (Eout.gt.0.) then
                  Rboundary=1.-0.5*(emin/(0.5*(emax-emin)))**2
                else
                  Rboundary=0.5*(emax/(0.5*(emax-emin)))**2
                endif
                emin=0.
              endif
              Eout=0.5*(emin+emax)
            else
c
c Type 4. Decay from continuum to discrete. The lowest possible mother
c excitation bin can not entirely decay to the discrete state.
c For the residual discrete state, it is checked whether the mother
c excitation bin is such a boundary case. This is done by adding the
c particle separation energy to the excitation energy of the residual
c discrete state. The correction is put in Rboundary.
c
c Exm: help variable
c
              Exm=Exout+SS
              if (Exm.le.Ex0plus.and.Exm.gt.Ex0min) then
                Rboundary=(Ex0plus-Exm)/dExinc
                Eout=0.5*(Ex0plus+Exm)-SS-Exout
              else
                Eout=Exinc-SS-Exout
              endif
            endif
          endif
c
c ********************** Determine level densities *********************
c
c The calculation of level densities can be done outside many loops of
c various quantum numbers performed in other subroutines. Therefore,
c we store the level density as function of residual nucleus (type),
c excitation energy (nexout), spin (Ir) and parity (Pprime) in the
c array rho0.
c
c Pprime  : parity
c parlev  : parity of level
c Ir,Rspin: residual spin
c jdis    : spin of level
c rho0    : integrated level density
c maxJ    : maximal J-value
c rhogrid : integrated level density
c
c For discrete states, the level density is set to Rboundary.
c
          if (nexout.le.NL) then
            Pprime=parlev(Zix,Nix,nexout)
            Ir=int(jdis(Zix,Nix,nexout))
            rho0(type,nexout,Ir,Pprime)=Rboundary
          else
c
c For decay to the continuum we use a spin and parity dependent level
c density.
c
            do 110 Pprime=-1,1,2
              do 120 Ir=0,maxJ(Zix,Nix,nexout)
                rho0(type,nexout,Ir,Pprime)=
     +            Rboundary*rhogrid(Zix,Nix,nexout,Ir,Pprime)
  120         continue
  110       continue
          endif
c
c Optional adjustment factors
c
c tljadjust: logical for energy-dependent Tlj adjustment
c tladjust : adjustable factor for Tlj (default 1.)
c adjust   : subroutine for energy-dependent parameter adjustment
c factl    : multiplication factor
c
          if (tljadjust(type)) then
            do 130 l=0,numl
              key='tljadjust'
              call adjust(Exout,key,0,0,type,l,factor)
              factl(l)=factor*tladjust(type,l)
  130       continue
          else
            do 140 l=0,numl
              factl(l)=1.
  140       continue
          endif
c
c ************* Interpolation of transmission coefficients *************
c
c 1. Gamma transmission coefficients
c
c lmaxhf   : maximal l-value for transmission coefficients
c gammax   : number of l-values for gamma multipolarity
c Egamma   : gamma energy
c irad     : variable to indicate M(=0) or E(=1) radiation
c Tgam     : gamma transmission coefficients
c twopi    : 2.*pi
c gamadjust: logical for energy-dependent gamma adjustment
c adjust   : subroutine for energy-dependent parameter adjustment
c factor   : multiplication factor
c gnorm    : gamma normalization factor
c fstrength: gamma ray strength function
c
          if (type.eq.0) then
            lmaxhf(0,nexout)=gammax
            Egamma=Exinc-Exout
            do 210 l=0,gammax
              do 210 irad=0,1
                Tgam(nexout,l,irad)=0.
  210       continue
            if (Egamma.le.0) goto 20
            if (gamadjust(Zcomp,Ncomp)) then
              key='gnorm'
              call adjust(Ecomp,key,0,0,0,0,factor)
              gn0=factor*gnorm
            else
              gn0=gnorm
            endif
            do 220 l=1,gammax
              do 220 irad=0,1
                Tgam(nexout,l,irad)=twopi*(Egamma**(2*l+1))*gn0*
     +            fstrength(Zcomp,Ncomp,Efs,Egamma,irad,l)*factl(l)
  220       continue
          else
c
c 2. Particle transmission coefficients
c
c updown    : spin index for transmission coefficient
c Tjl,Tjlnex: transmission coefficients as a function of
c             particle type, energy, spin and l-value
c Tl,Tlnex  : transmission coefficients as a function of
c             particle type, energy and l-value (averaged over spin)
c egrid     : outgoing energy grid
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c locate    : subroutine to find value in ordered table
c na,Ea,ta,.: help variables
c pol2      : subroutine for interpolation of second order
c tint      : help variable
c transeps  : absolute limit for transmission coefficient
c
            do 230 updown=-1,1
              do 230 l=0,numl
                Tjlnex(type,nexout,updown,l)=0.
  230       continue
            do 240 l=0,numl
              Tlnex(type,nexout,l)=0.
  240       continue
            lmaxhf(type,nexout)=0
            if (Eout.lt.egrid(ebegin(type))) goto 20
            if (ebegin(type).ge.eend(type)) goto 20
c
c To get the transmission coefficients on the excitation energy grid,
c Tjlnex, from those on the emission energy grid, Tjl, we use
c interpolation of the second order.
c
            call locate(egrid,ebegin(type),eend(type),Eout,nen)
            if (nen.gt.ebegin(type)+1) then
              na=nen-1
              nb=nen
              nc=nen+1
            else
              na=nen
              nb=nen+1
              nc=nen+2
            endif
            Ea=egrid(na)
            Eb=egrid(nb)
            Ec=egrid(nc)
            do 260 updown=-1,1
              do 260 l=0,lmax(type,nen)
                ta=Tjl(type,na,updown,l)
                tb=Tjl(type,nb,updown,l)
                tc=Tjl(type,nc,updown,l)
                call pol2(Ea,Eb,Ec,ta,tb,tc,Eout,tint)
                if (tint.lt.transeps) tint=0.
                Tjlnex(type,nexout,updown,l)=factl(l)*tint
  260       continue
            if (.not.flagfullhf) then
              do 270 l=0,lmax(type,nen)
                ta=Tl(type,na,l)
                tb=Tl(type,nb,l)
                tc=Tl(type,nc,l)
                call pol2(Ea,Eb,Ec,ta,tb,tc,Eout,tint)
                if (tint.lt.transeps) tint=0.
                Tlnex(type,nexout,l)=factl(l)*tint
  270         continue
            endif
c
c The maximal l-values needed in the compound nucleus calculations are
c determined.
c
c lmaxinc    : maximal l-value for transmission coefficients
c flaginitpop: flag for initial population distribution
c
            lmaxhf(type,nexout)=lmax(type,nen)
          endif
   20   continue
        if (nexmax(type).gt.0) lmaxhf(type,nexmax(type))=
     +    lmaxhf(type,nexmax(type)-1)
   10 continue
      if (flaginitpop) then
        lmaxhf(k0,0)=gammax
      else
        lmaxhf(k0,0)=lmaxinc
      endif
c
c **** Calculate fission level densities and Hill-Wheeler terms ********
c
c flagfission       : flag for fission
c nfisbar           : number of fission barrier parameters
c exfis             : help variable
c fecont            : start of continuum energy
c nbintfis,numbinfis: number of bins
c dExmin,dEx        : energy bin
c elowest,elow,emid : help variables
c dExhalf           : help variable
c eintfis           : excitation energy for fission
c parity            : parity
c numJ              : maximal J-value
c rhofis            : integrated level density
c density           : level density
c
      if ((flagfission).and.(idfis.eq.1)) then
        A=AA(Zcomp,Ncomp,0)
        odd=mod(A,2)
        Rodd=0.5*odd
        if (nfisbar(Zcomp,Ncomp).ne.0) then
          ldmod=ldmodel(Zcomp,Ncomp)
          do 310 ibar=1,nfisbar(Zcomp,Ncomp)
            if (primary) then
              exfis=Exinc-fecont(Zcomp,Ncomp,ibar)
            else
              exfis=Ex(Zcomp,Ncomp,maxex(Zcomp,Ncomp))-
     +          fecont(Zcomp,Ncomp,ibar)
            endif
            nbintfis(ibar)=numbinfis/2
            if (exfis.le.0.) goto 310
            dExmin=0.01
            dEx=exfis/nbintfis(ibar)
            if (dEx.lt.dExmin) then
              nbintfis(ibar)=max(int(exfis/dExmin),1)
              dEx=exfis/nbintfis(ibar)
            endif
            ibk=1
            dExhalf=0.5*dEx
            elowest=fecont(Zcomp,Ncomp,ibar)
            do 320 ibin=0,nbintfis(ibar)-1
              elow=elowest+ibin*dEx
              emid=elow+dExhalf
              eintfis(ibk,ibar)=elow
              eintfis(ibk+1,ibar)=emid
              do 330 parity=-1,1,2
                do 340 J=0,numJ
                  Rspin=J+Rodd
                  rhofis(ibk,J,parity,ibar)=
     +              density(Zcomp,Ncomp,elow,Rspin,parity,ibar,ldmod)
                  rhofis(ibk+1,J,parity,ibar)=
     +              density(Zcomp,Ncomp,emid,Rspin,parity,ibar,ldmod)
  340           continue
  330         continue
              ibk=ibk+2
  320       continue
            nbintfis(ibar)=ibk-1
            eintfis(nbintfis(ibar),ibar)=exfis+elowest
  310     continue
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
