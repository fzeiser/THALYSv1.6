      subroutine exgrid(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 10, 2013
c | Task  : Set excitation energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,Zdeep,Ndeep,type,Zix,Nix,Zmother,
     +                 Nmother,NL,i,A,odd,nex,Aix,nexbins,ldmod,Pbeg,
     +                 Pprime,Ir
      real             Edif,Eup(0:numex),eb,ee,dEx,Rodd,Exout,
     +                 Ex1min,Ex1plus,ald,spincut,Rspin
      double precision rho1,rho2,rho3,density,r1log,r2log,r3log
c
c ********************* Set maximum excitation energy ******************
c
c flagomponly: flag to execute ONLY an optical model calculation
c flagcomp   : flag for compound nucleus calculation
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c Zdeep      : charge number index for lightest residual nucleus
c Ndeep      : neutron number index for lightest residual nucleus
c parskip    : logical to skip outgoing particle
c Zindex     : charge number index for residual nucleus
c Nindex     : neutron number index for residual nucleus
c Exmax      : maximum excitation energy for residual nucleus
c Zmother    : charge number index for mother nucleus
c parZ       : charge number of particle
c Nmother    : neutron number index for mother nucleus
c parN       : neutron number of particle
c Exmax0     : maximum excitation energy for residual nucleus (including
c              negative energies)
c S          : separation energy per particle
c
c All possible routes to all residual nuclei are followed, to
c determine the maximum possible excitation energy for each nucleus,
c given the incident energy.
c
      if (flagomponly.and..not.flagcomp) return
      Zdeep=Zcomp
      Ndeep=Ncomp
      do 10 type=1,6
        if (parskip(type)) goto 10
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        Zdeep=max(Zdeep,Zix)
        Ndeep=max(Ndeep,Nix)
   10 continue
      do 20 Zix=Zcomp,Zdeep
        do 30 Nix=0,Ndeep
          if (Zix.eq.Zcomp.and.Nix.eq.Ncomp) goto 30
          if (Exmax(Zix,Nix).ne.0.) goto 30
          do 40 type=1,6
            if (parskip(type)) goto 40
            Zmother=Zix-parZ(type)
            if (Zmother.lt.0) goto 40
            Nmother=Nix-parN(type)
            if (Nmother.lt.0) goto 40
            Exmax0(Zix,Nix)=Exmax0(Zmother,Nmother)-
     +        S(Zmother,Nmother,type)
            Exmax(Zix,Nix)=max(Exmax0(Zix,Nix),0.)
   40     continue
   30   continue
   20 continue
c
c ******* Define excitation energy grid for residual nuclei ************
c
c The excitation energies are given by the array Ex. The first NL values
c of Ex correspond to the discrete level excitation energies of the
c residual nucleus specified by (Zix,Nix). The NL+1th value corresponds
c to the first continuum energy bin. The continuum bins have a width of
c deltaEx. The continuum part of the target and initial compound nucleus
c are divided into nbins equidistant energy bins where nbins was given
c in the input file or set by default. The continuum parts of the
c other residual nuclides are also divided in equidistant bins. For the
c first generation of nuclides (within 4 mass units of the initial
c compound nucleus), the number of bins is equal to that of target. For
c nuclides more than 8 mass units away, half the number of bins is
c chosen. For intermediate nuclides, an interpolated number is adopted.
c In sum, for each nucleus the excitation energy range is completely
c filled by equidistant bins. The bin widths thus gradually change.
c
c Nlast,NL   : last discrete level
c maxex      : maximum excitation energy bin for compound nucleus
c deltaEx,dEx: excitation energy bin for population arrays
c Qres       : Q-value for residual nucleus
c targetE    : energy of target
c Ltarget    : excited level of target
c k0         : index of incident particle
c Etotal     : total energy of compound system (target + projectile)
c Ethresh    : threshold incident energy for residual nucleus
c edis       : energy of level
c specmass   : specific mass for target nucleus
c
c The Q-value for residual nuclides is determined, both for the ground
c state and for possible isomers.
c
      do 110 type=0,6
        if (parskip(type)) goto 110
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        NL=Nlast(Zix,Nix,0)
        if (maxex(Zix,Nix).ne.0) goto 110
        do 120 nex=0,numex
          deltaEx(Zix,Nix,nex)=0.
  120   continue
        if (Qres(Zix,Nix,0).eq.0.) then
          Edif=Exmax0(Zix,Nix)-Etotal
          Qres(Zix,Nix,0)=S(0,0,k0)+targetE+Edif
        endif
        if (Ltarget.ne.0.and.Zix.eq.parZ(k0).and.Nix.eq.parN(k0))
     +    Qres(Zix,Nix,0)=targetE
        do 130 nex=0,NL
          if (Ethresh(Zix,Nix,nex).eq.0.) then
            Qres(Zix,Nix,nex)=Qres(Zix,Nix,0)-edis(Zix,Nix,nex)
            Ethresh(Zix,Nix,nex)=-(Qres(Zix,Nix,nex)/
     +        specmass(parZ(k0),parN(k0),k0))
            Ethresh(Zix,Nix,nex)=max(Ethresh(Zix,Nix,nex),0.d0)
          endif
  130   continue
c
c The highest possible excitation energy could be a discrete state.
c
c Ex: excitation energy
c
        do 140 nex=0,NL
          if (nex.gt.0.and.edis(Zix,Nix,nex).ge.Exmax(Zix,Nix)) then
            maxex(Zix,Nix)=nex-1
            goto 180
          endif
          Ex(Zix,Nix,nex)=edis(Zix,Nix,nex)
          if (nex.gt.0) deltaEx(Zix,Nix,nex)=
     +      0.5*(edis(Zix,Nix,min(NL,nex+1))-edis(Zix,Nix,nex-1))
  140   continue
c
c Division of the continuum into bins.
c
c Aix          : mass number index for residual nucleus
c nexbins,nbins: number of continuum excitation energy bins
c flagequi     : flag to use equidistant excitation bins instead of
c                logarithmic bins
c nexmax       : maximum excitation energy bin for residual nucleus
c
        Aix=Zix+Nix
        if (Aix.le.4) then
          nexbins=nbins
        else
          if (Aix.le.8) then
            nexbins=int(real(1-0.1*(Aix-4))*nbins)
          else
            nexbins=nbins/2
          endif
        endif
        nexbins=max(nexbins,2)
        eb=Ex(Zix,Nix,NL)
        ee=max(Exmax(Zix,Nix),eb+0.001)
        if (flagequi.or.eb.eq.0.) then
          do 150 i=0,nexbins
            Eup(i)=eb+real(i)/nexbins*(ee-eb)
  150     continue
        else
          eb=log(eb)
          ee=log(ee)
          do 160 i=0,nexbins
            Eup(i)=exp(eb+real(i)/nexbins*(ee-eb))
  160     continue
        endif
        do 170 i=1,nexbins
          nex=NL+i
          Ex(Zix,Nix,nex)=0.5*(Eup(i-1)+Eup(i))
          deltaEx(Zix,Nix,nex)=Eup(i)-Eup(i-1)
  170   continue
        maxex(Zix,Nix)=NL+nexbins
        if (Zix.eq.0.and.Nix.eq.0) Ex(0,0,maxex(0,0)+1)=Etotal
  180   nexmax(type)=maxex(Zix,Nix)
c
c ****** Determine level densities on basic excitation energy grid *****
c
c odd     : odd (1) or even (0) nucleus
c AA,A    : mass number of residual nucleus
c Rodd    : term to determine integer or half-integer spins
c ald     : level density parameter
c ldmodel : level density model
c Ex,Exout: excitation energy
c Ex1min  : lower boundary of residual bin
c Ex1plus : upper boundary of residual bin
c maxJ    : maximal J-value
c spincut : spin cutoff factor
c numJ    : maximal J-value
c
c The calculation of level densities can be done outside many loops of
c various quantum numbers performed in other subroutines. Therefore,
c we store the level density as function of residual nucleus, excitation
c energy (nex), spin (Ir) and parity (Pprime) in the array rhogrid.
c
        A=AA(Zcomp,Ncomp,type)
        odd=mod(A,2)
        Rodd=0.5*odd
        ald=real(A)/8.
        ldmod=ldmodel(Zix,Nix)
        do 210 nex=NL+1,maxex(Zix,Nix)
          dEx=deltaEx(Zix,Nix,nex)
          Exout=Ex(Zix,Nix,nex)
          Ex1min=Exout-0.5*dEx
          Ex1plus=Exout+0.5*dEx
c
c Here we define the maximum J that can be reached in a given excitation
c energy bin. By default, this maxJ value is set to 3*sigma where
c sigma is the square root of the spin cut-off parameter of the level
c density.
c
         maxJ(Zix,Nix,nex)=max(3.*sqrt(spincut(Zix,Nix,ald,Exout,0)),5.)
         maxJ(Zix,Nix,nex)=min(maxJ(Zix,Nix,nex),numJ)
c
c In the compound nucleus subroutines, the particle widths are
c determined by means of products of level densities and transmission
c coefficients. Instead of taking this product exactly at the middle of
c the excitation energy bins, we get a better numerical result by
c performing a logarithmic average over the bin for the level density,
c using the middle, top and bottom.
c
c flagparity: flag for non-equal parity distribution
c Pbeg      : help variable
c Pprime    : parity
c Ir,Rspin  : residual spin
c rho1-3    : help variables
c density   : level density
c r1log,..  : help variables
c rhogrid   : integrated level density
c
c For an equiprobable parity distribution for level densities, the loop
c over Pprime only needs to be performed once and the result for the
c level density is equal for both parities.
c
          if (flagparity) then
            Pbeg=-1
          else
            Pbeg=1
          endif
          do 220 Pprime=Pbeg,1,2
            do 230 Ir=0,maxJ(Zix,Nix,nex)
              Rspin=real(Ir)+Rodd
              rho1=density(Zix,Nix,Ex1min,Rspin,Pprime,0,ldmod)*
     +          (1.+1.d-10)
              rho2=density(Zix,Nix,Exout,Rspin,Pprime,0,ldmod)
              rho3=density(Zix,Nix,Ex1plus,Rspin,Pprime,0,ldmod)*
     +          (1.+1.d-10)
              r1log=log(rho1)
              r2log=log(rho2)
              r3log=log(rho3)
              if (r2log.ne.r1log.and.r2log.ne.r3log) then
                rhogrid(Zix,Nix,nex,Ir,Pprime)=0.5*dEx*
     +            ((rho1-rho2)/(r1log-r2log)+
     +            (rho2-rho3)/(r2log-r3log))
              else
                rhogrid(Zix,Nix,nex,Ir,Pprime)=dEx*rho2
              endif
              if (.not.flagparity)
     +          rhogrid(Zix,Nix,nex,Ir,-1)=rhogrid(Zix,Nix,nex,Ir,1)
  230       continue
  220     continue
  210   continue
  110 continue
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
