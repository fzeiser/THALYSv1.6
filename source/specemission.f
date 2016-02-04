      subroutine specemission(Zcomp,Ncomp,nex,idorg,type,nexout)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 21, 2012
c | Task  : Exclusive emission cross sections for continuum
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,nex,idorg,type,nexout,Zix,Nix,nexout2,nen,NL,
     +        nenbeg,nenend
      real    Ex0plus,Ex0min,SS,Exm,dEx,Exmin,Exout,Ex1min,Ex1plus,emax,
     +        emin,Eout,dE,dEhalf,feed0,xso00,feedm,xsom0,feedp,xsop0,
     +        xso0m,xsopm,xsomp,xso0p,edist,term,term1,term2,
     +        weight(numen),fracbot,fractop,sumweight,factor
c
c ****************** Smearing into emission spectra ********************
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c specemis   : exclusive emission contribution
c popexcl    : population cross section of bin just before decay
c speceps    : limit for cross section spectra
c Exinc      : excitation energy of entrance bin
c Ex,Exout   : excitation energy
c dExinc     : excitation energy bin for mother nucleus
c deltaEx,dEx: excitation energy bin for population arrays
c Ex0plus    : upper boundary of entrance bin
c Ex0min     : lower boundary of entrance bin
c SS,S       : separation energy
c Exm        : maximal attainable energy
c egrid      : outgoing energy grid
c ebegin     : first energy point of energy grid
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c maxex      : maximum excitation energy bin for compound nucleus
c Exmin      : help variable
c nexmax     : maximum excitation energy bin for residual nucleus
c
c The decay from compound to residual nuclei is converted from the
c excitation energy grid to emission energies.
c The spectrum is obtained by spreading the decay over the mother
c bin and, in the case of continuum-continuum transitions, the
c residual bin.
c
c 1. Loop over excitation energy bins of compound nucleus.
c
c For each mother excitation energy bin, determine the highest
c possible excitation energy bin nexmax for the residual nuclei.
c As reference, we take the top of the mother bin. The maximal
c residual excitation energy Exm is obtained by subtracting the
c separation energy from this. The bin in which this Exm falls
c is denoted by nexmax.
c
      do 10 nen=0,numen
        specemis(nen)=0.
   10 continue
      if (popexcl(Zcomp,Ncomp,nex).le.speceps) return
      Exinc=Ex(Zcomp,Ncomp,nex)
      dExinc=deltaEx(Zcomp,Ncomp,nex)
      Ex0plus=Exinc+0.5*dExinc
      Ex0min=Exinc-0.5*dExinc
      SS=S(Zcomp,Ncomp,type)
      Exm=Ex0plus-SS
      if (type.gt.1) Exm=Exm-egrid(ebegin(type))
      Zix=Zindex(Zcomp,Ncomp,type)
      Nix=Nindex(Zcomp,Ncomp,type)
      do 20 nexout2=maxex(Zix,Nix),0,-1
        dEx=deltaEx(Zix,Nix,nexout2)
        Exmin=Ex(Zix,Nix,nexout2)-0.5*dEx
        if (Exmin.lt.Exm) then
          nexmax(type)=nexout2
          goto 100
        endif
   20 continue
      return
c
c 2. Determine the widths over which the decay must be spread.
c
c NL,Nlast: last discrete level
c Ex1min  : lower boundary of residual bin
c Ex1plus : upper boundary of residual bin
c emax    : maximal emission energy within bin decay
c emin    : minimal emission energy within bin decay
c Eout    : emission energy
c
  100 Exout=Ex(Zix,Nix,nexout)
      dEx=deltaEx(Zix,Nix,nexout)
      NL=Nlast(Zix,Nix,0)
c
c Decay from continuum to continuum. For most residual continuum
c bins, no special care needs to be taken and the emission energy Eout
c that characterizes the transition is simply the average between the
c highest energetic transition that is possible (emax, from the top of
c the mother bin to the bottom of the residual bin) and the lowest
c (emin). However, the highest residual bin (nexout=nexmax) is
c characterized by different energies (Ex1plus is the maximal residual
c excitation energy and Eout is again the average between emin and
c emax).
c
      if (nexout.gt.NL) then
        Ex1min=Exout-0.5*dEx
        if (nexout.eq.nexmax(type).and.type.ge.1) then
          Ex1plus=Ex0plus-SS
          Exout=0.5*(Ex1plus+Ex1min)
        else
          Ex1plus=Exout+0.5*dEx
        endif
        emax=Ex0plus-SS-Ex1min
        emin=max(Ex0min-SS-Ex1plus,0.)
        Eout=0.5*(emin+emax)
      else
c
c Decay from continuum to discrete. The lowest possible mother
c excitation bin can not entirely decay to the discrete state.
c For the residual discrete state, it is checked whether the mother
c excitation bin is such a boundary case. This is done by adding the
c particle separation energy to the excitation energy of the residual
c discrete state.
c
        Exm=Exout+SS
        if (Exm.le.Ex0plus.and.Exm.gt.Ex0min) then
          Eout=0.5*(Ex0plus+Exm)-SS-Exout
          emin=0.
        else
          Eout=Exinc-SS-Exout
          emin=Ex0min-SS-Exout
        endif
        emax=Ex0plus-SS-Exout
      endif
c
c 3. Redistribution of decay from population on emission energy grid.
c
c After the bottom (emin) and top (emax) of the excitation energy
c bin have been determined. Then the begin and end points for the
c spectrum are located.
c
c locate       : subroutine to find value in ordered table
c Ebottom      : bottom of outgoing energy bin
c eend         : last energy point of energy grid
c nenbeg,nenend: help variables
c dE           : emission bin
c dEhalf       : half of emission bin
c
      call locate(Ebottom,ebegin(type),eend(type),emin,nenbeg)
      call locate(Ebottom,ebegin(type),eend(type),emax,nenend)
      nenbeg=max(nenbeg,1)
      if (nenend.lt.nenbeg) return
      dE=emax-emin
      dEhalf=0.5*dE
c
c For the interpolation, we first determine the feeding from adjacent
c mother bins. Note that we interpolate the whole product of terms
c from the exclusive cross section that depend on excitation energy.
c
c feed0,.......: help variable
c xsexcl       : exclusive cross section per excitation energy
c xso00,....   : cross section on excitation energy grid
c feedexcl     : feeding terms from compound excitation energy bin to
c                residual excitation energy bin
c feedm,feedp..: help variables
c
      feed0=0.
      if (popexcl(Zcomp,Ncomp,nex).ne.0.)
     + feed0=xsexcl(idorg,nex)/popexcl(Zcomp,Ncomp,nex)
      xso00=feedexcl(Zcomp,Ncomp,type,nex,nexout)*feed0
      if (xso00.eq.0.) return
      feedm=0.
      if (popexcl(Zcomp,Ncomp,nex-1).ne.0.)
     +  feedm=xsexcl(idorg,nex-1)/popexcl(Zcomp,Ncomp,nex-1)
      xsom0=feedexcl(Zcomp,Ncomp,type,nex-1,nexout)*feedm
      feedp=0.
      if (popexcl(Zcomp,Ncomp,nex+1).ne.0.)
     +  feedp=xsexcl(idorg,nex+1)/popexcl(Zcomp,Ncomp,nex+1)
      xsop0=feedexcl(Zcomp,Ncomp,type,nex+1,nexout)*feedp
c
c Decay from continuum to continuum.
c We need to interpolate between adjacent residual bins. They are also
c determined. For the lowest emission energies we ensure that the limit
c is zero for zero outgoing energy.
c
      if (nexout.gt.NL) then
        xso0m=feedexcl(Zcomp,Ncomp,type,nex,nexout-1)*feed0
        xsopm=feedexcl(Zcomp,Ncomp,type,nex+1,nexout-1)*feedp
        xsomp=feedexcl(Zcomp,Ncomp,type,nex-1,nexout+1)*feedm
        xso0p=feedexcl(Zcomp,Ncomp,type,nex,nexout+1)*feed0
        do 120 nen=nenbeg,nenend
c
c A. Emission energy lower than the average emission energy.
c    Interpolate with lower bins.
c
c edist       : help variable
c term,term1..: help variables
c
          if (egrid(nen).le.Eout) then
            edist=Eout-egrid(nen)
            if (emin.le.0.0001) then
              term=xso00*(1.-edist/dEhalf)
            else
              term1=xso00+0.5*edist/dE*(xsom0-xso00)
              term2=xso0p+0.5*edist/dE*(xsomp-xso0p)
              term=term1+0.5*edist/dE*(term2-term1)
            endif
          else
c
c B. Emission energy higher than the average emission energy.
c    Interpolate with higher bins.
c
            edist=egrid(nen)-Eout
            if (xsop0.eq.0.) then
              term1=xso00*(1.-0.5*edist/dEhalf)
              term2=0.
            else
              term1=xso00+0.5*edist/dE*(xsop0-xso00)
              term2=xso0m+0.5*edist/dE*(xsopm-xso0m)
            endif
            term=term1+0.5*edist/dE*(term2-term1)
          endif
c
c Each emission energy that is possible within the transition from bin
c to bin gets a weight.
c
c weight: weight of emission energy bin
c deltaE: energy bin around outgoing energies
c
          term=max(term,0.)
          weight(nen)=term*deltaE(nen)
  120   continue
      else
c
c Decay from continuum to discrete.
c We only need to interpolate between adjacent mother bins.
c
        do 130 nen=nenbeg,nenend
          if (egrid(nen).le.Eout) then
            edist=Eout-egrid(nen)
            if (Exm.le.Ex0plus.and.Exm.gt.Ex0min) then
              term=xso00*(1.-edist/dEhalf)
            else
              term=xso00+edist/dExinc*(xsom0-xso00)
            endif
          else
            edist=egrid(nen)-Eout
            if (xsop0.eq.0.) then
              term=xso00*(1.-edist/dEhalf)
            else
              term=xso00+edist/dExinc*(xsop0-xso00)
            endif
          endif
          term=max(term,0.)
          weight(nen)=term*deltaE(nen)
  130   continue
      endif
c
c The end points of the bin are taken into account and this is
c corrected in weight. There are also cases where the excitation
c energy bin falls completely in the emission energy bin
c (nenbeg=nenend).
c
c fracbot,fractop: help variables
c Etop           : top of outgoing energy bin
c sumweight      : help variable
c
      fracbot=(emin-Ebottom(nenbeg))/deltaE(nenbeg)
      fractop=(Etop(nenend)-emax)/deltaE(nenend)
      if (nexout.eq.NL+1) fractop=0.
      if (nenend.eq.nenbeg) then
        weight(nenbeg)=weight(nenbeg)*(1.-fracbot-fractop)
      else
        weight(nenbeg)=weight(nenbeg)*(1.-fracbot)
        weight(nenend)=weight(nenend)*(1.-fractop)
      endif
c
c The weights are summed to ensure that we get exact flux conservation
c by normalizing with sumweight.
c
      sumweight=0.
      do 140 nen=nenbeg,nenend
        sumweight=sumweight+weight(nen)
  140 continue
      if (sumweight.eq.0.) return
c
c The contribution for mother bin --> residual bin is added to the
c spectrum.
c
c factor: help variable
c
      do 150 nen=nenbeg,nenend
        factor=xso00/sumweight
        specemis(nen)=factor*weight(nen)/deltaE(nen)
  150 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
