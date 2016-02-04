      subroutine compemission(Zcomp,Ncomp,nex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : July 26, 2013
c | Task  : Compound emission spectra
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,nex,type,Zix,Nix,nexout,NL,nen,nenbeg,nenend
      real    Ex0plus,Ex0min,dEx,Exm,Exmin,SS,xsgrid(numen),
     +        xsmpegrid(numen),Exout,Ex1min,Ex1plus,emin,emax,Eout,
     +        dE,dEhalf,xso00,xsom0,xsop0,xso0m,xsopm,xsomp,xso0p,edist,
     +        term,term1,term2,weight(numen),fracbot,fractop,sumweight,
     +        fracmpe,factor,emisfac
c
c ****************** Smearing into emission spectra ********************
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c maxex      : maximum excitation energy bin for compound nucleus
c Nlast,NL   : last discrete level
c Exinc      : excitation energy of entrance bin
c Ex,Exout   : excitation energy
c Ex0plus    : upper boundary of entrance bin
c dExinc     : excitation energy bin for mother nucleus
c Ex0min     : lower boundary of entrance bin
c parskip    : logical to skip outgoing particle
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c Exm        : maximal attainable energy
c SS,S       : separation energy
c egrid      : outgoing energy grid
c ebegin     : first energy point of energy grid
c deltaEx,dEx: excitation energy bin for population arrays
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
      do 10 nex=maxex(Zcomp,Ncomp),Nlast(Zcomp,Ncomp,0)+1,-1
        Exinc=Ex(Zcomp,Ncomp,nex)
        dExinc=deltaEx(Zcomp,Ncomp,nex)
        Ex0plus=Exinc+0.5*dExinc
        Ex0min=Exinc-0.5*dExinc
        do 20 type=0,6
          if (parskip(type)) goto 20
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          Exm=Ex0plus-S(Zcomp,Ncomp,type)
          if (type.gt.1) Exm=Exm-egrid(ebegin(type))
          do 30 nexout=maxex(Zix,Nix),0,-1
            dEx=deltaEx(Zix,Nix,nexout)
            Exmin=Ex(Zix,Nix,nexout)-0.5*dEx
            if (Exmin.lt.Exm) then
              nexmax(type)=nexout
              goto 20
            endif
   30     continue
          nexmax(type)=-1
   20   continue
        nexmax(0)=nex-1
c
c 2. Determine the widths over which the decay must be spread.
c
c eend         : last energy point of energy grid
c numen        : maximum number of outgoing energies
c xsgrid,xsemis: emission spectrum from compound nucleus
c xsmpegrid    : multiple-preequilibrium emission spectrum from
c                compound nucleus
c
        do 110 type=0,6
          if (parskip(type)) goto 110
          if (nexmax(type).lt.0) goto 110
          if (ebegin(type).ge.eend(type)) goto 110
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          NL=Nlast(Zix,Nix,0)
          SS=S(Zcomp,Ncomp,type)
c
c end of recoil initialisation
c
          do 120 nen=1,numen
            xsgrid(nen)=0.
            xsmpegrid(nen)=0.
  120     continue
c
c Loop over outgoing excitation energy bins
c
c mcontrib: contribution to emission spectrum
c speceps : limit for cross section spectra
c
          do 130 nexout=0,nexmax(type)
            dEx=deltaEx(Zix,Nix,nexout)
            if (nexout.eq.0.and.NL.eq.0) goto 130
            if (mcontrib(type,nex,nexout).le.speceps) goto 130
            Exout=Ex(Zix,Nix,nexout)
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
c Ex1min : lower boundary of residual bin
c Ex1plus: upper boundary of residual bin
c emax   : maximal emission energy within bin decay
c emin   : minimal emission energy
c Eout   : emission energy
c
            if (nexout.gt.NL) then
              Ex1min=Exout-0.5*dEx
              if (nexout.eq.nexmax(type).and.type.ge.1) then
                Ex1plus=Ex0plus-SS
                Exout=0.5*(Ex1plus+Ex1min)
              else
                Ex1plus=Exout+0.5*dEx
              endif
              emin=Ex0min-SS-Ex1plus
              emin=max(emin,0.)
              emax=Ex0plus-SS-Ex1min
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
c locate       : subroutine to find value in ordered table
c Ebottom      : bottom of outgoing energy bin
c xso00,....   : cross section on excitation energy grid
c nenbeg,nenend: help variables
c dE           : emission bin
c dEhalf       : half of emission bin
c
c After the bottom (emin) and top (emax) of the excitation energy
c bin have been determined. Then the begin and end points for the
c spectrum are located.
c
            call locate(Ebottom,ebegin(type),eend(type),emin,nenbeg)
            call locate(Ebottom,ebegin(type),eend(type),emax,nenend)
            nenbeg=max(nenbeg,1)
            if (nenend.lt.nenbeg) goto 130
            dE=emax-emin
            dEhalf=0.5*dE
c
c For the interpolation, we first determine the feeding from adjacent
c mother bins.
c
            xso00=mcontrib(type,nex,nexout)
            xsom0=mcontrib(type,nex-1,nexout)
            xsop0=mcontrib(type,nex+1,nexout)
c
c Decay from continuum to continuum.
c We need to interpolate between adjacent residual bins. They are also
c determined. For the lowest emission energies we ensure that the limit
c is zero for zero outgoing energy.
c
            if (nexout.gt.NL) then
              xso0m=mcontrib(type,nex,nexout-1)
              xsopm=mcontrib(type,nex+1,nexout-1)
              xsomp=mcontrib(type,nex-1,nexout+1)
              xso0p=mcontrib(type,nex,nexout+1)
              do 140 nen=nenbeg,nenend
c
c A. Emission energy lower than the average emission energy.
c    Interpolate with lower bins.
c
c edist        : help variable
c term,term1...: help variables
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
  140         continue
            else
c
c Decay from continuum to discrete.
c We only need to interpolate between adjacent mother bins.
c
              do 150 nen=nenbeg,nenend
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
  150         continue
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
            do 160 nen=nenbeg,nenend
              sumweight=sumweight+weight(nen)
  160       continue
            if (sumweight.eq.0.) goto 130
c
c The contribution for mother bin --> residual bin is added to the
c spectrum.
c
c fracmpe   : multiple pre-equilibrium fraction
c mpecontrib: contribution to multiple pre-equilibrium emission spectrum
c factor    : help variable
c emisfac   : help variable
c compspect : compound part of spectrum
c preeqspect: multiple pre-equilibrium part of spectrum
c flagrecoil: flag for calculation of recoils
c comprecoil: subroutine for recoils from compound decay
c
c Note that below 20 MeV, generally fracmpe=0 (no multiple
c pre-equilibrium emission).
c
            fracmpe=mpecontrib(type,nex,nexout)/xso00
            do 170 nen=nenbeg,nenend
              factor=xso00/sumweight
              emisfac=factor*weight(nen)/deltaE(nen)
              compspect(nen)=(1.-fracmpe)*emisfac
              xsgrid(nen)=xsgrid(nen)+compspect(nen)
              preeqspect(nen)=fracmpe*emisfac
              xsmpegrid(nen)=xsmpegrid(nen)+preeqspect(nen)
  170       continue
            if (flagrecoil)
     +        call comprecoil(Zcomp,Ncomp,nex,type,nexout,nenbeg,nenend)
  130     continue
c
c Create final spectra
c
c xsmpeemis: multiple-preequilibrium emission spectrum from
c            compound nucleus
c
          do 180 nen=ebegin(type),eend(type)
            xsemis(type,nen)=xsemis(type,nen)+xsgrid(nen)
            xsmpeemis(type,nen)=xsmpeemis(type,nen)+xsmpegrid(nen)
  180     continue
  110   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
