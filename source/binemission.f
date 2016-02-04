      subroutine binemission
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 30, 2011
c | Task  : Compound emission cross sections for binary reaction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Zix,Nix,NL,nexout,nen,na,nb,nenbeg,nenend
      real    SS,dEx,Exout,Eo(0:numex),xsMeV(0:numex+1),Eout,dE,xs,Ea,
     +        Eb,xsa,xsb,emissum,frac,emin,emax,fracbot,fractop
c
c ******************** Calculation of binary spectra *******************
c
c parskip    : logical to skip outgoing particle
c binnorm    : normalization factor for binary spectra
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c Nlast,NL   : last discrete level
c SS,S       : separation energy
c maxex      : maximum excitation energy bin for compound nucleus
c deltaEx,dEx: excitation energy bin for population arrays
c Ex,Exout   : excitation energy
c Eo         : outgoing energy grid based on excitation energy
c Exinc      : excitation energy of entrance bin
c xsMeV      : cross section on outgoing energy grid
c contrib    : contribution to emission spectrum
c
c The decay from the primary compound nucleus to residual nuclei is
c converted from the excitation energy grid to emission energies.
c Decay from the primary compound nucleus to discrete states is
c already taken into account and the emission energies for these
c transitions are exactly known. The emission energies Eo correspond
c exactly with the excitation energy grid points. The contribution
c per MeV, xsMeV, is constructed.
c
      do 10 type=0,6
        if (parskip(type)) goto 10
        binnorm(type)=0.
        if (ebegin(type).ge.eend(type)) goto 10
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        NL=Nlast(Zix,Nix,0)
        SS=S(0,0,type)
        if (maxex(Zix,Nix).le.NL) goto 10
        do 20 nexout=NL+1,maxex(Zix,Nix)
          dEx=deltaEx(Zix,Nix,nexout)
          Exout=Ex(Zix,Nix,nexout)
          Eo(nexout)=Exinc-SS-Exout
          xsMeV(nexout)=contrib(type,nexout)/dEx
   20   continue
c
c To avoid unphysical interpolations, the contribution for the
c last discrete level is temporarily set to that of the first
c continuum bin.
c
c Eout : outgoing energy
c egrid: outgoing energy grid
c
        Eo(NL)=Exinc-SS-Ex(Zix,Nix,NL)
        xsMeV(NL)=xsMeV(NL+1)
        xsMeV(maxex(Zix,Nix)+1)=0.
        do 30 nen=ebegin(type),eend(type)
          Eout=egrid(nen)
c
c Transitions to discrete states are skipped. For the lowest emission
c energies we ensure that the limit is zero for zero outgoing energy.
c
c dE,xs: help variables
c
          if (Eout.ge.Eo(NL)) goto 30
          if (Eout.le.Eo(maxex(Zix,Nix))) then
            nexout=maxex(Zix,Nix)
            dE=Eout/Eo(nexout)
            xs=xsMeV(nexout)*dE
            goto 40
          endif
c
c To obtain the emission energies belonging to the excitation energy
c grid points, we use first order interpolation.
c
c locate : subroutine to find value in ordered table
c na,nb  : help variables
c Ea,Eb  : help variables
c xsa,xsb: help variables
c pol1   : subroutine for interpolation of first order
c speceps: limit for cross section spectra
c xsemis : cross section for emission from compound nucleus
c
          call locate(Eo,NL,maxex(Zix,Nix),Eout,nexout)
          na=nexout
          nb=nexout+1
          Ea=Eo(na)
          Eb=Eo(nb)
          xsa=xsMeV(na)
          xsb=xsMeV(nb)
          call pol1(Ea,Eb,xsa,xsb,Eout,xs)
          if (xs.lt.speceps) xs=0.
   40     xsemis(type,nen)=xs
   30   continue
c
c Normalisation of binary spectra to binary cross sections
c
c eoutdis        : outgoing energy of discrete state reaction
c nendisc        : last discrete bin
c emissum        : integrated binary emission spectrum
c deltaE         : energy bin around outgoing energies
c frac           : help variable
c xscompcont     : compound cross section for continuum
c flagchannels   : flag for exclusive channels calculation
c emax           : maximal emission energy within bin decay
c emin           : minimal emission energy
c Ebottom        : bottom of outgoing energy bin
c nenbeg,nenend  : help variables
c binemis        : emission spectra from initial compound nucleus
c Etop           : top of outgoing energy bin
c fraclow,fracbot: help variables
c xspreeq        : preequilibrium cross section per particle type and
c                  outgoing energy
c xsgr           : smoothed giant resonance cross section
c
c Due to interpolation errors, there is always a small difference
c between the binary continuum cross section and the integral of the
c spectrum over the continuum. The spectrum is accordingly normalized.
c The fraction of the first continuum bin is taken into account.
c
        emissum=0.
        do 110 nen=ebegin(type),nendisc(type)
          emissum=emissum+xsemis(type,nen)*deltaE(nen)
  110   continue
        if (eoutdis(type,NL).gt.0.) then
          frac=Etop(nendisc(type))-eoutdis(type,NL)
          emissum=emissum-xsemis(type,nendisc(type))*frac
        endif
        if (emissum.ne.0.) binnorm(type)=xscompcont(type)/emissum
        do 120 nen=ebegin(type),eend(type)
          xsemis(type,nen)=binnorm(type)*xsemis(type,nen)
  120   continue
c
c Exclusive channels.
c The formula for exclusive reaction channels shows that we
c need to store the emission spectrum per excitation energy bin for
c each outgoing particle type. This is done in array binemis.
c First the bottom (emin) and top (emax) of the excitation energy
c bin are determined. Then the begin and end points for the spectrum
c are located. The preequilibrium contribution is also added.
c
        if (flagchannels) then
          do 130 nexout=NL+1,maxex(Zix,Nix)
            dEx=deltaEx(Zix,Nix,nexout)
            emin=Eo(nexout)-0.5*dEx
            emax=Eo(nexout)+0.5*dEx
            call locate(Ebottom,ebegin(type),eend(type),emin,nenbeg)
            call locate(Ebottom,ebegin(type),eend(type),emax,nenend)
            nenbeg=max(nenbeg,1)
            if (nenend.lt.nenbeg) goto 130
            do 140 nen=nenbeg,nenend
              if (egrid(nen).ge.Eo(NL)) goto 140
              binemis(type,nexout,nen)=xsemis(type,nen)+
     +          xspreeq(type,nen)+xsgr(type,nen)
  140       continue
c
c The end points of the bin are taken into account and this is
c corrected in binemis. There are also cases where the excitation
c energy bin falls completely in the emission energy bin
c (nenbeg=nenend).
c
            fracbot=(emin-Ebottom(nenbeg))/deltaE(nenbeg)
            fractop=(Etop(nenend)-emax)/deltaE(nenend)
            if (nexout.eq.NL+1) fractop=0.
            if (nenend.eq.nenbeg) then
              binemis(type,nexout,nenbeg)=binemis(type,nexout,nenbeg)*
     +          (1.-fracbot-fractop)
            else
              binemis(type,nexout,nenbeg)=binemis(type,nexout,nenbeg)*
     +          (1.-fracbot)
              binemis(type,nexout,nenend)=binemis(type,nexout,nenend)*
     +          (1.-fractop)
            endif
  130     continue
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
