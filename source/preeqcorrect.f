      subroutine preeqcorrect
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 10, 2010
c | Task  : Correct pre-equilibrium cross sections for direct effects
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,NL,i,nen1,nen2,nen,p,ppi,pnu,parity,J
      real    esd,esd2,esd1,xs1,Elast,Rboundary
c
c ******************** Pre-equilibrium to direct ***********************
c
c 1. Correction of pre-equilibrium cross sections for direct discrete
c    cross sections. If the cross sections for discrete states have NOT
c    been calculated  by a direct reaction model, we collapse the
c    continuum pre-equilibrium cross sections in the high energy region
c    on the associated discrete states.
c
c parskip    : logical to skip outgoing particle
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c Ltarget    : excited level of target
c Nlast,NL   : last discrete level
c parZ       : charge number of particle
c parN       : neutron number of particle
c k0         : index for incident particle
c xsdirdisc  : direct cross section for discrete state
c eoutdis,esd: outgoing energy of discrete state reaction
c esd1,esd2  : help variables
c
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (ebegin(type).ge.eend(type)) goto 10
        if (Ltarget.ne.0) goto 100
        NL=Nlast(parZ(type),parN(type),0)
        do 20 i=0,NL
          if (xsdirdisc(type,i).ne.0.) goto 20
          esd=eoutdis(type,i)
          if (esd.lt.0.) goto 100
          if (i.eq.0.or.(type.eq.k0.and.i.eq.1)) then
            esd2=eoutdis(type,0)
          else
            esd2=0.5*(esd+eoutdis(type,i-1))
          endif
          if (i.eq.NL) then
            esd1=eoutdis(type,NL)
          else
            if (eoutdis(type,i+1).gt.0.) then
              esd1=0.5*(esd+eoutdis(type,i+1))
            else
              esd1=0.
            endif
          endif
c
c Find the part of the continuum spectrum that corresponds with the
c discrete states.
c
c xseps         : limit for cross sections
c locate        : subroutine to find value in ordered table
c egrid         : outgoing energy grid
c nendisc       : last discrete bin
c xs1           : help variable
c xspreeq       : preequilibrium cross section per particle type and
c                 outgoing energy
c xspreeqdisc   : preequilibrium cross section for discrete state
c xspreeqdisctot: preequilibrium cross section summed over discrete
c                 states
c xspreeqdiscsum: total preequilibrium cross section for discrete states
c dorigin       : origin of direct cross section (Direct or Preeq)
c
          if (type.eq.k0) then
            xs1=xseps
          else
            call locate(egrid,nendisc(type),eend(type),esd1,nen1)
            call locate(egrid,nendisc(type),eend(type),esd2,nen2)
            xs1=0.5*(xspreeq(type,nen1)+xspreeq(type,nen2))*(esd2-esd1)
          endif
          xspreeqdisc(type,i)=xs1
          xspreeqdisctot(type)=xspreeqdisctot(type)+xs1
          xspreeqdiscsum=xspreeqdiscsum+xs1
          dorigin(type,i)='Preeq '
   20   continue
c
c 2. Set pre-equilibrium cross section in discrete energy region to
c    zero.
c
c Elast    : help variable
c Rboundary: correction factor
c Etop     : top of outgoing energy bin
c deltaE   : energy bin around outgoing energies
c xspreeqps: preequilibrium cross section per particle type and
c            outgoing energy for pickup and stripping
c xspreeqki: preequilibrium cross section per particle type and
c            outgoing energy for knockout and inelastic
c xspreeqbu: preequilibrium cross section per particle type and
c            outgoing energy for breakup
c maxpar   : maximal particle number
c xsstep   : preequilibrium cross section per particle type, stage
c            and outgoing energy
c flag2comp: flag for two-component pre-equilibrium model
c ppi0     : initial proton number
c pnu0     : initial neutron number
c xsstep2  : two-component preequilibrium cross section
c parity   : parity
c maxJph   : maximal spin for particle-hole states
c xspreeqJP: preequilibrium cross section per particle type,
c            outgoing energy, spin and parity
c
c The first continuum outgoing energy bin is only partially depleted by
c the last discrete level. This is corrected using Rboundary.
c
        nen=nendisc(type)
        Elast=eoutdis(type,NL)
        Rboundary=(Etop(nen)-Elast)/deltaE(nen)
        if (abs(Rboundary).gt.1.) Rboundary=0.
        xspreeq(type,nen)=xspreeq(type,nen)*(1.-Rboundary)
        xspreeqps(type,nen)=xspreeqps(type,nen)*(1.-Rboundary)
        xspreeqki(type,nen)=xspreeqki(type,nen)*(1.-Rboundary)
        xspreeqbu(type,nen)=xspreeqbu(type,nen)*(1.-Rboundary)
        do 110 p=1,maxpar
          xsstep(type,p,nen)=xsstep(type,p,nen)*(1.-Rboundary)
  110   continue
        if (flag2comp) then
          do 120 ppi=ppi0,maxpar
            do 120 pnu=pnu0,maxpar
              xsstep2(type,ppi,pnu,nen)=xsstep2(type,ppi,pnu,nen)*
     +          (1.-Rboundary)
  120     continue
        endif
        do 130 parity=-1,1,2
          do 130 J=0,maxJph
            xspreeqJP(type,nen,J,parity)=xspreeqJP(type,nen,J,parity)*
     +        (1.-Rboundary)
  130   continue
c
c The pre-equilibrium spectrum for energies corresponding to discrete
c transitions is set to zero.
c
  100   do 140 nen=nendisc(type)+1,eend(type)
          xspreeq(type,nen)=0.
          xspreeqps(type,nen)=0.
          xspreeqki(type,nen)=0.
          xspreeqbu(type,nen)=0.
          do 150 p=1,maxpar
            xsstep(type,p,nen)=0.
  150     continue
          if (flag2comp) then
            do 160 ppi=ppi0,maxpar
              do 160 pnu=pnu0,maxpar
                xsstep2(type,ppi,pnu,nen)=0.
  160       continue
          endif
          do 170 parity=-1,1,2
            do 170 J=0,maxJph
              xspreeqJP(type,nen,J,parity)=0.
  170     continue
  140   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
