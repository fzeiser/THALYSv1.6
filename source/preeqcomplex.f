      subroutine preeqcomplex
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 21, 2012
c | Task  : Pre-equilibrium complex particle emission
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,parity,J
      real    pecompsum,factor,xspecomp
c
c ********************** Various components ****************************
c
c stripping: subroutine for contribution of stripping and pickup
c            reactions
c knockout : subroutine for contribution of knockout reactions
c k0       : index of incident particle
c breakup  : subroutine for contribution of breakup reactions
c
      call stripping
      call knockout
      if (k0.gt.2) call breakup
c
c *************************** Corrections ******************************
c
c Prevent complex particle pre-equilibrium to exceed the reaction cross
c section.
c
c pecompsum: help variable
c ebegin   : first energy point of energy grid
c eend     : last energy point of energy grid
c xspecomp : pre-equilibrium complex particle cross section
c xspreeqps: preequilibrium cross section per particle type and
c            outgoing energy for pickup and stripping
c xspreeqki: preequilibrium cross section per particle type and
c            outgoing energy for knockout and inelastic
c xspreeqbu: preequilibrium cross section per particle type and
c            outgoing energy for breakup
c deltaE   : energy bin around outgoing energies
c factor   : help variable
c xsflux   : cross section flux
c
      pecompsum=0.
      do 10 type=1,6
        do 20 nen=ebegin(type),eend(type)
        xspecomp=xspreeqps(type,nen)+xspreeqki(type,nen)+
     +    xspreeqbu(type,nen)
        pecompsum=pecompsum+xspecomp*deltaE(nen)
   20   continue
   10 continue
      if (pecompsum.gt.xsflux) then
        factor=xsflux/pecompsum
        do 30 type=1,6
          do 40 nen=ebegin(type),eend(type)
            xspreeqps(type,nen)=xspreeqps(type,nen)*factor
            xspreeqki(type,nen)=xspreeqki(type,nen)*factor
            xspreeqbu(type,nen)=xspreeqbu(type,nen)*factor
   40     continue
   30   continue
      endif
c
c If the pre-equilibrium spin distribution is chosen, we assume that the
c spin distribution for pickup, stripping and knockout is the same as in
c the exciton model.
c
c pespinmodel: model for pre-equilibrium spin distribution or compound
c              spin distribution for pre-equilibrium cross section
c xspreeq    : preequilibrium cross section per particle type and
c              outgoing energy
c parity     : parity
c maxJph     : maximal spin for particle-hole states
c xspreeqJP  : preequilibrium cross section per particle type,
c              outgoing energy, spin and parity
c
      do 110 type=1,6
        do 120 nen=ebegin(type),eend(type)
          xspecomp=xspreeqps(type,nen)+xspreeqki(type,nen)+
     +      xspreeqbu(type,nen)
          if (pespinmodel.eq.3.and.xspreeq(type,nen).ne.0.) then
            do 130 parity=-1,1,2
              do 130 J=0,maxJph
                factor=xspreeqJP(type,nen,J,parity)/xspreeq(type,nen)
                xspreeqJP(type,nen,J,parity)=
     +            xspreeqJP(type,nen,J,parity)+factor*xspecomp
  130       continue
          endif
          xspreeq(type,nen)=xspreeq(type,nen)+xspecomp
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
