      subroutine msdplusmsc
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 13, 2004
c | Task  : Total quantum-mechanical pre-equilibrium cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,ns,parity,J,iang
c
c *********** Angle-integrated pre-equilibrium cross sections **********
c
c parskip: logical to skip outgoing particle
c ebegin : first energy point of energy grid
c eend   : last energy point of energy grid
c xspreeq: preequilibrium cross section per particle type and
c          outgoing energy
c msdtot : multi-step direct cross section summed over steps
c maxmsd : number of MSD steps
c xsstep : preequilibrium cross section per particle type, stage
c          and outgoing energy
c msdstep: continuum n-step direct cross section
c
      do 10 type=1,2
        if (parskip(type)) goto 10
        do 20 nen=ebegin(type),eend(type)
          xspreeq(type,nen)=msdtot(type,nen)
   20   continue
        do 30 ns=1,maxmsd
          do 40 nen=ebegin(type),eend(type)
            xsstep(type,ns,nen)=msdstep(type,ns,nen)
c
c Spin distribution the same as in exciton model. This will be
c changed in the true MSD model.
c
c Create J-dependent pre-equilibrium cross sections using the spin
c distribution. The result is normalized with the sum of RnJ over J.
c
c parity   : parity
c maxJph   : maximal spin for particle-hole states
c xspreeqJP: preequilibrium cross section per particle type,
c            outgoing energy, spin and parity
c RnJ      : spin distribution for particle-hole states
c RnJsum   : (2J+1)*sum over spin distributions
c
            do 50 parity=-1,1,2
              do 50 J=0,maxJph
                xspreeqJP(type,nen,J,parity)=
     +            xspreeqJP(type,nen,J,parity)+
     +            xsstep(type,ns,nen)*0.5*(2.*J+1.)*RnJ(ns*2+1,J)/
     +            RnJsum(ns*2+1)
   50       continue
   40     continue
   30   continue
   10 continue
c
c **************** Pre-equilibrium angular distributions ***************
c
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c xspreeqad : preequilibrium angular distribution per particle type
c             and outgoing energy
c msdtotad  : multi-step direct angular distribution summed over steps
c
      if (.not.flagddx) return
      do 110 type=1,2
        if (parskip(type)) goto 110
        do 120 iang=0,nanglecont
          do 130 nen=ebegin(type),eend(type)
            xspreeqad(type,nen,iang)=msdtotad(type,nen,iang)
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
