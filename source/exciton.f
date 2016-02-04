      subroutine exciton
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 21, 2012
c | Task  : Exciton model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,p,h,n,type,nen,parity,J
      real    factor,xs
c
c ******************** Never-come-back solution ************************
c
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c Ecomp       : total energy of composite system
c Etotal      : total energy of compound system (target + projectile)
c p           : particle number
c p0          : initial particle number
c maxpar      : maximal particle number
c h           : hole number
c emissionrate: subroutine for emission rate
c lifetime    : subroutine for calculation of lifetime of exciton state
c
c For each exciton number, we subsequently calculate the emission
c rates and the lifetime of the exciton state according to the
c never-come-back approximation.
c
      Zcomp=0
      Ncomp=0
      Ecomp=Etotal
      do 10 p=p0,maxpar
        h=p-p0
        call emissionrate(Zcomp,Ncomp,p,h)
        call lifetime(Zcomp,Ncomp,p,h)
   10 continue
c
c ************* Calculation of pre-equilibrium cross section ***********
c
c n         : exciton number
c factor,xs: help variables
c xsflux   : cross section flux
c tauexc   : mean lifetime
c parskip  : logical to skip outgoing particle
c preeqmode: designator for pre-equilibrium model
c ebegin   : first energy point of energy grid
c eend     : last energy point of energy grid
c wemission: emission rate
c xsstep   : preequilibrium cross section per particle type, stage
c            and outgoing energy
c xspreeq  : preequilibrium cross section per particle type and
c            outgoing energy
c
c The lifetimes and emission rates are processed into pre-equilibrium
c cross sections.
c
      do 110 p=p0,maxpar
        h=p-p0
        n=p+h
        factor=xsflux*tauexc(p,h)
        do 120 type=0,6
          if (parskip(type)) goto 120
c
c Neutron and proton emission can also be calculated with MSD/MSC,
c if requested.
c
          if (preeqmode.eq.4.and.(type.eq.1.or.type.eq.2)) goto 120
          do 130 nen=ebegin(type),eend(type)
            xs=factor*wemission(type,p,h,nen)
            xsstep(type,p,nen)=xs
            xspreeq(type,nen)=xspreeq(type,nen)+xs
c
c Create J-dependent pre-equilibrium cross sections using the spin
c distribution. The result is normalized with the sum of RnJ over J.
c As an alternative, the Hauser-Feshbach spin distribution is adopted.
c
c pespinmodel: model for pre-equilibrium spin distribution or compound
c              spin distribution for pre-equilibrium cross section
c parity     : parity
c J          : total angular momentum
c maxJph     : maximal spin for particle-hole states
c xspreeqJP  : preequilibrium cross section per particle type,
c              outgoing energy, spin and parity
c RnJ        : spin distribution for particle-hole states
c RnJsum     : (2J+1)*sum over spin distributions
c
            if (pespinmodel.eq.3) then
              do 140 parity=-1,1,2
                do 140 J=0,maxJph
                  xspreeqJP(type,nen,J,parity)=
     +              xspreeqJP(type,nen,J,parity)+
     +              xs*0.5*(2.*J+1.)*RnJ(n,J)/RnJsum(n)
  140         continue
            endif
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
