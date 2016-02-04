      subroutine exciton2
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 21, 2012
c | Task  : Two-component exciton model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,p,ppi,hpi,pnu,hnu,h,n,type,nen,parity,J
      real    factor,xs
c
c ******************** Never-come-back solution ************************
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c Ecomp    : total energy of composite system
c Etotal   : total energy of compound system (target + projectile)
c exchange2: subroutine for calculation of two-component exchange terms
c p        : particle number
c maxpar   : maximal particle number
c ppi      : proton particle number
c ppi0     : initial proton number
c hpi      : proton hole number
c pnu      : neutron particle number
c pnu0     : initial neutron number
c hnu      : neutron hole number
c lifetime2: subroutine for calculation of lifetime of two-component
c            exciton state
c
c For each exciton number, we subsequently calculate the emission rates
c and exchange terms (both in subroutine exchange2) and the lifetime of
c the exciton state.
c
      Zcomp=0
      Ncomp=0
      Ecomp=Etotal
      call exchange2(Zcomp,Ncomp)
      do 10 p=p0,maxpar
        do 20 ppi=ppi0,maxpar
          hpi=ppi-ppi0
          do 20 pnu=pnu0,maxpar
            hnu=pnu-pnu0
            if (ppi+pnu.eq.p) call lifetime2(ppi,hpi,pnu,hnu)
   20   continue
   10 continue
c
c ************* Calculation of pre-equilibrium cross section ***********
c
c h         : hole number
c n         : exciton number
c factor,xs : help variables
c xsflux    : cross section flux
c Spre      : time-integrated strength of two-component exciton state
c parskip   : logical to skip outgoing particle
c preeqmode : designator for pre-equilibrium model
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c wemission2: two-component emission rate
c xsstep    : preequilibrium cross section per particle type, stage
c             and outgoing energy
c xsstep2   : two-component preequilibrium cross section
c xspreeq   : preequilibrium cross section per particle type and
c             outgoing energy
c
c The lifetimes and emission rates are processed into pre-equilibrium
c cross sections.
c
      do 110 ppi=ppi0,maxpar
        hpi=ppi-ppi0
        do 120 pnu=pnu0,maxpar
          hnu=pnu-pnu0
          p=ppi+pnu
          if (p.gt.maxpar) goto 110
          h=hpi+hnu
          n=p+h
          factor=xsflux*Spre(ppi,hpi,pnu,hnu)
          do 130 type=0,6
            if (parskip(type)) goto 130
            if (preeqmode.eq.4.and.(type.eq.1.or.type.eq.2)) goto 130
            do 140 nen=ebegin(type),eend(type)
              xs=factor*wemission2(type,ppi,hpi,pnu,hnu,nen)
              xsstep(type,p,nen)=xsstep(type,p,nen)+xs
              xsstep2(type,ppi,pnu,nen)=xs
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
                do 150 parity=-1,1,2
                  do 150 J=0,maxJph
                    xspreeqJP(type,nen,J,parity)=
     +                xspreeqJP(type,nen,J,parity)+
     +                xs*0.5*(2.*J+1.)*RnJ(n,J)/RnJsum(n)
  150           continue
              endif
  140       continue
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
