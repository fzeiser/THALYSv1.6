      subroutine urr
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Gilles Noguere
c | Date  : August 15, 2013
c | Task  : Unresolved resonance range parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,J,type,i,istop
      real    sum,sumgJ,gJ,err,Rsig0,Rsig1,Rsig
c
c Channel widths and neutron strength function
c
c xselasinc   : total elastic cross section (neutrons only) for incident
c               channel
c fourpi      : 4.*pi
c dtheory     : subroutine for theoretical calculation of average
c               neutron spacings
c lurr        : maximal orbital angular momentum for URR calculation
c lminU,lmaxU : minimal and maximal orbital angular momentum
c Einc        : incident energy in MeV
c JminU,JmaxU : minimal and maximal total angular momentum
c strengthfunc: function for (l,j) neutron strength function for URR
c strengthlj  : (l,j) neutron strength function
c strengthl   : l neutron strength function
c targetspin2 : 2 * spin of target
c Turrljinc   : incident channel (l,j) transmission coefficient for
c               URR calculation
c urrwidth    : channel width in URR
c Dlj         : mean resonance spacing per J,l value
c Turrlj      : (l,j) transmission coefficients
c
      call dtheory(0,0,Einc)
      do 10 l=lminU,min(lmaxU,lurr)
        sum=0.
        sumgJ=0.
        do 20 J=JminU(l),JmaxU(l)
          if (nulj(0,l,J).gt.0)
     +      Turrljinc(l,J)=Turrljinc(l,J)/real(nulj(0,l,J))
          call strengthfunc(Turrljinc(l,J),l,J)
          gJ=(2*J+1.)/(2.*(targetspin2+1.))
          sumgJ=sumgJ+gJ
          sum=sum+gJ*strengthlj(l,J)
          do 30 type=-1,6
            urrwidth(type,l,J)=Dlj(l,J)*Turrlj(type,l,J)/twopi
   30     continue
          urrwidth(3,l,J)=Dlj(l,J)*strengthlj(l,J)*1.e-4
   20   continue
        if (sumgJ.gt.0.) strengthl(l)=sum/sumgJ
   10 continue
c
c Cross section after correction with NJOY formalism (lmaxU=2)
c
c flagurrnjoy : normalization of URR parameters with NJOY method
c Rprime(0,U) : potential scattering radius
c
      if (RprimeU.eq.0.) RprimeU=Rprime
      if (flagurrnjoy) then
        Rprime0=0.1*RprimeU
        err=0.005
        Rsig=0.
        Rsig0=1.
        Rsig1=1.
        istop=0
        do while (Rsig.gt.(1.+err).or.Rsig.lt.(1.-err))
          istop=istop+1
          do 110 i=1,4
            xsurrN(i)=0.
            xsurrT(i)=0.
  110     continue
          do 120 l=lminU,min(2,lurr)
            do 130 J=JminU(l),JmaxU(l)
              urrwidth(1,l,J)=urrwidth(1,l,J)/Rsig0
              urrwidth(-1,l,J)=urrwidth(-1,l,J)*Rsig1
              call csunr2(Rprime0,l,J)
              xsurrN(3)=xsurrN(3)+sigurrf(l,J)*1000.
              xsurrN(4)=xsurrN(4)+sigurrc(l,J)*1000.
              xsurrT(1)=xsurrT(1)+xsbinarylj(1,l,J)
              xsurrT(3)=xsurrT(3)+xsbinarylj(-1,l,J)
              xsurrT(4)=xsurrT(4)+xsbinarylj(0,l,J)
  130       continue
  120     continue
          if (xsurrN(4).gt.0.) then
            Rsig0=xsurrT(4)/xsurrN(4)
          else
            Rsig0=1.
          endif
          if (xsurrN(3).gt.0.) then
            Rsig1=xsurrT(3)/xsurrN(3)
          else
            Rsig1=1.
          endif
          if (xsurrT(4).gt.xsurrT(3)) then
            if (xsurrT(1).gt.0.) then
              Rsig=Rsig0
              Rsig1=1.
            else
              Rsig=1.
              Rsig1=1.
            endif
          else
            if (Rsig1.gt.(1.+err).or.Rsig1.lt.(1.-err)) then
              Rsig=Rsig1
              Rsig0=1.
            else
              if (xsurrT(1).gt.0.) then
                Rsig=Rsig0
                Rsig1=1.
              else
                Rsig=1.
                Rsig1=1.
              endif
            endif
          endif
          if (istop.gt.100) Rsig=1.
        enddo
      endif
c
c Final cross sections
c
      do 140 i=1,4
        xsurrN(i)=0.
        xsurrT(i)=0.
  140 continue
      do 150 l=lminU,min(2,lurr)
        do 160 J=JminU(l),JmaxU(l)
          if (flagurrnjoy) then
            call csunr2(Rprime0,l,J)
            xsurrN(1)=xsurrN(1)+xsbinarylj(1,l,J)
            xsurrN(2)=xsurrN(2)+sigurrs(l,J)*1000.
            xsurrN(3)=xsurrN(3)+sigurrf(l,J)*1000.
            xsurrN(4)=xsurrN(4)+sigurrc(l,J)*1000.
          endif
          xsurrT(1)=xsurrT(1)+xsbinarylj(1,l,J)
          xsurrT(3)=xsurrT(3)+xsbinarylj(-1,l,J)
          xsurrT(4)=xsurrT(4)+xsbinarylj(0,l,J)
  160   continue
        if (flagurrnjoy) xsurrN(2)=xsurrN(2)+spot(l)*1000.
  150 continue
      if (flagurrnjoy) xsurrN(1)=xsurrN(1)+xsurrN(2)+xsurrN(3)+xsurrN(4)
      xsurrT(2)=xselastot
      xsurrT(1)=xsurrT(1)+xsurrT(2)+xsurrT(3)+xsurrT(4)
c
c Output of URR
c
c urrout: subroutine for output of unresolved resonance parameters
c         in separate files
c
      call urrout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
