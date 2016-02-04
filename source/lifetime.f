      subroutine lifetime(Zcomp,Ncomp,p,h)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 9, 2004
c | Task  : Calculation of lifetime of exciton state
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,p,h
      real    tplus,lambdaplus,term1,term2
c
c ******************** Never-come-back solution ************************
c
c Zcomp           : charge number index for compound nucleus
c Ncomp           : neutron number index for compound nucleus
c p               : particle number
c h               : hole number
c p0              : initial particle number
c depletion       : depletion factor at each stage
c lambdaplus,tplus: transition rate for n --> n+2
c term1,term2     : help variables
c wemistot        : total emission rate per exciton number
c tauexc          : mean lifetime
c
c The lifetime of the exciton state is calculated, see the manual.
c
      if (p.eq.p0) then
        depletion(p,h)=1.
      else
        tplus=lambdaplus(Zcomp,Ncomp,p-1,h-1)
        term1=tplus+wemistot(p-1,h-1)
        if (term1.ne.0.) then
          depletion(p,h)=depletion(p-1,h-1)*tplus/term1
        else
          depletion(p,h)=0.
        endif
      endif
      term2=lambdaplus(Zcomp,Ncomp,p,h)+wemistot(p,h)
      if (term2.ne.0.) then
        tauexc(p,h)=depletion(p,h)/term2
      else
        tauexc(p,h)=0.
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
