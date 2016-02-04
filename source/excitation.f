      subroutine excitation
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : December 30, 2011
c | Task  : Excitation energy population
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,A,odd,NL,nex0,nex,nen,parity,J
      real    Eex,dEx,ald,ignatyuk,Ea,Eb,Pa,Pb,Rspin,Probex,Prob,spindis
c
c ******************** Fill energy bins with population ****************
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c AA,A     : mass number of residual nucleus
c odd      : odd (1) or even (0) nucleus
c Nlast,NL : last discrete level
c Ex,Eex   : excitation energy
c dEx      : excitation energy bin
c edis     : energy of level
c deltaEx  : excitation energy bin for population arrays
c locate   : subroutine to find value in ordered table
c Exdist   : excitation energy of population distribution
c npopbins : number of excitation energy bins for population
c            distribution
c Ea,Eb,.. : help variables
c ignatyuk : function for energy dependent level density parameter a
c maxex    : maximum excitation energy bin for compound nucleus
c maxJ     : maximal J-value
c jdis     : spin of level
c parlev   : parity of level
c xspop    : population cross section
c npopJ    : number of spins for population distribution
c Pdistex  : population distribution, spin-independent
c Pdist    : population distribution per spin and parity
c spindis  : spin distribution
c pardis   : parity distribution
c xspopex  : population cross section summed over spin and parity
c xspopnuc : population cross section per nucleus
c feedexcl : feeding terms from compound excitation energy bin to
c            residual excitation energy bin
c xsinitpop: initial population cross section
c popexcl  : population cross section of bin just before decay
c
      Zcomp=0
      Ncomp=0
      A=AA(Zcomp,Ncomp,0)
      odd=mod(A,2)
      NL=Nlast(Zcomp,Ncomp,0)
      nex0=maxex(Zcomp,Ncomp)+1
      do 10 nex=1,maxex(Zcomp,Ncomp)
        Eex=Ex(Zcomp,Ncomp,nex)
        if (nex.le.NL) then
          dEx=0.5*(edis(Zcomp,Ncomp,min(nex+1,NL))-
     +      edis(Zcomp,Ncomp,nex-1))
        else
          dEx=deltaEx(Zcomp,Ncomp,nex)
        endif
        call locate(Exdist,0,npopbins,Eex,nen)
        nen=max(0,nen)
        Ea=Exdist(nen)
        Eb=Exdist(nen+1)
        if (npopJ.eq.0) then
          Pa=Pdistex(nen)
          Pb=Pdistex(nen+1)
          call pol1(Ea,Eb,Pa,Pb,Eex,Probex)
          ald=ignatyuk(Zcomp,Ncomp,Eex,0)
        endif
        do 20 parity=-1,1,2
          do 30 J=0,maxJ(Zcomp,Ncomp,nex)
            Rspin=real(J)+0.5*odd
            if (nex.le.NL.and.(jdis(Zcomp,Ncomp,nex).ne.J.or.
     +        parlev(Zcomp,Ncomp,nex).ne.parity)) goto 30
            if (npopJ.eq.0) then
              Prob=Probex*spindis(Zcomp,Ncomp,Eex,ald,Rspin,0)*
     +          pardis
            else
              Pa=Pdist(nen,J,parity)
              Pb=Pdist(nen+1,J,parity)
              call pol1(Ea,Eb,Pa,Pb,Eex,Prob)
            endif
            xspop(Zcomp,Ncomp,nex,J,parity)=Prob*dEx
            xspopex(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)+
     +        xspop(Zcomp,Ncomp,nex,J,parity)
   30     continue
   20   continue
        xspopnuc(Zcomp,Ncomp)=xspopnuc(Zcomp,Ncomp)+
     +    xspopex(Zcomp,Ncomp,nex)
        feedexcl(Zcomp,Ncomp,0,nex0,nex)=xspopex(Zcomp,Ncomp,nex)
   10 continue
      xsinitpop=xspopnuc(Zcomp,Ncomp)
      popexcl(Zcomp,Ncomp,nex0)=xsinitpop
      return
      end
