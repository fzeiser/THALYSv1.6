      subroutine dtheory(Zix,Nix,E)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 1, 2012
c | Task  : Theoretical calculation of average neutron spacings
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,l,J,Nres,L0,tspin2,tpar,parspin2i,J2b,
     +                 J2e,parity,J2,jj2beg,jj2end,jj2,l2beg,l2end,l2
      real             E,tspin
      double precision rho(0:numl,0:numJ),rhosum,density
c
c ********************** Level density parameters **********************
c
c Zix      : charge number index for residual nucleus
c Nix      : neutron number index for residual nucleus
c E        : incident energy
c l        : orbital angular momentum
c Dl       : theoretical resonance spacing per l value
c Dlj      : mean resonance spacing per J,l value
c Nres     : maximal neutron number index for residual nucleus
c numN     : maximal number of neutrons away from the initial compound
c            nucleus
c levels   : subroutine for discrete levels
c tspin    : target spin
c jdis     : spin of level
c tspin2   : 2 * target spin
c tpar     : target parity
c parlev   : parity of level
c parspin2i: 2 * particle spin for incident channel
c parspin  : spin of incident particle
c k0       : index of incident particle
c J2b      : 2 * start of J summation
c J2e      : 2 * end of J summation
c lmaxinc  : maximal l-value for transmission coefficients for
c            incident channel
c parity   : parity
c pardif   : difference between target and compound nucleus parity
c jj2beg   : 2 * start of j summation
c jj2end   : 2 * end of j summation
c J        : J-value
c l2beg    : 2 * start of l summation
c l2end    : 2 * end of l summation
c sgn      : +1 for even argument, -1 for odd argument
c density  : level density
c S        : separation energy per particle
c ldmodel  : level density model
c rhosum   : help variable
c
      do 10 l=0,numl
        Dl(l)=0.
        do 10 J=0,numJ
          Dlj(l,J)=0.
          rho(l,J)=0.
   10 continue
      Nres=min(numN,Nix+1)
      call levels(Zix,Nres)
      if (Zix.eq.0.and.Nix.eq.0) then
        L0=Ltarget
      else
        L0=0
      endif
      tspin=jdis(Zix,Nres,L0)
      tspin2=int(2.*jdis(Zix,Nres,L0))
      tpar=parlev(Zix,Nres,L0)
      parspin2i=1
      J2b=mod(int(2.*(tspin+0.5)),2)
      J2e=int(2*(lmaxinc+0.5+tspin))
      J2e=min(J2e,2*numJ)
      do 110 parity=-1,1,2
        pardif=abs(tpar-parity)/2
        do 120 J2=J2b,J2e,2
          J=J2/2
          jj2beg=abs(J2-tspin2)
          jj2end=J2+tspin2
          do 130 jj2=jj2beg,jj2end,2
            l2beg=abs(jj2-parspin2i)
            l2end=jj2+parspin2i
            l2end=min(l2end,2*lmaxinc)
            do 140 l2=l2beg,l2end,2
              l=l2/2
              if (mod(l,2).ne.pardif) goto 140
              rho(l,J)=density(Zix,Nix,max(0.,S(Zix,Nix,1)+E),0.5*J2,
     +          parity,0,ldmodel(Zix,Nix))
              Dlj(l,J)=real(1.e6/rho(l,J))
  140       continue
  130     continue
  120   continue
  110 continue
      do 150 l=0,numl
        rhosum=0.
        do 160 J=0,numJ
          rhosum=rhosum+rho(l,J)
  160   continue
        if (rhosum.ge.1.e-10) Dl(l)=real(1.e6/rhosum)
  150 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
