      subroutine preeqinit
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 1, 2008
c | Task  : Initialization of general pre-equilibrium parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer n,k,J,p,h,ppi,hpi,pnu,hnu,i
      real    sigma2ph,denom,expo,gs,Epp,factor,gsp,gsn,Eppi,
     +        factorp,Epnu,factorn
c
c ********************** General initializations ***********************
c
c maxexc: maximal exciton number
c numexc: maximal exciton number (set in talys.cmb)
c maxpar: maximal particle number
c Efermi: depth of Fermi well
c
      maxexc=numexc
      maxpar=maxexc/2
      Efermi=38.
c
c ******************* Factorials and combinatorials ********************
c
c nfac : n!
c ncomb: n!/(n!(n-k)!)
c n    : exciton number
c
      nfac(0)=1.
      ncomb(0,0)=1.
      do 10 n=1,numexc
        nfac(n)=real(n)*nfac(n-1)
        ncomb(n,0)=1.
        do 10 k=1,n
          ncomb(n,k)=nfac(n)/(nfac(k)*nfac(n-k))
   10 continue
c
c ************* Spin distribution for particle-hole states *************
c
c maxJph   : maximal spin for particle-hole states
c RnJsum   : (2J+1)*sum over spin distributions
c sigma2ph : spin cutoff factor for particle-hole states
c Atarget  : mass number of target nucleus
c twothird : 2/3
c J        : spin
c denom    : help variable
c sqrttwopi: sqrt(2.*pi)
c expo     : help variable
c RnJ      : spin distribution for particle-hole states
c
c Default spin cutoff factor for p-h density: Gruppelaar,
c Group Meeting on level densities, Brookhaven (1983), p. 143.
c
      maxJph=30
      do 110 n=1,maxexc
        RnJsum(n)=0.
        sigma2ph=0.24*n*Atarget**twothird
        do 120 J=0,maxJph
          denom=2.*sqrttwopi*sigma2ph**1.5
          expo=-(J+0.5)**2/(2.*sigma2ph)
          RnJ(n,J)=(2.*J+1)/denom*exp(expo)
          RnJsum(n)=RnJsum(n)+(2.*J+1)*RnJ(n,J)
  120   continue
  110 continue
c
c ******************** Pauli correction factors ************************
c
c gs,g   : single-particle level density parameter
c p      : particle number
c numparx: maximal particle number (set in talys.cmb)
c h      : hole number
c Apauli : Pauli blocking correction energy
c Epp    : Pauli energy
c factor : help variable
c
c 1. One component
c
      gs=g(0,0)
      do 210 p=-1,numparx+1
        do 210 h=-1,numparx+1
          if (p.eq.-1.or.h.eq.-1) then
            Apauli(p,h)=0.
          else
            Epp=max(p,h)**2/gs
            factor=(p*p+h*h+p+h)/(4.*gs)
            Apauli(p,h)=Epp-factor
          endif
  210 continue
c
c 2. Two components (Kalbach, Phys. Rev. C33, 818 (1986)).
c
c gp,gsp         : single-particle proton level density parameter
c gn,gsn         : single-particle neutron level density parameter
c ppi,hpi        : proton particle and hole number
c pnu,hnu        : neutron particle and hole number
c Eppi,..        : Pauli energy
c factorp,factorn: help variables
c Apauli2        : two-component Pauli blocking correction factor
c
      gsp=gp(0,0)
      gsn=gn(0,0)
      do 220 ppi=-1,numparx+1
        do 220 hpi=-1,numparx+1
          Eppi=max(ppi,hpi)**2/gsp
          factorp=(ppi*ppi+hpi*hpi+ppi+hpi)/(4.*gsp)
          do 220 pnu=-1,numparx+1
            do 220 hnu=-1,numparx+1
              if (ppi.eq.-1.or.hpi.eq.-1.or.pnu.eq.-1.or.hnu.eq.-1) then
                Apauli2(ppi,hpi,pnu,hnu)=0.
              else
                Epnu=max(pnu,hnu)**2/gsn
                factorn=(pnu*pnu+hnu*hnu+pnu+hnu)/(4.*gsn)
                Apauli2(ppi,hpi,pnu,hnu)=Eppi+Epnu-factorp-factorn
              endif
  220 continue
c
c ************* Potential for optical model exciton model **************
c
c bonetti: subroutine for determination of effective absorption optical
c          potential
c
      if (preeqmode.eq.3) call bonetti
c
c ************* Blann's R-factor for multiple pre-equilibrium **********
c
c flagmulpre: flag for multiple pre-equilibrium calculation
c flag2comp : flag for two-component pre-equilibrium model
c k0        : index of incident particle
c Rblann    : Blann's factor
c
c The R-factor takes into account the neutron-proton distinction in
c multiple preequilibrium. The indices for the R-factor are
c Rblann(itype,type,nstep) where type is the second emitted particle,
c itype is the first emitted particle and nstep is the stage of primary
c pre-equilibrium emission. For this factor it is assumed that N=Z and
c that an effective interaction ratio of  p-n  to p-p or n-n of a
c factor 3 applies, see Blann and Vonach, Phys. Rev. C28, 1475 (1983).
c Eqs. (9) and (10) and the paragraph below. The numbers are taken
c from GNASH.
c
      if (.not.(flagmulpre.and..not.flag2comp)) return
c
c 1. Incident neutron
c
      if (k0.eq.1) then
        Rblann(2,2,1)=0.
        Rblann(2,1,1)=1.
        Rblann(2,2,2)=0.25
        Rblann(2,1,2)=0.75
        Rblann(2,2,3)=0.375
        Rblann(2,1,3)=0.625
        Rblann(2,2,4)=0.438
        Rblann(2,1,4)=0.562
        Rblann(2,2,5)=0.469
        Rblann(2,1,5)=0.531
        Rblann(1,1,1)=0.25
        Rblann(1,2,1)=0.75
        Rblann(1,1,2)=0.375
        Rblann(1,2,2)=0.625
        Rblann(1,1,3)=0.438
        Rblann(1,2,3)=0.562
        Rblann(1,1,4)=0.469
        Rblann(1,2,4)=0.531
        Rblann(1,1,5)=0.484
        Rblann(1,2,5)=0.516
      else
c
c 1. Incident proton
c
        Rblann(1,1,1)=0.
        Rblann(1,2,1)=1.
        Rblann(1,1,2)=0.25
        Rblann(1,2,2)=0.75
        Rblann(1,1,3)=0.375
        Rblann(1,2,3)=0.625
        Rblann(1,1,4)=0.438
        Rblann(1,2,4)=0.562
        Rblann(1,1,5)=0.469
        Rblann(1,2,5)=0.531
        Rblann(2,1,1)=0.75
        Rblann(2,2,1)=0.25
        Rblann(2,1,2)=0.625
        Rblann(2,2,2)=0.375
        Rblann(2,1,3)=0.562
        Rblann(2,2,3)=0.438
        Rblann(2,1,4)=0.531
        Rblann(2,2,4)=0.469
        Rblann(2,1,5)=0.516
        Rblann(2,2,5)=0.484
      endif
      do 410 i=1,2
        do 410 J=1,2
          do 410 k=6,numparx
            Rblann(i,J,k)=0.5
  410 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
