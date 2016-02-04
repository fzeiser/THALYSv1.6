      subroutine compoundinit
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 9, 2004
c | Task  : Initialization of compound model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,k
c
c ************ Help variable for Hauser-Feshbach subroutine ************
c
c spin2  : 2 * spin of particle (usually)
c parspin: spin of particle
c
c For alpha-particles and photons, spin2 is set to 1, to prevent
c division by zero in Hauser-Feshbach routines. This is fine, since in
c these cases the enumerator in the expression in which spin2 appears
c is always zero.
c
      spin2(0)=1
      do 10 type=1,5
        spin2(type)=int(2.*parspin(type))
   10 continue
      spin2(6)=1
c
c *********** Initialization for width fluctuation calculation *********
c
c enincmin       : minimum incident energy
c ewfc           : off-set incident energy for width fluctuation
c wpower         : power used for rho*(t**wpower)
c wmode          : designator for width fluctuation model
c nmold          : number of points for Gauss-Laguerre integration
c gaulag         : subroutine for calculation of Gauss-Laguerre arrays
c xmold,wmold    : variables for Gauss-Laguerre integration
c ngoep,ngoes,...: number of points for Gauss-Legendre integration
c xgoep,wgoep,...: variables for Gauss-Legendre integration
c gauleg         : subroutine for calculation of Gauss-Legendre arrays
c
c Generate weight and nodes for Gauss-Legendre/Laguerre integration
c
      if (enincmin.le.ewfc) then
        wpower=1
c
c 1. Moldauer
c
        if (wmode.eq.1) then
          nmold=32
          call gaulag(nmold,xmold,wmold)
        endif
c
c 2. HRTW
c
        if (wmode.eq.2) wpower=2
c
c 3. GOE
c
        if (wmode.eq.3) then
          wpower=5
          ngoep=50
          ngoes=50
          ngoet=50
          call gauleg(ngoep,xgoep,wgoep)
          call gauleg(ngoes,xgoes,wgoes)
          call gauleg(ngoet,xgoet,wgoet)
        endif
      endif
c
c ***** Initialization for compound nucleus angular distributions ******
c
c flagang: flag for output of angular distributions
c numfact: number of terms for factorial logarithm
c logfact: factorial logarithm
c
      if (flagang) then
        logfact(1) = 0.
        logfact(2) = 0.
        do 110 k=3,numfact
          logfact(k)=logfact(k-1)+log(float(k-1))
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
