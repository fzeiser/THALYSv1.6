      subroutine basicinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 25, 2009
c | Task  : Initialization of arrays for basic cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nen,type,l,ispin
c
c *************************** Initialization ***************************
c
c numen     : maximum number of outgoing energies
c lmax      : maximal l-value for transmission coefficients
c gammax    : number of l-values for gamma multipolarity
c xstot     : total cross section (neutrons only) for
c xsreac    : reaction cross section
c xsopt     : optical model reaction cross section
c xselas    : total elastic cross section (neutrons only)
c numl      : maximum l-value (set in talys.cmb)
c Tl        : transmission coefficients as a function of particle type,
c             energy and l-value (averaged over spin)
c Tjl       : transmission coefficients as a function of particle type,
c             energy, spin and l-value
c
      do 10 nen=0,numen
        do 10 type=0,6
          if (type.eq.0) then
            lmax(type,nen)=0
          else
            lmax(type,nen)=gammax
          endif
          xstot(type,nen)=0.
          xsreac(type,nen)=0.
          xsopt(type,nen)=0.
          xselas(type,nen)=0.
   10 continue
      do 20 l=0,numl
        do 20 nen=0,numen
          do 20 type=0,6
            Tl(type,nen,l)=0.
   20 continue
      do 30 l=0,numl
        do 30 ispin=-1,1
          do 30 nen=0,numen
            do 30 type=0,6
              Tjl(type,nen,ispin,l)=0.
   30 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
