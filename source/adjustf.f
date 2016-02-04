      subroutine adjustF(E,nrange,en1,en2,Dr,sr,factor)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : February 1, 2012
c | Task  : Local parameter adjustment
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nrange,nr
      real    E,en1(numrange),en2(numrange),Dr(numrange),sr(numrange),
     +        factor,elow,eup,emid,D,sigma,expo,offset
c
c ************************* Woods-Saxon function ***********************
c
c factor     : Woods-Saxon multiplication factor
c E          : incident energy
c nrange     : number of energy ranges for local adjustment
c en1        : start energy of local adjustment
c en2        : end energy of local adjustment
c Dr         : depth of local adjustment
c sr         : variance of local adjustment
c offset     : offset to guarantee continuity
c elow,eup,..: help variables
c
      factor=1.
      do 10 nr=1,nrange
        elow=en1(nr)
        eup=en2(nr)
        if (E.gt.elow.and.E.lt.eup) then
          emid=0.5*(elow+eup)
          D=0.01*Dr(nr)
          sigma=sr(nr)
          if (sigma.eq.0.) sigma=(eup-emid)/2.
          expo=(eup-emid)**2/(2.*sigma**2)
          offset=0.
          if (expo.le.80.) offset=-D*exp(-expo)
          expo=(E-emid)**2/(2.*sigma**2)
          factor=1.+offset
          if (expo.le.80.) factor=1.+D*exp(-expo)+offset
          return
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
