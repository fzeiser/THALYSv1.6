      subroutine opticalalpha
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : March 21, 2012
c | Task  : Other optical potential for alphas
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c Alternative options for alpha OMP
c
c S. Goriely: inclusion of the alpha OMP of Mc Fadden & Satchler
c for alphaomp=2, folding model for alphaomp=3,4,5
c
c Overwrite some of the previous values.
c
c v,rv,...: optical model parameters
c
      v=185.0
      rv=1.40
      av=0.52
      w=25.0
      rw=rv
      aw=av
      vd=0.
      wd=0.
      vso=0.
      wso=0.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
