      function thill(ehw,bhw,whw)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : August 25, 2004
c | Task  : Hill-Wheeler penetrability
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real thill,ehw,bhw,whw,expo
c
c ************************* Hill-Wheeler formula ***********************
c
c thill: Hill-Wheeler penetrability
c ehw  : E(compound nucleus) - E(transition state)
c bhw  : fission barrier height
c whw  : fission barrier width
c expo : exponent
c twopi: 2.*pi
c
      thill=0.
      expo=twopi*(ehw-bhw)/whw
      if (expo.gt.-80.) thill=1./(1.+exp(-expo))
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
