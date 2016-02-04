      function evap(rn)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : evap=0 determines the number of evaporated neutrons for a
c |         given nucleus with mass and charges numbers amm and zee
c |         with excitation energy ess.
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This subroutine is based on the subroutine evap originally developed
c by U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real rn,evap,b,dumm,dum,an,bn,an1,bn1
c
c **********************************************************************
c
      call bdef(amm,zee,0.,dum,dumm,b)
      an=amm-rn
      call bdef(an,zee,0.,dum,dumm,bn)
      an1=an-1.
      call bdef(an1,zee,0.,dum,dumm,bn1)
      evap=ess-0.5*(bn+bn1)+b-0.621*rn*sqrt(rn+1.)
      return
      end
