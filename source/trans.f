      subroutine trans(Zix,Nix,transm)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : November 28, 2013
c | Task  : Transmission coefficients per fission mode
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,i,icount
      real             xa,xb
      double precision transm,funcfismode,snew,sold
c
c **********************************************************************
c
      xa=0.0
      xb=excfis
      if (excfis.eq.0) then
        transm=funcfismode(Zix,Nix,0.)
      else
        icount=0
        i=0
        snew=10.0
        sold=1.0
   10   if (abs(snew/sold).le.0.99.or.abs(snew/sold).ge.1.01)
     +    then
          sold=snew
          i=i+1
          if (i.gt.10) goto 20
          call trapzd(Zix,Nix,xa,xb,snew,i)
          goto 10
        else
          icount=icount+1
          if (icount.eq.2) goto 20
          sold=snew
          i=i+1
          if (i.gt.10) goto 20
          call trapzd(Zix,Nix,xa,xb,snew,i)
          goto 10
        endif
 20     transm=snew
      endif
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
