      subroutine trapzd(Zix,Nix,abeg,bend,snew,nk)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Integration
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c  (C) Copr. 1986-92 Numerical Recipes Software %#..5#2P.
c
c ****************** Declarations and common blocks ********************
c
      integer          nk,Zix,Nix,it,jk
      real             abeg,bend,x,del
      double precision sum,tnm,snew,funcfismode
c
c **********************************************************************
c
      snew=0.0
      if (nk.eq.1) then
        snew=0.5*(bend-abeg)*(funcfismode(Zix,Nix,abeg)+
     +    funcfismode(Zix,Nix,bend))
      else
        it=2**(nk-2)
        tnm=it
        del=(bend-abeg)/tnm
        x=abeg+0.5*del
        sum=0.
        do 10 jk=1,it
          sum=sum+funcfismode(Zix,Nix,x)
          x=x+del
 10    continue
        snew=0.5*(snew+(bend-abeg)*sum/tnm)
      endif
      return
      end
