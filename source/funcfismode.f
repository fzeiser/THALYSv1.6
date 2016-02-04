      function funcfismode(Zix,Nix,Epar)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 13, 2005
c | Task  : Transmission coefficient per fission mode
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix
      real             Epar,temps(9),ald,tmp,bft,hwt,ignatyuk
      double precision density,funcfismode,expo
      data temps /0.0,0.3,0.6,0.9,1.2,1.6,2.0,2.5,3.0/
c
c **********************************************************************
c
      ald=ignatyuk(Zix,Nix,Epar,0)
      tmp=sqrt(Epar/ald)
      call splint(temps,bf,bfsplin,numtemp,tmp,bft)
      call splint(temps,hw,hwsplin,numtemp,tmp,hwt)
      expo=2.*pi*(bft+Epar-excfis)/hwt
      if (expo.le.80.) then
        funcfismode=1./(1.+exp(expo))
      else
        funcfismode=exp(-80.)/(exp(-80.)+exp(expo-80.))
      endif
      if (Epar.gt.0.)
     +  funcfismode=funcfismode*density(Zix,Nix,Epar,0.,1,0,1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
