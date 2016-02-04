      subroutine cascade(Zcomp,Ncomp,nex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 17, 2011
c | Task  : Gamma-ray cascade
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,nex,J,parity,i,k,Jres,Pres
      double precision xsJP,intens,xsgamma
c
c ******************* Gamma-ray cascade ********************************
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c nex        : excitation energy bin of compound nucleus
c jdis,J,Jres: spin of level
c parity,Pres: parity
c parlev     : parity of level
c xsJP,xspop : population cross section for spin and parity
c branchratio: gamma-ray branching ratio to level
c branchlevel: level to which branching takes place
c intens     : total gamma intensity
c xspopex    : population cross section summed over spin and parity
c xspartial  : emitted cross section flux per energy bin
c numZchan   : maximal number of outgoing proton units in individual
c              channel description
c numNchan   : maximal number of outgoing neutron units in individual
c              channel description
c mcontrib   : contribution to emission spectrum
c
      J=int(jdis(Zcomp,Ncomp,nex))
      parity=parlev(Zcomp,Ncomp,nex)
      xsJP=xspop(Zcomp,Ncomp,nex,J,parity)
      do 10 i=1,nbranch(Zcomp,Ncomp,nex)
        k=branchlevel(Zcomp,Ncomp,nex,i)
        Jres=int(jdis(Zcomp,Ncomp,k))
        Pres=parlev(Zcomp,Ncomp,k)
        intens=xsJP*branchratio(Zcomp,Ncomp,nex,i)
        xspop(Zcomp,Ncomp,k,Jres,Pres)=
     +    xspop(Zcomp,Ncomp,k,Jres,Pres)+intens
        xspopex(Zcomp,Ncomp,k)=xspopex(Zcomp,Ncomp,k)+intens
        xspopex(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)-intens
        xspartial(0,nex)=xspartial(0,nex)+intens
        if (Zcomp.le.numZchan.and.Ncomp.le.numNchan)
     +    mcontrib(0,nex,k)=intens
c
c ************ Storage of discrete gamma line intensities **************
c
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flagelectron: flag for application of electron conversion coefficient
c xsgamma     : help variable
c conv        : conversion coefficient
c xsgamdis    : discrete gamma-ray cross section
c xsgamdistot : total discrete gamma-ray cross section
c
        if (flaggamdis) then
          if (flagelectron) then
            xsgamma=intens/(1.+conv(Zcomp,Ncomp,nex,i))
          else
            xsgamma=intens
          endif
          xsgamdis(Zcomp,Ncomp,nex,k)=xsgamma
          xsgamdistot(Zcomp,Ncomp)=xsgamdistot(Zcomp,Ncomp)+xsgamma
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
