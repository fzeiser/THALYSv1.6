      subroutine onestepA(type)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 1, 2008
c | Task  : Unnormalized one-step direct cross sections for outgoing
c |         energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Zix,Nix,h,p,nen,J,nen2,na,nb,nc,iang
      real    gs,Eout,ignatyuk,rJ,omegaJ(0:numJmsd),omega,Ea,Eb,Ec,
     +        total,xsa,xsb,xsc,xsi,xs
c
c ************* Calculate continuum one-step cross sections ************
c
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c h          : hole number
c p          : particle number
c g,gs       : single-particle level density parameter
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c Eout       : outgoing energy
c egrid      : outgoing energy grid
c eninccm    : center-of-mass incident energy in MeV
c Exmsd      : excitation energy for MSD energy grid
c Q          : Q-value for target nucleus
c flaggshell : flag for energy dependence of single particle level
c              density parameter g
c ignatyuk   : function for energy dependent level density parameter a
c alev       : level density parameter
c maxJmsd    : maximal spin for MSD calculation
c rJ         : help variable
c omegaJ     : help variable
c omega      : particle-hole state density
c locate     : subroutine to find value in ordered table
c Emsd       : MSD energy grid
c msdbins2   : number of energy points for MSD calculation
c na,nb,nc   : help variables
c Ea,Eb,Ec   : help variables
c total        help variable
c xsa,xsb,xsc: help variables
c xs,xsi     : help variables
c xsdwin     : DWBA cross section as a function of incident energy,
c              outgoing energy and angular momentum
c pol2       : subroutine for polynomial interpolation of second order
c msdstep1   : continuum one-step direct cross section (unnormalized)
c flagddx    : flag for output of double-differential cross sections
c nanglecont : number of angles for continuum
c xsdw       : DWBA angular distribution as a function of incident
c              energy, outgoing energy, angular momentum and angle
c msdstepad1 : continuum one-step direct angular distribution
c              (unnormalized)
c
      Zix=Zindex(0,0,type)
      Nix=Nindex(0,0,type)
      h=1
      p=1
      gs=g(Zix,Nix)
      do 10 nen=ebegin(type),eend(type)
        Eout=egrid(nen)
        if (Eout.gt.eninccm) goto 10
        Exmsd=eninccm-Eout+Q(type)
        if (flaggshell) gs=g(Zix,Nix)*ignatyuk(Zix,Nix,Exmsd,0)/
     +    alev(Zix,Nix)
        do 20 J=0,maxJmsd
          rJ=real(J)
          omegaJ(J)=omega(Zix,Nix,p,h,gs,Exmsd,rJ)
   20   continue
        call locate(Emsd,0,msdbins2,Eout,nen2)
        if (nen2.gt.1) then
          na=nen2-1
          nb=nen2
          nc=nen2+1
        else
          na=nen2
          nb=nen2+1
          nc=nen2+2
        endif
        Ea=Emsd(na)
        Eb=Emsd(nb)
        Ec=Emsd(nc)
        total=0.
        do 30 J=0,maxJmsd
          xsa=log(max(xsdwin(0,na,J,0),1.e-30))
          xsb=log(max(xsdwin(0,nb,J,0),1.e-30))
          xsc=log(max(xsdwin(0,nc,J,0),1.e-30))
          call pol2(Ea,Eb,Ec,xsa,xsb,xsc,Eout,xsi)
          xs=exp(xsi)
          if (xs.lt.1.e-30) xs=0.
          total=total+omegaJ(J)*xs
   30   continue
        msdstep1(type,nen)=total
        if (flagddx) then
          do 40 iang=0,nanglecont
            total=0.
            do 50 J=0,maxJmsd
              xsa=log(max(xsdw(0,na,J,iang,0),1.e-30))
              xsb=log(max(xsdw(0,nb,J,iang,0),1.e-30))
              xsc=log(max(xsdw(0,nc,J,iang,0),1.e-30))
              call pol2(Ea,Eb,Ec,xsa,xsb,xsc,Eout,xsi)
              xs=exp(xsi)
              if (xs.lt.1.e-30) xs=0.
              total=total+omegaJ(J)*xs
   50       continue
            msdstepad1(type,nen,iang)=total
   40     continue
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
