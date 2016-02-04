      subroutine multistepB
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : Multi-step direct cross sections on outgoing energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,nen2,na,nb,nc,ns,iang
      real    Eout,Ea,Eb,Ec,xsa,xsb,xsc,xs
c
c ********** Interpolate continuum multi-step cross sections ***********
c
c parskip    : logical to skip outgoing particle
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c Eout       : outgoing energy
c egrid      : outgoing energy grid
c eninccm    : center-of-mass incident energy in MeV
c locate     : subroutine to find value in ordered table
c Emsd       : MSD energy grid
c msdbins2   : number of energy points for MSD calculation
c na,nb,nc   : help variables
c Ea,Eb,Ec   : help variables
c maxmsd     : number of MSD steps
c xsa,xsb,xsc: help variables
c msdstep0   : n-step cross section for MSD
c pol2       : subroutine for polynomial interpolation of second order
c xs         : help variable
c msdstep    : continuum n-step direct cross section
c flagddx    : flag for output of double-differential cross sections
c nanglecont : number of angles for continuum
c msdstepad0 : n-step angular distribution for MSD
c msdstepad  : continuum n-step direct angular distribution
c
      do 10 type=1,2
        if (parskip(type)) goto 10
        do 20 nen=ebegin(type),eend(type)
          Eout=egrid(nen)
          if (Eout.gt.eninccm) goto 20
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
          do 30 ns=2,maxmsd
            xsa=max(msdstep0(type,ns,na),1.e-30)
            xsb=max(msdstep0(type,ns,nb),1.e-30)
            xsc=max(msdstep0(type,ns,nc),1.e-30)
            call pol2(Ea,Eb,Ec,xsa,xsb,xsc,Eout,xs)
            if (xs.lt.1.e-30) xs=0.
            msdstep(type,ns,nen)=xs
            if (flagddx) then
              do 40 iang=0,nanglecont
                xsa=max(msdstepad0(type,ns,na,iang),1.e-30)
                xsb=max(msdstepad0(type,ns,nb,iang),1.e-30)
                xsc=max(msdstepad0(type,ns,nc,iang),1.e-30)
                call pol2(Ea,Eb,Ec,xsa,xsb,xsc,Eout,xs)
                if (xs.lt.1.e-30) xs=0.
                msdstepad(type,ns,nen,iang)=xs
   40         continue
            endif
   30     continue
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
