      subroutine msdtotal
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 10, 2004
c | Task  : Total multi-step direct cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,ns,iang
c
c ************** Angle-integrated multi-step cross sections ************
c
c msdall    : total multi-step direct cross section
c msdsum    : multi-step direct cross section summed over steps and
c             integrated over energy
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c msdtot    : multi-step direct cross section summed over steps
c maxmsd    : number of MSD steps
c msdstepint: n-step direct cross section integrated over energy
c msdstep   : continuum n-step direct cross section
c deltaE    : energy bin around outgoing energies
c flagddx   : flag for output of double-differential cross sections
c
      msdall=0.
      do 10 type=1,2
        msdsum(type)=0.
        do 20 nen=ebegin(type),eend(type)
          msdtot(type,nen)=0.
   20   continue
        do 30 ns=1,maxmsd
          msdstepint(type,ns)=0.
          do 40 nen=ebegin(type),eend(type)
            msdtot(type,nen)=msdtot(type,nen)+msdstep(type,ns,nen)
            msdstepint(type,ns)=msdstepint(type,ns)+
     +        msdstep(type,ns,nen)*deltaE(nen)
   40     continue
          msdsum(type)=msdsum(type)+msdstepint(type,ns)
   30   continue
        msdall=msdall+msdsum(type)
   10 continue
      if (.not.flagddx) return
c
c ******************* Multi-step angular distributions *****************
c
c nanglecont  : number of angles for continuum
c msdtotintad : multi-step direct angular distribution summed over steps
c               and integrated over energy
c msdtotad    : multi-step direct angular distribution summed over steps
c msdstepintad: n-step direct angular distribution integrated over
c               energy
c msdstepad   : continuum n-step direct angular distribution
c
c Total multi-step angular distributions
c
      do 110 type=1,2
        do 120 iang=0,nanglecont
          msdtotintad(type,iang)=0.
  120   continue
        do 130 nen=ebegin(type),eend(type)
          do 130 iang=0,nanglecont
            msdtotad(type,nen,iang)=0.
  130   continue
        do 140 ns=1,maxmsd
          do 150 iang=0,nanglecont
            msdstepintad(type,ns,iang)=0.
  150     continue
          do 160 iang=0,nanglecont
            do 170 nen=ebegin(type),eend(type)
              msdtotad(type,nen,iang)=msdtotad(type,nen,iang)+
     +          msdstepad(type,ns,nen,iang)
              msdstepintad(type,ns,iang)=msdstepintad(type,ns,iang)+
     +          msdstepad(type,ns,nen,iang)*deltaE(nen)
  170       continue
            msdtotintad(type,iang)=msdtotintad(type,iang)+
     +        msdstepintad(type,ns,iang)
  160     continue
  140   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
