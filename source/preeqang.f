      subroutine preeqang
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 10, 2004
c | Task  : Pre-equilibrium angular distribution
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,iang,NL,i
      real    Eout,ang,kalbach,xs
c
c ************** Kalbach angular distribution for exciton model ********
c
c parskip   : logical to skip outgoing particle
c preeqmode : designator for pre-equilibrium model
c xspreeqtot: preequilibrium cross section per particle type
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c xspreeq   : preequilibrium cross section per particle type
c             and outgoing energy
c egrid,Eout: outgoing energy grid
c nanglecont: number of angles for continuum
c ang       : angle
c anglecont : angle in degrees for continuum
c deg2rad   : conversion factor for degrees to radians
c xspreeqad : preequilibrium angular distribution per particle type
c kalbach   : Kalbach function
c Einc      : incident energy in MeV
c
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (preeqmode.eq.4.and.(type.eq.1.or.type.eq.2)) goto 10
        if (xspreeqtot(type).eq.0.) goto 10
        do 20 nen=ebegin(type),eend(type)
          if (xspreeq(type,nen).eq.0.) goto 20
          Eout=egrid(nen)
          do 30 iang=0,nanglecont
            ang=anglecont(iang)*deg2rad
            xspreeqad(type,nen,iang)=xspreeq(type,nen)*
     +        kalbach(type,Einc,Eout,ang)
   30     continue
   20   continue
   10 continue
c
c ************ Pre-equilibrium cross sections for direct states ********
c
c Correction of pre-equilibrium cross sections for direct discrete
c cross sections. If the cross sections for discrete states have NOT
c been calculated  by a direct reaction model, we collapse the
c continuum pre-equilibrium cross sections in the high energy region
c on the associated discrete states.
c
c Nlast,NL   : last discrete level
c parZ       : charge number of particle
c parN       : neutron number of particle
c xs         : help variable
c xspreeqdisc: preequilibrium cross section for discrete state
c eoutdis    : outgoing energy of discrete state reaction
c nangle     : number of angles
c angle      : angle in degrees
c deg2rad    : conversion factor for degrees to radians
c directad   : direct angular distribution

      do 110 type=0,6
        if (parskip(type)) goto 110
        NL=Nlast(parZ(type),parN(type),0)
        do 120 i=0,NL
          xs=xspreeqdisc(type,i)
          if (xs.eq.0.) goto 120
          Eout=eoutdis(type,NL)
          do 130 iang=0,nangle
            ang=angle(iang)*deg2rad
            directad(type,i,iang)=directad(type,i,iang)+
     +        xs*kalbach(type,Einc,Eout,ang)
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
