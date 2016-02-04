      subroutine multistepA
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : Multi-step direct cross sections for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,ns,nenout,nenint,iang,iangint,iangout
      real    Cmulti,total,factor1,factor2,dang,tot1,fac1,fac2,angint
c
c ************** Angle-integrated multi-step cross sections ************
c
c Cmulti     : constant for second and higher steps
c specmass   : specific mass for target nucleus and all particles
c parZ       : charge number of particle
c k0         : index for incident particle
c parN       : neutron number of particle
c amu4pi2h2c2: amu/(4*pi*pi*clight*clight*hbar**2) in mb**-1.MeV**-1
c parskip    : logical to skip outgoing particle
c maxmsd     : number of MSD steps
c msdbins2   : number of energy points for MSD calculation
c msdstep0   : n-step cross section for MSD
c factor1-2  : help variables
c Emsd       : MSD energy grid
c xscont     : continuum one-step direct cross section for MSD
c parmass    : mass of particle in a.m.u.
c dEmsd      : energy bin for MSD
c
      Cmulti=real(specmass(parZ(k0),parN(k0),k0)*amu4pi2h2c2)
      do 10 type=1,2
        if (parskip(type)) goto 10
        do 20 ns=2,maxmsd
          do 30 nenout=1,msdbins2
            msdstep0(type,ns,nenout)=0.
            total=0.
            do 40 nenint=1,nenout-1
              factor1=msdstep0(type,ns-1,nenout)*Emsd(nenint)*
     +          xscont(type,type,nenint,nenout)
              factor2=msdstep0(k0,ns-1,nenout)*Emsd(nenint)*
     +          xscont(k0,type,nenint,nenout)
              total=total+0.5*(factor1+factor2)
   40       continue
            msdstep0(type,ns,nenout)=total*real(parmass(type))*Cmulti*
     +        dEmsd
   30     continue
   20   continue
   10 continue
c
c ******************* Multi-step angular distributions *****************
c
c flagddx    : flag for output of double-differential cross sections
c dang       : angle step
c pi         : pi
c nanglecont : number of angles for continuum
c msdstepad0 : n-step angular distribution for MSD
c iangout    : outgoing angle index
c iangint    : intermediate angle index
c fac1-2,tot1: help variables
c nangleint  : number of possibilities to link intermediate angle to
c              final angle
c angint     : intermediate angle
c deg2rad    : conversion factor for degrees to radians
c xscontad   : continuum one-step direct angular distribution for MSD
c anglecont  : angle in degrees for continuum
c
      if (.not.flagddx) return
      dang=pi/nanglecont
      do 110 type=1,2
        if (parskip(type)) goto 110
        do 120 ns=2,maxmsd
          do 130 nenout=1,msdbins2
            do 140 iang=0,nanglecont
              msdstepad0(type,ns,nenout,iang)=0.
  140       continue
            do 150 iangout=0,nanglecont
              total=0.
              do 160 nenint=1,nenout-1
                do 170 iangint=0,nanglecont
                  tot1=0.
                  fac1=msdstepad0(type,ns-1,nenint,iangint)*
     +              Emsd(nenint)
                  fac2=msdstepad0(k0,ns-1,nenint,iangint)*
     +              Emsd(nenint)
                  do 180 iang=0,nanglecont
                    factor1=fac1*nangleint(iangout,iangint,iang)*
     +                xscontad(type,type,nenint,nenout,iang)
                    factor2=fac2*nangleint(iangout,iangint,iang)*
     +                xscontad(k0,type,nenint,nenout,iang)
                    tot1=tot1+0.5*(factor1+factor2)
  180             continue
                  angint=anglecont(iangint)*deg2rad
                  total=total+tot1*dang*dang*sin(angint)
  170           continue
                msdstepad0(type,ns,nenout,iangout)=total*
     +            real(parmass(type))*Cmulti*dEmsd
  160         continue
  150       continue
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
