      function tripathi(zproj,aproj,iz,ia,e)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning (adapted from R.K. Tripathi)
c | Date  : December 18, 2013
c | Task  : Semi-empirical reaction cross section of Tripathi et al.
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer zproj,aproj,iz,ia,xzt,xat
      real    tripathi,e,pi,onethird,fourthird,rp,radius,rt,vp,vt,
     +        dens,const,const1,t1,gcm,bcm,plab,ecmp,ecmt,rela,ecm,
     +        bigr,bigb,xm,x1,sl,ce,term1,delta,beta,twxsec,xabs,expo,
     +        expo1
c
c ******************* Reaction cross section calculation ***************
c
c zproj         : charge number of projectile
c aproj         : mass number of projectile
c iz            : charge number of target nucleus
c ia            : mass number of target nucleus
c e             : incident energy in MeV per nucleon
c pi            : pi
c onethird      : 1/3
c fourthird     : 4/3
c rp            : projectile radius
c rt            : target radius
c all other var.: ask Tripathi
c
      data pi /3.14159265358979323/
      onethird=1./3.
      fourthird=4./3.
      tripathi=0.
      rp=radius(real(aproj))
      rt=radius(real(ia))
      vp=fourthird*pi*rp**3
      vt=fourthird*pi*rt**3
      dens=0.5*((aproj/vp)+(ia/vt))
      const=1.75*dens/8.824728e-02
      if (iz.lt.zproj) then
        xzt=zproj
        xat=aproj
        zproj=iz
        aproj=ia
        zproj=xzt
        aproj=xat
      endif
      if (aproj.eq.1) const=2.05
      if (zproj.eq.2.and.aproj.eq.4)
     +  const1=2.77-ia*8.0e-03+(ia*ia)*1.8e-05
      if (zproj.eq.3) const=const*onethird
      t1=40.
      if (zproj.eq.0) then
        if (ia.ge.11.and.ia.lt.40) t1=30.
        if (iz.eq.14) t1=35.
        if (iz.eq.26) t1=30.
      endif
      gcm=(aproj*(1.+e/938.)+ia)/
     +  sqrt(aproj**2+ia**2+2.*aproj*(e+938.)*ia/938.)
c
c A. Koning: Safety
c
      if (gcm.le.1.) return
      bcm=sqrt(1.-1./gcm**2)
      plab=aproj*sqrt(2.*938.*e+e*e)
      ecmp=gcm*(e+938.)*aproj-bcm*gcm*plab-aproj*938.
      ecmt=gcm*938.*ia-ia*938.
      rela=ecmp+ecmt
      ecm=rela
      bigr=rp+rt+1.2*(aproj**onethird+ia**onethird)/(ecm**onethird)
      bigb=1.44*zproj*iz/bigr
      if (zproj.eq.1.and.ia.gt.56) bigb=0.90*bigb
      if (aproj.gt.56.and.iz.eq.1) bigb=0.90*bigb
      if (aproj.eq.1.and.ia.eq.12) bigb=3.5*bigb
      if (aproj.eq.1) then
        if (ia.le.16.and.ia.ge.13) bigb=(ia/7.)*bigb
        if (iz.eq.12) bigb=1.8*bigb
        if (iz.eq.14) bigb=1.4*bigb
        if (iz.eq.20) bigb=1.3*bigb
      endif
      if (aproj.eq.1.and.ia.lt.4) bigb=21.0*bigb
      if (aproj.lt.4.and.ia.eq.1) bigb=21.0*bigb
      if (aproj.eq.1.and.ia.eq.4) bigb=27.0*bigb
      if (aproj.eq.4.and.ia.eq.1) bigb=27.0*bigb
      if (zproj.eq.0.or.iz.eq.0) bigb=0.0
      xm=1.
      if (zproj.eq.0) then
        if (ia.lt.200.) then
          x1=2.83-3.1e-02*ia+1.7e-04*ia*ia
          if (x1.le.1) x1=1.
          sl=1.0
          if (ia.eq.12) sl=1.6
          if (ia.lt.12) sl= 0.6
          xm=(1-x1*exp(-e/(sl*x1)))
        else
          xm=(1.-0.3*exp(-(e-1.)/15.))*(1.-exp(-(e-0.9)))
        endif
      endif
      if (zproj.eq.2.and.aproj.eq.4) const=const1-
     +  0.8/(1.+exp((250.-e)/75.))
      expo=min((e-20)/10.,80.)
      if (zproj.eq.1.and.aproj.eq.1) then
        if (ia.gt.45) t1=40.+ia/3.
        if (ia.lt.4) t1=55.
        const= 2.05 - 0.05/(1.+exp((250.-e)/75.))
        if (ia.lt.4) const=1.7
        if (iz.eq.12) then
          t1=40.
          const = 2.05 -3.0/(1.+exp(expo))
        endif
        if (iz.eq.14) then
          t1=40.
          const=2.05-1.75/(1.+exp(expo))
        endif
        if (iz.eq.18) then
          t1=40.
          const=2.05-2.0/(1.+exp(expo))
        endif
        if (iz.eq.20) then
          t1=40.
          expo1=min((e-40)/10.,80.)
          const=2.05-1.0/(1.+exp(expo1))
        endif
      endif
      if (zproj.eq.0.and.aproj.eq.1) then
        const=2.*(0.134457/dens)
        if (ia.gt.140.and.ia.lt.200) const=const-1.5*(ia-2.*iz)/ia
        if (ia.lt.60) const=const-1.5*(ia-2.*iz)/ia
        if (ia.le.40) const=const+0.25/(1.+exp(-(170.-e)/100.))
        if (iz.gt.82) const=const-real(iz)/(ia-iz)
        if (iz.ge.82) then
          expo1=min((e-20)/20.,80.)
          const=const-2.0/(1.+exp(expo))
        endif
        if (iz.le.20.and.iz.ge.10)
     +    const=const-1.0/(1.+exp(expo))
      endif
      ce=const*(1.-exp(-e/t1))-0.292*exp(-e/792)*cos(0.229*e**0.453)
      term1=(ia*aproj)**onethird/(ia**onethird+aproj**onethird)
      delta=1.615*term1-0.873*ce
      delta=delta+0.140*term1/ecm**onethird
      delta=delta+0.794*(ia-2.*iz)*zproj/(ia*aproj)
      delta=-delta
      beta=1.
      twxsec=10.*pi*1.26*1.26*beta*(0.873*aproj**onethird+
     +  0.873*ia**onethird-delta)**2
      xabs=twxsec*(1.-bigb/ecm)*xm
      if (xabs.lt.0) xabs=0.
      tripathi=xabs
      return
      end
