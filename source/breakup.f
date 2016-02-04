      subroutine breakup
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 4, 2012
c | Task  : Contribution of breakup reactions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      integer      type,type2,nen
      real         Nab(numpar,numpar),r0,Deff,Ca,Cb,Ecent,arg,step,Sab,
     +             F,H,Feff,width,fac1,Emax,wplus,wmin,Epk,Tpren,factor,
     +             Ehalf,xsbreakup,Eout,wi,fac2,TE,PE
c
c ************************** Kalbach model *****************************
c
c Break-up model by Kalbach: FENDL-3 document, unpublished.
c
c Nab    : break-up normalization constants
c parskip: logical to skip outgoing particle
c k0     : index of incident particle
c
c We first determine the possible break-up channels and the other
c particle resulting from the break-up.
c
      Nab(3,1)=5.4
      Nab(3,2)=5.4
      Nab(4,1)=5.0
      Nab(4,2)=1.25
      Nab(4,3)=1.22
      Nab(5,1)=1.25
      Nab(5,2)=5.0
      Nab(5,3)=1.22
      Nab(6,1)=1.07
      Nab(6,2)=1.15
      Nab(6,3)=0.32
      Nab(6,4)=0.31
      Nab(6,5)=0.73
      do 10 type=1,6
        if (parskip(type)) goto 10
c
c (d,n) and (d,p)
c
        if (k0.eq.3) then
          if (type.ne.1.and.type.ne.2) goto 10
          if (type.eq.1) type2=2
          if (type.eq.2) type2=1
        endif
c
c (t,d) and (t,p)
c
        if (k0.eq.4) then
          if (type.ne.3) goto 10
          if (type.eq.1) type2=3
          if (type.eq.2) type2=1
          if (type.eq.3) type2=1
        endif
c
c (h,d) and (h,p)
c
        if (k0.eq.5) then
          if (type.gt.3) goto 10
          if (type.eq.1) type2=2
          if (type.eq.2) type2=3
          if (type.eq.3) type2=2
        endif
c
c (a,n), (a,p), (a,d), (a,t) and (a,h)
c
        if (k0.eq.6) then
          if (type.eq.6) goto 10
          if (type.eq.1) type2=5
          if (type.eq.2) type2=4
          if (type.eq.3) type2=3
          if (type.eq.4) type2=2
          if (type.eq.5) type2=1
        endif
c
c Calculation of terms independent of emission energy.
c
c r0       : effective radius parameter
c Einc     : incident energy in MeV
c Deff     : effective target-projectile separation
c Atarget  : mass number of target nucleus
c onethird : 1/3
c Ca,Cb    : effective Coulomb barrier
c parZ     : charge number of particle
c Ztarget  : charge number of target nucleus
c ZZ       : charge number of residual nucleus
c Ecent    : centroid energy for emission spectrum
c parA     : mass number of particle
c
c Calculation of terms independent of emission energy.
c
c Centroid energy
c
        r0=1.2+5./(1.+exp(Einc/30.))
        Deff=r0*(Atarget**onethird)+1.2
        Ca=1.44*parZ(k0)*Ztarget/Deff
        Cb=1.44*parZ(type)*ZZ(0,0,type)/Deff
        Ecent=parA(type)/real(parA(k0))*(Einc-Ca)+Cb
c
c Full width at half maximum
c
c arg      : help variable
c step     : step function
c Sab      : separation energy for projectile
c parmass  : mass of particle in a.m.u.
c amu      : atomic mass unit in MeV
c F        : full width at half maximum
c
        arg=parA(k0)-parA(type)-1.5
        if (arg.lt.0.) then
          step=0.
        else
          step=1.
        endif
        Sab=(parmass(type)+parmass(type2)-parmass(k0))*amu
        F=62.*(1.-1./exp(Einc/173.))*(1.-Atarget/(155.*Sab*Sab))-3.*step
c
c Effective full width at half maximum (asymmetric peaks)
c
c H         : help variable
c Feff      : effective full width at half maximum
c width     : width of break-up peak in emission spectrum
c fac1,fac2 : help variables
c sqrttwopi : sqrt(2.*pi)
c Emax      : maximal emission energy for particle channel
c eninccm   : center-of-mass incident energy in MeV
c Q         : Q-value
c wplus,wmin: half widths
c Epk       : peak energy
c
        H=0.5*F
        Emax=eninccm+Q(type)
        Feff=H+min(H,0.6*(Emax-Ecent))
        width=max(Feff/2.35,0.1*Einc)
        fac1=1./(width*sqrttwopi)
        wplus=max(0.,min(H,0.6*(Emax-Ecent)))/2.35
        wmin=max(0.,H-max(0.,0.6*(Ecent-Emax)))/2.35
        if (Emax.ge.Ecent-1.67*H.and.Emax.le.Ecent) then
          Epk=Emax
        else
          Epk=Ecent
        endif
c
c Breakup cross section
c
c Tpren    : barrier-penetrability factor
c Cbreak   : adjustable parameter for break-up reactions
c adjust   : subroutine for energy-dependent parameter adjustment
c factor   : multiplication factor
c Ehalf    : help variable
c twothird : 2/3
c xsbreakup: break-up cross section
c
        Ehalf=42.*(parA(k0)-parA(type))**twothird
        Tpren=1./(1.+exp((Ehalf-Einc)/14.))
        key='cbreak'
        call adjust(Einc,key,0,0,type,0,factor)
        xsbreakup=factor*Cbreak(type)*Nab(k0,type)*Deff*Deff*
     +    exp(Einc/170.)*Tpren
c
c Break-up term that depends on emission energy.
c
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c Eout,egrid: outgoing energy
c TE        : barrier penetrability
c PE        : spectrum function
c xspreeqbu : preequilibrium cross section per particle type and
c             outgoing energy for breakup
c
        do 110 nen=ebegin(type),eend(type)
          Eout=egrid(nen)
          if (Eout.le.Epk) then
            wi=max(wmin,0.05*Einc)
          else
            wi=max(wplus,0.05*Einc)
          endif
          fac2=1./(2.*(2.*wi)**2)
          if (type.gt.1) then
            TE=1./(1.+exp(3.*(Cb-Eout)/Cb))
          else
            TE=1.
          endif
          PE=fac1*exp(-(Eout-Epk)**2*fac2)*TE
          xspreeqbu(type,nen)=xsbreakup*PE
  110   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
