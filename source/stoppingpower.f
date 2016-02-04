      subroutine stoppingpower(E,dEdx)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 23, 2012
c | Task  : Calculate stopping power
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real E,dEdx,pmass,enum,denom,beta,Imean,gam,eta,Wmax,term1,term2
c
c ********************* Velocity of projectile *************************
c
c E         : beam energy in MeV
c dEdx      : stopping power in MeV/cm
c pmass     : rest mass of projectile in MeV (m_0*c^2)
c parmass   : mass of particle in a.m.u.
c k0        : index of incident particle
c amu       : atomic mass unit in MeV
c enum,denom: help variables
c beta      : v/c
c
      pmass=parmass(k0)*amu
      enum=E*(E+2.*pmass)
      denom=(E+pmass)**2
      beta=sqrt(enum/denom)
c
c ********************* Bethe-Bloch formula ****************************
c
c The formula for the stopping power can be found on page 24 of
c W.R. Leo, "Techniques for nuclear and particle physics experiments",
c Springer Verlag, Berlin, 1994
c
c Imean      : mean excitation potential in MeV
c Ztarget    : charge number of target nucleus
c gam        : gamma
c eta        : beta*gam
c Wmax       : maximum energy transfer in knock-on collision in MeV
c emass      : electron mass in MeV/c^2
c rhotarget  : target density
c Atarget    : mass number of target nucleus
c parZ       : charge number of particle
c term1,term2: help variables
c
      Imean=Ztarget*(9.76+58.8*Ztarget**(-1.19))*1.e-6
      gam=1./sqrt(1.-beta**2)
      eta=beta*gam
      Wmax=2.*emass*(eta**2)
      term1=0.1535*rhotarget*Ztarget/real(Atarget)*
     +  (parZ(k0)**2)/(beta**2)
      term2=log(2.*emass*(eta**2)*Wmax/(Imean**2))-2.*beta**2
      dEdx=term1*term2
      return
      end
Copyright (C) 2010  A.J. Koning
