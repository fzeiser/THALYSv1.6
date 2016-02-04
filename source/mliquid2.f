      double precision function mliquid2(Z,A)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 3, 2007
c | Task  : Goriely liquid drop mass
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Z,A,N
      double precision avol,as,asym,ass,ac,Mn,Mh,rA,rZ,rN,factor,Ev,
     +                 Esur,Esym,Ecoul,Eldm
c
c ****************** Spherical Goriely parameters **********************
c
c mliquid2: function for liquid drop mass
c rZ,Z    : charge number of residual nucleus
c rA,A    : mass number of residual nucleus
c avol,as.: Goriely parameters
c Mn,Mh   : mass excess of neutron and hydrogen
c rN,N    : neutron number of residual nucleus
c factor  : help variable
c Ev      : volume energy
c Esur    : surface energy
c twothird: 2/3
c Esym    : symmetry energy
c onethird: 1/3
c Ecoul   : Coulomb energy
c Eldm    : liquid drop energy
c amu     : atomic mass unit in MeV
c
c Goriely model: S. Goriely, ND2001, Tsukuba, Japan, J. Nuc. Sci Techn.
c August 2002, suppl 2. p.536 (2002)
c
      data avol,as,asym,ass,ac /-15.6428,17.5418,27.9418,-25.3440,0.70/
      data Mn,Mh/ 8.07132281,7.2889694/
      N=A-Z
      rA=real(A)
      rZ=real(Z)
      rN=real(N)
      factor=(rN-rZ)/rA
      Ev=avol*rA
      Esur=as*(rA**twothird)
      Esym=(asym+ass*(rA**(-onethird)))*rA*factor**2
      Ecoul=ac*rZ**2/(rA**onethird)
      Eldm=Z*Mh+N*Mn+Ev+Esur+Esym+Ecoul-1.4333e-5*rZ**2.39
      mliquid2=A+Eldm/amu
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
