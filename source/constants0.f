      block data constants0
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 23, 2012
c | Task  : Constants and basic properties of particles
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer indx
c
c ****************** General properties of particles *******************
c
c Indices: fission =-1
c          photon  = 0
c          neutron = 1
c          proton  = 2
c          deuteron= 3
c          triton  = 4
c          helium-3= 5
c          alpha   = 6
c
c parname: name of particle
c parsym : symbol of particle
c parZ   : charge number of particle
c parN   : neutron number of particle
c parA   : mass number of particle
c parmass: mass of particle in a.m.u.
c excmass: mass excess of particle in a.m.u.
c parspin: spin of particle
c
      data (parname(indx),indx=-1,6) /'fission ','gamma   ','neutron ',
     +  'proton  ','deuteron','triton  ','helium-3','alpha   '/
      data (parsym(indx),indx=-1,6)  /'f','g','n','p','d','t','h','a'/
      data (parZ(indx),indx=0,6)    /0,0,1,1,1,2,2/
      data (parN(indx),indx=0,6)    /0,1,0,1,2,1,2/
      data (parA(indx),indx=0,6)    /0,1,1,2,3,3,4/
      data (parmass(indx),indx=0,6) /0.,1.008664904,1.007276470,
     +  2.013553214,3.016049268,3.016029310,4.002603250/
      data (excmass(indx),indx=0,6) /0.,8.664923e-3,7.825032e-3,
     +  1.4101778e-2,1.6049268e-2,1.6029310e-2,2.603250e-3/
      data (parspin(indx),indx=0,6) /0.,0.5,0.5,1.,0.5,0.5,0./
c
c ************************ Nuclear symbols *****************************
c
c nuc    : symbol of nucleus
c numelem: number of elements
c magic  : magic numbers
c
      data (nuc(indx),indx=1,numelem) /
     +  'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     +  'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     +  'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     +  'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     +  'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     +  'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     +  'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +  'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     +  'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     +  'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     +  'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/
      data (magic(indx),indx=1,8) /2,8,20,28,50,82,126,184/
c
c **************** Set character symbol for parities *******************
c
c The '+' parity will have the value 1 and the '-' parity the value -1
c
c cparity: parity (character)
c
      data (cparity(indx),indx=-1,1,2) /'-','+'/
c
c *********************** Fundamental constants ************************
c
c pi      : pi
c amu     : atomic mass unit in MeV
c e2      : square of elementary charge in MeV.fm
c hbar    : Planck's constant / 2.pi in MeV.s
c clight  : speed of light in vacuum in m/s
c kT      : energy kT expressed in MeV corresponding to a
c           temperature T9=1
c emass   : electron mass in MeV/c^2
c avogadro: Avogadro's number
c qelem   : elementary charge in C
c
      data pi /3.14159265358979323/
      data amu,e2 /931.49386,1.439965161/
      data hbar,clight /6.5821220e-22,2.99792458e8/
      data kT,emass /0.086173,0.510999/
      data avogadro,qelem /6.0221367e23,1.6021773e-19/
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
