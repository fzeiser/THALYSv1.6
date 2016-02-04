      subroutine constants
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 27, 2012
c | Task  : Constants and initialization
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer i
c
c ************************* Derived constants **************************
c
c twopi      : 2.*pi
c pi         : pi
c pi2        : pi**2
c sqrttwopi  : sqrt(2.*pi)
c fourpi     : 4.*pi
c deg2rad    : conversion factor for degrees to radians
c rad2deg    : conversion factor for radians to degrees
c onethird   : 1/3
c twothird   : 2/3
c twopihbar  : 2*pi/hbar
c hbar       : Planck's constant / 2.pi in MeV.s
c hbarc      : hbar.c in MeV.fm
c clight     : speed of light in vacuum in m/s
c pi2h2c2    : 1/(pi*pi*clight*clight*hbar**2) in mb**-1.MeV**-2
c pi2h3c2    : 1/(pi*pi*clight*clight*hbar**3) in mb**-1.MeV**-3.s**-1
c amupi2h3c2 : amu/(pi*pi*clight*clight*hbar**3) in mb**-1.MeV**-2.s**-1
c amu        : atomic mass unit in MeV
c amu4pi2h2c2: amu/(4*pi*pi*clight*clight*hbar**2) in mb**-1.MeV**-1
c
      twopi=2.*pi
      pi2=pi*pi
      sqrttwopi=sqrt(2.*pi)
      fourpi=4.*pi
      deg2rad=pi/180.
      rad2deg=180./pi
      onethird=1./3.
      twothird=2.*onethird
      twopihbar=twopi/hbar
      hbarc=hbar*clight*1.e15
      pi2h2c2=0.1/(pi2*hbarc*hbarc)
      pi2h3c2=pi2h2c2/hbar
      amupi2h3c2=real(amu)*pi2h3c2
      amu4pi2h2c2=real(amu)*0.25*pi2h2c2
c
c ************************** Set signs *********************************
c
c sgn   : +1 for even argument, -1 for odd argument
c numl  : maximum l-value (set in talys.cmb)
c pardis: parity distribution
c
c sgn is used in level density and compound nucleus calculations
c
      sgn(0)=1.
      do 10 i=2,2*numl,2
        sgn(i)=1.
        sgn(i-1)=-1.
   10 continue
c
c ATTENTION: As long as we use an equidistant parity distribution, we
c            set it equal to 0.5, instead of calling a subroutine.
c            As soon as an analytical non-equidistant parity
c            distribution is used in TALYS, pardis should be changed
c            into a function pardis(J,Ex) throughout TALYS.
c
      pardis=0.5
c
c ************************ Set counter for isotope *********************
c
c iso      : counter for isotope
c numiso   : maximum number of isotopes per element
c natstring: string extension for file names
c
      iso=1
      do 110 i=1,numiso
        natstring(i)='    '
  110 continue
c
c *********************** Set default mass for fission *****************
c
c fislim: mass above which nuclide fissions
c
      fislim=215
c
c ***************** Set maximum acceptable energy for TALYS ************
c
c Emaxtalys: maximum acceptable energy for TALYS
c
      Emaxtalys=1000.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
