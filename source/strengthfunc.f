      subroutine strengthfunc(Tinc,l,j)
c
c +---------------------------------------------------------------------
c | Author: Gilles Noguere
c | Date  : December 18, 2011
c | Task  : (l,j) neutron strength function for URR calculations
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,j,i
      real    Tinc,ac,P(0:numl),Shift(0:numl),ri,den,num,vl
c
c ******************** Initialization parameters ***********************
c
c strengthlj: (l,j) neutron strength function
c Tinc      : neutron transmission coefficient
c ac        : channel radius in ENDF-6 convention
c tarmass   : mass of target nucleus
c onethird  : 1/3
c wavenum   : wave number
c P(0)      : neutron penetrability factor for s-wave
c Shift(0)  : hard sphere shift factor for s-wave
c
      ac=1.23*tarmass**onethird+0.8
      P(0)=wavenum*ac
      Shift(0)=0.
c
c *********************  Penetrability factor **************************
c
c num,den : help variables
c P(l)    : l-dependent penetrability factor
c Shift(l): l-dependent hard sphere shift factor
c
      do 10 i=1,l
        ri=real(i)
        num=P(i-1)*P(0)**2
        den=P(i-1)**2+(ri-Shift(i-1))**2
        P(i)=num/den
        num=(ri-Shift(i-1))*P(0)**2
        den=P(i-1)**2+(ri-Shift(i-1))**2
        Shift(i)=num/den-ri
   10 continue
c
c *********************** Strength function ****************************
c
c vl   : penetrability factor
c nulj : degree of freedom
c Einc : incident energy in MeV
c twopi: 2.*pi
c
      vl=P(l)/(P(0)*10.)
      strengthlj(l,j)=Tinc/(twopi*vl*sqrt(Einc))
      return
      end
