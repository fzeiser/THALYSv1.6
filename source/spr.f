      subroutine spr
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : S, P and R' resonance parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real    Efac,Rpot,r2k2
c
c **************** S, P and R' resonance parameters ********************
c
c Rprime     : potential scattering radius
c xselasinc  : total elastic cross section (neutrons only) for incident
c              channel
c fourpi     : 4.*pi
c Efac       : factor with incident energy
c Einc       : incident energy in MeV
c twopi      : 2.*pi
c Rpot       : standard value for R
c Atarget    : mass number of target nucleus
c onethird   : 1/3
c r2k2       : help variable
c wavenum    : wave number
c Sstrength  : s,p,d,etc-wave strength function
c Tlinc      : transmission coefficients as a function of l for the
c              incident channel, averaged over spin
c flagendf   : flag for information for ENDF-6 file
c flagendfdet: flag for detailed ENDF-6 information per channel
c nin        : counter for incident energy
c numinclow  : number of incident energies below Elow
c Ztarget    : charge number of target nucleus
c
      Rprime=10.*sqrt(0.001*max(xselasinc,0.)/fourpi)
      Efac=1./(sqrt(1.e6*Einc)*twopi)
      Rpot=1.35*Atarget**onethird
      r2k2=Rpot*Rpot*wavenum*wavenum
      Sstrength(0)=Tlinc(0)*Efac
      Sstrength(1)=Tlinc(1)*Efac*(1.+r2k2)/r2k2
      Sstrength(2)=Tlinc(2)*Efac*(9.+3.*r2k2+r2k2*r2k2)/(r2k2*r2k2)
      if (flagendf.and.flagendfdet.and.
     +  (Einc.le.0.1.or.nin.eq.numinclow+1)) then
        open (unit=1,status='unknown',file='spr.opt')
        write(1,'(2i4,3f8.4)') Atarget,Ztarget,Sstrength(0)*1.e4,
     +    Sstrength(1)*1.e4,Rprime
        close (unit=1)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
