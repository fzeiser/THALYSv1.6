      subroutine incidentout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 25, 2012
c | Task  : Reaction output for incident channel
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,iang
c
c *********** Total cross sections for incident channel ****************
c
c k0       : index of incident particle
c xstotinc : total cross section (neutrons only) for incident channel
c xsreacinc: reaction cross section for incident channel
c xselasinc: total elastic cross section (neutrons only) for incident
c            channel
c
      write(*,'(/" Optical model results"/)')
      if (k0.eq.1) then
        write(*,'(" Total cross section   :",1p,e11.4," mb")') xstotinc
      endif
      write(*,'(" Reaction cross section:",1p,e11.4," mb")') xsreacinc
      if (k0.eq.1) then
        write(*,'(" Elastic cross section :",1p,e11.4," mb")') xselasinc
      endif
c
c For low energy neutrons we give the resonance parameters.
c
c Einc     : incident energy in MeV
c Atarget  : mass number of target nucleus
c Sstrength: s,p,d,etc-wave strength function
c Rprime   : potential scattering radius
c
      if (k0.eq.1.and.Einc.le.0.1) then
        write(*,'(/" S-wave and P-wave strength functions and ",
     +    "potential scattering radius"/)')
        write(*,'("      A      Value"/)')
        write(*,'(" S0:",i4,f8.4," .e-4")') Atarget,Sstrength(0)*1.e4
        write(*,'(" S1:",i4,f8.4," .e-4")') Atarget,Sstrength(1)*1.e4
        write(*,'(" R :",i4,f8.4," fm")') Atarget,Rprime
      endif
c
c *********** Transmission coefficients for incident channel ***********
c
c lmaxinc: maximal l-value for transmission coefficients for incident
c          channel
c parname: name of particle
c Tjlinc : transmission coefficients as a function of spin and l for
c          the incident channel
c Tlinc  : transmission coefficients as a function of l for the
c          incident channel, averaged over spin
c
      if (lmaxinc.eq.-1) return
      write(*,'(/" Transmission coefficients for incident ",a8,
     +  " at ",f7.3," MeV"/)') parname(k0),Einc
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
      if (k0.eq.1.or.k0.eq.2.or.k0.eq.4.or.k0.eq.5) then
        write(*,'("  L    T(L-1/2,L)   T(L+1/2,L)    Tav(L)"/)')
        do 10 l=0,lmaxinc
          write(*,'(1x,i3,1p,3e13.5)') l,Tjlinc(-1,l),Tjlinc(1,l),
     +      Tlinc(l)
   10   continue
      endif
c
c 2. Spin 1 particles: Deuterons
c
      if (k0.eq.3) then
        write(*,'("  L     T(L-1,L)      T(L,L)      T(L+1,L)",
     +    "     Tav(L)"/)')
        do 20 l=0,lmaxinc
          write(*,'(1x,i3,1p,4e13.5)') l,Tjlinc(-1,l),Tjlinc(0,l),
     +      Tjlinc(1,l),Tlinc(l)
   20   continue
      endif
c
c 3. Spin 0 particles: Alpha-particles
c
      if (k0.eq.6) then
        write(*,'("  L     T(L)"/)')
        do 30 l=0,lmaxinc
          write(*,'(1x,i3,1p,e13.5)') l,Tjlinc(0,l)
   30   continue
      endif
c
c 4. Photons
c
c gammax: number of l-values for gamma multipolarity
c
      if (k0.eq.0) then
        write(*,'("  L     T(L)"/)')
        do 40 l=1,gammax
          write(*,'(1x,i3,1p,e13.5)') l,Tjlinc(0,l)
   40   continue
      endif
c
c *********** Shape elastic scattering angular distribution ************
c
c flagang : flag for output of angular distributions
c nangle  : number of angles
c angle   : angle
c Ltarget : excited level of target
c directad: direct angular distribution
c
      if (flagang.and.k0.gt.0) then
        write(*,'(/" Shape elastic scattering angular distribution"/)')
        write(*,'(" Angle    Cross section"/)')
        do 110 iang=0,nangle
          write(*,'(1x,f5.1,1p,e16.5)') angle(iang),
     +      directad(k0,Ltarget,iang)
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
