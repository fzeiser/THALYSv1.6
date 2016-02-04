      subroutine sumrules
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 1, 2004
c | Task  : Giant resonance sum rules
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,Zix,Nix,i
      real    Sgmr,Sgqr,Sgor,Sleor,Sheor,betasq,deform1,edis1,jdis1
c
c ******************* Strength of giant resonances *********************
c
c Sgmr    : sum rule for giant monopole resonance
c Atarget : mass number of target nucleus
c Egrcoll : energy of giant resonance
c Ggrcoll : width of giant resonance
c Sgqr    : sum rule for giant quadrupole resonance
c onethird: 1/3
c twothird: 2/3
c Sgor    : sum rule for giant octupole resonance
c Sleor   : sum rule for low energy octupole resonance
c Sheor   : sum rule for high energy octupole resonance
c
c Using sum rules for giant resonances, the collective strength in the
c continuum is determined.
c
      Sgmr=23.*Atarget**(-5./3.)
      Egrcoll(0,1)=18.7-0.025*Atarget
      Ggrcoll(0,1)=3.
      Sgqr=575.*Atarget**(-5./3.)
      Egrcoll(2,1)=65.*Atarget**(-onethird)
      Ggrcoll(2,1)=85.*Atarget**(-twothird)
      Sgor=1208.*Atarget**(-5./3.)
      Sleor=0.3*Sgor
      Egrcoll(3,1)=31.*Atarget**(-onethird)
      Ggrcoll(3,1)=5.
      Sheor=0.7*Sgor
      Egrcoll(3,2)=115.*Atarget**(-onethird)
      Ggrcoll(3,2)=9.3-Atarget/48.
c
c Subtract the deformation parameters of the collective low lying
c states from the sum.
c
c Zcomp  : charge number index for compound nucleus
c Ncomp  : neutron number index for compound nucleus
c Zindex : charge number index for residual nucleus
c Nindex : neutron number index for residual nucleus
c k0     : index of incident particle
c numlev2: maximum number of levels
c deform : deformation parameter
c deftype: deformation length (D) or parameter (B)
c edis   : energy of level
c jdis   : spin of level
c
      Zcomp=0
      Ncomp=0
      Zix=Zindex(Zcomp,Ncomp,k0)
      Nix=Nindex(Zcomp,Ncomp,k0)
      do 10 i=1,numlev2
        deform1=deform(Zix,Nix,i)
        if (deform1.ne.0.) then
          if (deftype(Zix,Nix).eq.'D')
     +      deform1=deform1/(1.24*Atarget**onethird)
          betasq=deform1*deform1
          edis1=edis(Zix,Nix,i)
          jdis1=int(jdis(Zix,Nix,i))
          if (jdis1.eq.0) Sgmr=Sgmr-betasq*edis1
          if (jdis1.eq.2) Sgqr=Sgqr-betasq*edis1
          if (jdis1.eq.3) Sleor=Sleor-betasq*edis1
        endif
   10 continue
c
c Determine final GR deformation parameters.
c
c betagr: deformation parameter for giant resonance
c
      if (Sgmr.gt.0.) betagr(0,1)=sqrt(Sgmr/Egrcoll(0,1))
      if (Sgqr.gt.0.) betagr(2,1)=sqrt(Sgqr/Egrcoll(2,1))
      if (Sleor.gt.0.) betagr(3,1)=sqrt(Sleor/Egrcoll(3,1))
      betagr(3,2)=sqrt(Sheor/Egrcoll(3,2))
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
