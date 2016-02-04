      subroutine soukhovitskii(k,Z,A,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : March 21, 2012
c | Task  : Global optical model parameters for actinides by
c |         Soukhovitskii et al.
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer k,Z,A
      real    eopt,asym,eferm,f,Cviso,V0r,Var,vrdisp,v1r,v2r,lambdaR,
     +        viso,Ccoul,phicoul,rr,Crr,widr,w1loc,w2loc,wddisp,Cwiso,
     +        Wad,d1loc,d2loc,d3loc,vso1loc,vso2loc,wso1loc,wso2loc
c
c *** Parameters of Soukhovitskii et al, J. Phys. G30, p. 905 (2004) ***
c
c k         : designator for particle
c Z         : charge number of residual nucleus
c A         : mass number of residual nucleus
c eopt      : incident energy
c asym      : asymmetry parameter
c eferm     : Fermi energy
c f         : eopt-efer
c Cviso,....: optical potential parameters
c w1loc.....: help variables
c Fv1,......: adjustable factors for OMP
c v1,v2,v3  : components for V
c w1,w2     : components for W
c d1,d2,d3  : components for Wd
c mw,md     : powers for W and Wd
c vso1,vso2 : components for Vso
c wso1,wso2 : components for Wso

      asym=(A-2.*Z)/real(A)
      if (k.eq.1) then
        eferm=-11.2814+0.02646*A
      else
        eferm=-8.4075+0.01378*A
      endif
      f=eopt-eferm
      Cviso=10.5
      V0r=-41.45
      Var=-0.06667
      vrdisp=92.44
      v1r=0.03
      v2r=2.05e-4
      lambdaR=3.9075e-3
      viso=1.+((-1)**k)*Cviso*asym/(V0r+Var*(A-232.)+vrdisp)
      v=(V0r+Var*(A-232.)+v1r*f+v2r*(f**2)+vrdisp*exp(-lambdaR*f))*viso
      if (k.eq.2) then
        Ccoul=0.9
        phicoul=(lambdaR*vrdisp*exp(-lambdaR*f)-v1r-2.*v2r*f)*viso
        v=v+Ccoul*Z/(A**onethird)*phicoul
      endif
      v=Fv1*v
      rr=1.245
      Crr=0.05
      widr=100.
      rv=Frv*rr*(1.-Crr*f**2/(f**2+widr**2))
      av=Fav*(0.660+2.53e-4*eopt)
      w1loc=Fw1*14.74
      w2loc=Fw2*81.63
      w=w1loc*f**2/(f**2+w2loc**2)
      rw=Frw*1.2476
      aw=Faw*0.594
      vd=0.
      rvd=Frvd*1.2080
      avd=Favd*0.614
      wddisp=17.38
      Cwiso=24.
      Wad=0.03833
      d1loc=Fd1*(wddisp+Wad*(A-232.)+((-1)**k)*Cwiso*asym)
      d2loc=Fd2*0.01759
      d3loc=Fd3*11.79
      wd=d1loc*f**2*exp(-d2loc*f)/(f**2+d3loc**2)
      rwd=Frwd*1.2080
      awd=Fawd*0.614
      vso1loc=Fvso1*5.86
      vso2loc=Fvso2*0.0050
      vso=vso1loc*exp(-vso2loc*f)
      rvso=Frvso*1.1213
      avso=Favso*0.59
      wso1loc=-3.1*Fwso1
      wso2loc=Fwso2*160.
      wso=wso1loc*f**2/(f**2+wso2loc**2)
      rwso=Frwso*1.1213
      awso=Fawso*0.59
      if (k.eq.1) then
        rc=0.
      else
        rc=Frc*1.2643
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
