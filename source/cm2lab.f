      subroutine cm2lab(type)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 12, 2004
c | Task  : Recoil for binary reaction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type
      real    coeff,vlabx,vlaby,vlab2,vlab
c
c EJECTILE TREATMENT
c
c The ejectile LAB angles and energies corresponding to the CM points
c are deduced
c
c Eejlab11   : LAB ejectile energy corresponding to (Eejcm1,angejcm1)
c cosejlab11 : LAB ejectile angle cosine corresponding to
c              (Eejcm1,angejcm1)
c Eejlab12   : LAB ejectile energy corresponding to (Eejcm1,angejcm2)
c cosejlab12 : LAB ejectile angle cosine corresponding to
c              (Eejcm1,angejcm2)
c Eejlab21   : LAB ejectile energy corresponding to (Eejcm2,angejcm1)
c cosejlab21 : LAB ejectile angle cosine corresponding to
c              (Eejcm2,angejcm1)
c Eejlab22   : LAB ejectile energy corresponding to (Eejcm2,angejcm2)
c cosejlab22 : LAB ejectile angle cosine corresponding to
c              (Eejcm2,angejcm2)
c
      if (type.eq.0.) then
        Eejlab11=Eejcm1
        Eejlab12=Eejcm1
        Eejlab21=Eejcm2
        Eejlab22=Eejcm2
        cosejlab11=cosejcm1
        cosejlab12=cosejcm2
        cosejlab21=cosejcm1
        cosejlab22=cosejcm2
        sinejlab11=sinejcm1
        sinejlab12=sinejcm2
        sinejlab21=sinejcm1
        sinejlab22=sinejcm2
      else
        coeff=0.5*ejectmass
        vlabx=vcm+vejcm1*cosejcm1
        vlaby=vejcm1*sinejcm1
        vlab2=vlabx**2+vlaby**2
        vlab=max(sqrt(vlab2),1.e-20)
        cosejlab11=vlabx/vlab
        sinejlab11=vlaby/vlab
        Eejlab11=coeff*vlab2
        vlabx=vcm+vejcm1*cosejcm2
        vlaby=vejcm1*sinejcm2
        vlab2=vlabx**2+vlaby**2
        vlab=max(sqrt(vlab2),1.e-20)
        cosejlab12=vlabx/vlab
        sinejlab12=vlaby/vlab
        Eejlab12=coeff*vlab2
        vlabx=vcm+vejcm2*cosejcm1
        vlaby=vejcm2*sinejcm1
        vlab2=vlabx**2+vlaby**2
        vlab=max(sqrt(vlab2),1.e-20)
        cosejlab21=vlabx/vlab
        sinejlab21=vlaby/vlab
        Eejlab21=coeff*vlab2
        vlabx=vcm+vejcm2*cosejcm2
        vlaby=vejcm2*sinejcm2
        vlab2=vlabx**2+vlaby**2
        vlab=max(sqrt(vlab2),1.e-20)
        cosejlab22=vlabx/vlab
        sinejlab22=vlaby/vlab
        Eejlab22=coeff*vlab2
      endif
c
c RECOIL TREATMENT
c
c The Recoil LAB angles and energies corresponding to the CM points
c are deduced
c
c vlabx      : LAB recoil velocity projection on x-axis // to CM
c              velocity
c vlaby      : LAB recoil velocity projection on y-axis
c vlab       : LAB recoil velocity
c cosreclab11: LAB recoil angle cosine corresponding to
c              (Eejcm1,angejcm1)
c Ereclab11  : LAB recoil energy corresponding to (Eejcm1,angejcm1)
c cosreclab12: LAB recoil angle cosine corresponding to
c              (Eejcm1,angejcm2)
c Ereclab12  : LAB recoil energy corresponding to (Eejcm1,angejcm2)
c cosreclab21: LAB recoil angle cosine corresponding to
c              (Eejcm2,angejcm1)
c Ereclab21  : LAB recoil energy corresponding to (Eejcm2,angejcm1)
c cosreclab22: LAB recoil angle cosine corresponding to
c              (Eejcm2,angejcm2)
c Ereclab22  : LAB recoil energy corresponding to (Eejcm2,angejcm2)
c
      coeff=0.5*recoilmass
      vlabx=vcm-vreccm1*cosejcm1
      vlaby=-vreccm1*sinejcm1
      vlab2=vlabx**2+vlaby**2
      vlab=max(sqrt(vlab2),1.e-20)
      cosreclab11=vlabx/vlab
      sinreclab11=-sinejcm1
      Ereclab11=coeff*vlab2
      vlabx=vcm-vreccm1*cosejcm2
      vlaby=-vreccm1*sinejcm2
      vlab2=vlabx**2+vlaby**2
      vlab=max(sqrt(vlab2),1.e-20)
      cosreclab12=vlabx/vlab
      sinreclab12=-sinejcm2
      Ereclab12=coeff*vlab2
      vlabx=vcm-vreccm2*cosejcm1
      vlaby=-vreccm2*sinejcm1
      vlab2=vlabx**2+vlaby**2
      vlab=max(sqrt(vlab2),1.e-20)
      cosreclab21=vlabx/vlab
      sinreclab21=-sinejcm1
      Ereclab21=coeff*vlab2
      vlabx=vcm-vreccm2*cosejcm2
      vlaby=-vreccm2*sinejcm2
      vlab2=vlabx**2+vlaby**2
      vlab=max(sqrt(vlab2),1.e-20)
      cosreclab22=vlabx/vlab
      sinreclab22=-sinejcm2
      Ereclab22=coeff*vlab2
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
