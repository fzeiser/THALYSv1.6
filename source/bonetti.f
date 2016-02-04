      subroutine bonetti
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : October 31, 2013
c | Task  : Determination of effective absorption optical potential
c +---------------------------------------------------------------------
c
c R. Bonetti, S. Galbiati, and L. Milazzo-Colli, Lett. Nuov. Cim. 18,
c 557 (1977)
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real    dr,emax,e,sum1,sum2,radv,radd,rr,expo,fwsvol,fwssurf,
     +        term1,term2,wjlm,Rjlm
      integer Zix,Nix,nrbins,Z,nenend,k,nen,i
c
c ********** Calculation of volume absorptive optical potential ********
c
c Zix,Zindex : charge number index for residual nucleus
c Nix,Nindex : neutron number index for residual nucleus
c k0         : index of incident particle
c dr         : integration bin width
c emax       : maximum excitation energy
c enincmax   : maximum incident energy
c S          : separation energy
c nenend,expo: help variables
c jlmexist   : flag for existence of tabulated radial matter density
c nrbins     : number of continuum excitation energy bins
c e          : energy
c ZZ,Z       : charge number of residual nucleus
c parZ       : charge number of particle
c prodZ      : product of charges of projectile and target nucleus
c mom        : subroutine for microscopic optical model (Eric Bauge)
c optical    : subroutine for determination of optical potential
c radv       : imaginary volume nuclear radius
c radd       : imaginary surface nuclear radius
c Atarget    : mass number of target nucleus
c onethird   : 1/3
c rhojlmn    : density for neutrons
c rhojlmp    : density for protons
c normjlm    : JLM potential normalization factors
c potjlm     : JLM potential depth values
c w,rw,aw    : imaginary volume potential, radius, diffuseness
c wd,rwd,awd : imaginary surface potential, radius, diffuseness
c rr         : running variable in integration over the radius
c fwsvol     : imaginary volume optical potential form factor
c fwssurf    : derivative of the imaginary surface optical potential
c              form factor
c term1,term2: help variables
c wjlm       : help variable
c Rjlm       : JLM normalization factor
c sum1,sum2  : help variables
c wvol       : absorption part of the optical potential averaged over
c              the volume
c
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      emax=enincmax+S(0,0,k0)
      nenend=10*min(numen,int(emax+1.))
      Rjlm=1.
      do 10 k=1,2
c
c A. Phenomenological potential (always done, also for the case of
c    normalizing the JLM OMP)
c
        nrbins=50
        dr=20./nrbins
        do 20 nen=-80,nenend
          if (jlmexist(Zix,Nix,k).and.nen.gt.1) goto 20
          e=0.1*real(nen)
          call optical(Zix,Nix,k,e)
          radv=rw*Atarget**onethird
          radd=rwd*Atarget**onethird
          sum1=0.
          sum2=0.
          do 30 i=1,nrbins
            rr=(i-0.5)*dr
            expo=(rr-radv)/aw
            if (expo.le.80.) then
              fwsvol=1./(1.+exp(expo))
            else
              fwsvol=0.
            endif
            expo=(rr-radd)/awd
            if (expo.le.80.) then
              fwssurf=-exp(expo)/(awd*(1.+exp(expo)**2))
            else
              fwssurf=0.
            endif
            term2=fwsvol*(rr**2)*dr
            term1=term2*(w*fwsvol-4.*awd*wd*fwssurf)
            sum1=sum1+term1
            sum2=sum2+term2
 30       continue
          if (sum2.ne.0.) wvol(k,nen)=sum1/sum2
 20     continue
c
c B. JLM potential
c
        if (jlmexist(Zix,Nix,k)) then
          nrbins=122
          dr=12.2/nrbins
          Z=ZZ(Zix,Nix,k)
          prodZ=Z*parZ(k)
          do 40 nen=nenend,-80,-1
c
c For negative energies, the phenomenological OMP is adopted, normalized
c at the first positive energy point.
c
            if (nen.lt.1) then
              wvol(k,nen)=Rjlm*wvol(k,nen)
              goto 40
            endif
            e=0.1*real(nen)
            call mom(Zix,Nix,dble(prodZ),dble(e))
            sum1=0.
            sum2=0.
            do 50 i=1,nrbins
              rr=(i-0.5)*dr
              if (k.eq.1) then
                term2=rhojlmn(Zix,Nix,i,1)*(rr**2)*dr
              else
                term2=rhojlmp(Zix,Nix,i,1)*(rr**2)*dr
              endif
              term1=-term2*normjlm(Zix,Nix,2)*potjlm(Zix,Nix,i,2)
              sum1=sum1+term1
              sum2=sum2+term2
 50         continue
            if (nen.eq.1) then
              Rjlm=1.
              if (sum2.ne.0.) then
                wjlm=sum1/sum2
                Rjlm=wjlm/wvol(k,nen)
              endif
            endif
            if (sum2.ne.0.) wvol(k,nen)=sum1/sum2
 40       continue
        endif
        do 60 nen=-200,-81
          wvol(k,nen)=wvol(k,-80)
 60     continue
 10   continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
