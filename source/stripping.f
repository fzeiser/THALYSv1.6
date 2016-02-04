      subroutine stripping
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Vivian Demetriou
c | Date  : April 4, 2012
c | Task  : Contribution of stripping and pickup reactions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      logical      surfwell
      integer      type,ndelta,ndeltapi,ndeltanu,ppi,hpi,pnu,hnu,A,Z,N,
     +             i,j,nen,k,l
      real         proj2sp1,ejecmass,ejec2sp1,term1,Kap,Va,term2,term3,
     +             term4,base,term5,factor,termps,V1well,surface,XNT,
     +             gsn,gsp,Ewell,Eout,factor1,Eres,P,preeqpair,omegaNT,
     +             omegaph,phdens2
c
c ************************** Kalbach model *****************************
c
c The stripping and pickup model is described in
c C. Kalbach, "Preequilibrium reactions with complex channels",
c Phys. Rev. C71, 034606 (2005).
c
c projmass : mass of projectile
c parmass  : mass of particle in a.m.u.
c k0       : index of incident particle
c proj2sp1 : 2*spin +1 of projectile
c parspin  : spin of particle
c parskip  : logical to skip outgoing particle
c
c Factors for pickup and stripping processes. For reactions involving
c only neutrons and protons this is thus not considered.
c
      projmass=parmass(k0)
      proj2sp1=2.*parspin(k0)+1.
      do 10 type=1,6
        if (parskip(type)) goto 10
        if (k0.le.2.and.type.le.2) goto 10
c
c Calculation of terms independent of emission energy.
c
c ejecmass: mass of ejectile
c ejec2sp1: 2*spin +1 of ejectile
c term1,..: help variables
c Kap     : alpha enhancement factor
c eninccm : center-of-mass incident energy in MeV
c Einc    : incident energy in MeV
c ndelta  : number of transferred particles
c parA    : mass number of particle
c ndeltapi: number of transferred protons
c parZ    : charge number of particle
c ndeltanu: number of transferred neutrons
c parN    : neutron number of particle
c Va      : potential drop
c
c Kap=5 for outgoing helions was introduced to better fit (n,h) data,
c i.e. it does not come from the original Kalback model.
c
        ejecmass=parmass(type)
        ejec2sp1=2.*parspin(type)+1.
        term1=ejec2sp1*ejecmass/(proj2sp1*projmass)
        Kap=1.
        if ((k0.eq.1.or.k0.eq.2)) then
          if (type.eq.6) Kap=12.
          if (type.eq.3.and.Einc.gt.80.) Kap=80./Einc
          if (type.eq.5) Kap=5.
        endif
        if (k0.eq.6.and.(type.eq.1.or.type.eq.2))
     +    Kap=12.-11.*max(eninccm-20.,0.)/eninccm
        ndelta=abs(parA(k0)-parA(type))
        ndeltapi=parZ(k0)-parZ(type)
        ndeltanu=parN(k0)-parN(type)
        Va=12.5*projmass
        term2=(Kap/projmass)*(projmass/(Einc+Va))**(2*ndelta)
c
c Initial configuration for pickup, stripping or t-h charge exchange
c
c ppi: proton particle number
c hpi: proton hole number
c pnu: neutron particle number
c hnu: neutron hole number
c
        ppi=max(ndeltapi,0)
        hpi=max(-ndeltapi,0)
        pnu=max(ndeltanu,0)
        hnu=max(-ndeltanu,0)
c
c Further terms
c
c AA,A   : mass number of residual nucleus
c ZZ,Z   : charge number of residual nucleus
c NN,N   : neutron number of residual nucleus
c base   : help variable
c Ztarget: charge number of target nucleus
c Atarget: mass number of target nucleus
c termps : term for pickup and stripping
c adjust : subroutine for energy-dependent parameter adjustment
c Cstrip : adjustable parameter for stripping/pick-up reactions
c
        A=AA(0,0,type)
        Z=ZZ(0,0,type)
        N=NN(0,0,type)
        if (k0.eq.1) then
          term3=(5500./real(A))**ndelta
        else
          term3=(3800./real(A))**ndelta
        endif
        if (parA(k0).lt.parA(type)) term4=1./(80.*eninccm)
        if (parA(k0).gt.parA(type)) term4=1./(580.*sqrt(eninccm))
        if (parA(k0).eq.parA(type)) term4=1./(1160.*sqrt(eninccm))
        base=real(2.*Ztarget)/real(Atarget)
        term5=base**(2*(parZ(k0)+2)*hpi+2*pnu)
        key='cstrip'
        call adjust(Einc,key,0,0,type,0,factor)
        termps=factor*Cstrip(type)*term1*term2*term3*term4*term5
c
c XNT function
c
c V1well,Ewell: depth of potential well
c XNT         : probability of exciting particle-hole pair
c gsn,gsp,gsa : single-particle level density parameter
c Ntarget     : neutron number of target nucleus
c Kph         : constant for single-particle level density parameter
c               (g=A/Kph)
c
        V1well=17.
        if (k0.eq.1) V1well=surface(1,Einc)
        if (k0.eq.2) V1well=surface(2,Einc)
        if (k0.eq.5.or.k0.eq.6) V1well=25.
        XNT=sqrt(min(Einc,100.)/projmass)*7./(V1well*Atarget*Atarget)*
     +    (pnu**2+ppi**2+hnu**2+1.5*hpi**2)
        gsn=N/Kph
        gsp=Z/Kph
        if (ndeltapi.eq.0) then
          Ewell=V1well*base
        else
          Ewell=V1well
        endif
c
c Calculation of stripping or pickup spectra.
c
c ebegin   : first energy point of energy grid
c eend     : last energy point of energy grid
c factor1  : help variable
c Eres     : total energy of residual system
c Etotal   : total energy of compound system (target + projectile)
c S        : separation energy per particle
c preeqpair: pre-equilibrium pairing energy
c pairmodel: model for preequilibrium pairing energy
c
        do 110 nen=ebegin(type),eend(type)
          Eout=egrid(nen)
          factor1=xsreac(type,nen)*Eout
          P=preeqpair(parZ(type),parN(type),ndelta,Etotal,pairmodel)
          Eres=Etotal-S(0,0,type)-Eout-P
c
c Check if outgoing energy exceeds maximal possible energy
c
          if (Eres.lt.0.) goto 110
c
c Stripping/pick-up terms that depend on emission energy.
c For reactions in which only one hole is left, e.g. (p,d),
c the finite well function leads to a discontinuity at the well
c depth. Therefore, for this case the well depth is set equal to
c the Fermi energy.
c
c omegaNT  : state density function for pickup and stripping
c surfwell : flag for surface effects in finite well
c omegaph  : particle-hole state density
c phdens2  : function for two-component particle-hole state density
c xspreeqps: preequilibrium cross section per particle type and
c            outgoing energy for pickup and stripping
c
          omegaNT=0.
          surfwell=.false.
          do 120 i=0,3
            do 120 j=0,3-i
              omegaph=phdens2(parZ(type),parN(type),ppi+i,hpi+i,pnu+j,
     +          hnu+j,gsp,gsn,Eres,Ewell,surfwell)
              omegaNT=omegaNT+XNT**(i+j)*omegaph
  120     continue
          do 130 i=0,ppi
            do 130 j=0,hpi
              do 130 k=0,pnu
                do 130 l=0,hnu
                  if (i+j+k+l.ne.0) omegaNT=omegaNT+
     +              phdens2(parZ(type),parN(type),ppi-i,hpi-j,pnu-k,
     +              hnu-l,gsp,gsn,Eres,Ewell,surfwell)
  130     continue
          xspreeqps(type,nen)=termps*omegaNT*factor1
  110   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
