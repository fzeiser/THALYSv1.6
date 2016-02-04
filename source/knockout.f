      subroutine knockout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Vivian Demetriou
c | Date  : April 4, 2012
c | Task  : Contribution of knockout and complex inelastic reactions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      logical      flagknock,flaginel
      integer      type,ndelta,type2,nen
      real         proj2sp1,ejecmass,ejec2sp1,termki,AKO,Ccl,phi,denom,
     +             Pn(6),gscomp(6),emax,dE,total,Eout,sigav,denomki(6),
     +             term1,factor,Ckn,termk0,termin0,factor1,Eres,P,
     +             preeqpair,U
c
c ************************** Kalbach model *****************************
c
c The knockout model is described in
c C. Kalbach, "Preequilibrium reactions with complex channels",
c Phys. Rev. C71, 034606 (2005).
c
c projmass: mass of projectile
c parmass : mass of particle in a.m.u.
c proj2sp1: 2*spin +1 of projectile
c k0      : index of incident particle
c parskip : logical to skip outgoing particle
c
      projmass=parmass(k0)
      proj2sp1=2.*parspin(k0)+1.
      do 10 type=1,6
        if (parskip(type)) goto 10
        if (k0.le.2.and.type.le.2) goto 10
c
c Factors for knockout and inelastic processes.
c
c ejecmass : mass of ejectile
c ejec2sp1 : 2*spin +1 of ejectile
c ndelta   : number of transferred particles
c parA     : mass number of particle
c flagknock: flag for knockout
c flaginel : flag for inelastic scattering
c termki   : term for knockout and inelastic
c AKO      : Pauli correction factor
c Ccl      : adjustment factor
c phi      : time fraction for alpha cluster
c Ztarget  : charge number of target nucleus
c Atarget  : mass number of target nucleus
c Ntarget  : neutron number of target nucleus
c denom,Pn : help variables
c gscomp   : single-particle level density parameter
c
c Calculation of terms independent of emission energy.
c
        ejecmass=parmass(type)
        ejec2sp1=2.*parspin(type)+1.
        ndelta=abs(parA(k0)-parA(type))
        flagknock=((k0.eq.1.or.k0.eq.2).and.type.eq.6)
        flaginel=(k0.eq.type)
        if (.not.(flagknock.or.flaginel)) goto 10
        termki=0.
        AKO=0.
        Ccl=1./14.
        phi=0.08
        if (Ntarget.gt.116.and.Ntarget.lt.126)
     +    phi=0.02+0.06*(126-Ntarget)/10.
        if (Ntarget.ge.126.and.Ntarget.lt.129)
     +    phi=0.02+0.06*(Ntarget-126)/3.
        denom=Atarget-2.*phi*Ztarget+0.5*phi*Ztarget
        Pn(1)=(Ntarget-phi*Ztarget)/denom
        Pn(2)=(Ztarget-phi*Ztarget)/denom
        Pn(6)=0.5*phi*Ztarget/denom
        gscomp(1)=Atarget/13.
        gscomp(2)=Atarget/13.
        gscomp(3)=Atarget/52.
        gscomp(4)=Atarget/156.
        gscomp(5)=Atarget/156.
        gscomp(6)=Atarget/208.
        if (flagknock)
     +    AKO=1./(2.*(gscomp(k0)**2))+1./(2.*(gscomp(6)**2))
c
c Denominator of cluster emission formula
c
c denomki   : denominator for knockout formula
c emax      : maximal emission energy for particle channel
c eninccm   : center-of-mass incident energy in MeV
c Q         : Q-value for target nucleus
c dE,total  : help variables
c coulbar   : Coulomb barrier
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c egrid,Eout: energies of basic energy grid in MeV
c sigav     : average cross section for emission channel
c xsreac    : reaction cross section
c deltaE    : energy bin around outgoing energies
c
        do 20 type2=1,6
          denomki(type2)=0.
          if (parskip(type2)) goto 20
          emax=eninccm+Q(type2)
          dE=emax-coulbar(type2)
          if (dE.gt.2.) then
            total=0.
            do 30 nen=ebegin(type2),eend(type2)
              Eout=egrid(nen)
              if (Eout.lt.coulbar(type2)) goto 30
              total=total+xsreac(type2,nen)*deltaE(nen)
   30       continue
            sigav=total/dE
          else
            sigav=xsreac(type2,eend(type2))
          endif
          denomki(type2)=(2.*parspin(type2)+1.)*sigav*
     +      (emax+2.*coulbar(type2))*max((emax-coulbar(type2))**2,1.)
   20   continue
c
c Knockout and inelastic terms
c
c adjust        : subroutine for energy-dependent parameter adjustment
c termk0,termin0: help variable
c Cknock        : adjustable parameter for knockout reactions
c xsreacinc     : reaction cross section for incident channel
c
        key='cknock'
        call adjust(Einc,key,0,0,type,0,factor)
        Ckn=factor*Cknock(type)
        if (flagknock) then
          termk0=projmass*denomki(k0)*
     +      gscomp(k0)*gscomp(6)**2/(6.*gscomp(k0))+
     +      ejecmass*denomki(6)*
     +      gscomp(k0)*gscomp(6)**2/(6.*gscomp(6))
          if (termk0.ne.0.) then
            term1=Ccl*xsreacinc*ejec2sp1*ejecmass
            termki=Ckn*term1*Pn(6)*gscomp(k0)*gscomp(6)/termk0
          endif
        else
          do 50 type2=1,6
            if (type2.ge.3.and.type2.le.5) goto 50
            termin0=projmass*denomki(k0)*
     +        gscomp(k0)*gscomp(type2)**2/(6.*gscomp(k0))+
     +        projmass*denomki(type2)*
     +        gscomp(k0)*gscomp(type2)**2/(6.*gscomp(type2))
            if (termin0.ne.0.) then
              term1=Ccl*xsreacinc*proj2sp1*projmass
              termki=termki+Ckn*term1*Pn(type2)*gscomp(type2)*
     +          gscomp(type2)/termin0
            endif
   50     continue
        endif
c
c Calculation of knockout or inelastic spectra.
c
c factor1  : help variable
c Eres     : total energy of residual system
c Etotal   : total energy of compound system (target + projectile)
c S        : separation energy per particle
c preeqpair: pre-equilibrium pairing energy
c parZ     : charge number of particle
c parN     : neutron number of particle
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
c U        : help variable
c xspreeqki: preequilibrium cross section per particle type and
c            outgoing energy for knockout and inelastic
c
c Knockout term that depends on emission energy.
c
          if (flagknock.or.flaginel) then
            U=max(Eres-AKO,0.)
            xspreeqki(type,nen)=termki*factor1*U
          endif
  110   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
