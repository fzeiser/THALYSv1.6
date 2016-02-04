      subroutine excitoninit
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : July 9, 2004
c | Task  : Initialization of exciton model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,type,Zix,Nix,zproj,nproj,zejec,nejec,aejec,p,
     +        ppi,hpi,npi,pnu,hnu,nnu
      real    ZoverA,NoverA,Qenum,Qdenom,factor1,factor2,factor3,
     +        factor4,factor5,factor6,factor7
c
c Factors for the emission rate.
c
c wfac      : factor for emission rate
c pi2h3c2   : 1/(pi*pi*clight*clight*hbar**3) in mb**-1.MeV**-3.s**-1
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c parskip   : logical to skip outgoing particle
c Zindex    : charge number index for residual nucleus
c Nindex    : neutron number index for residual nucleus
c amupi2h3c2: amu/(pi*pi*clight*clight*hbar**3) in mb**-1.MeV**-2.s**-1
c parspin   : spin of particle
c redumass  : reduced mass
c
      wfac(0)=pi2h3c2
      Zcomp=0
      Ncomp=0
      do 10 type=1,6
        if (parskip(type)) goto 10
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        wfac(type)=amupi2h3c2*(2.*parspin(type)+1.)*
     +    redumass(Zix,Nix,type)
   10 continue
c
c ******************* Q-factors for emission rate **********************
c
c flag2comp       : flag for two-component pre-equilibrium model
c ZoverA,NoverA   : help variables
c Zinit           : charge number of initial compound nucleus
c Ninit           : neutron number of initial compound nucleus
c Ainit           : mass number of initial compound nucleus
c zproj,zejec,parZ: charge number of particle
c nproj,nejec,parN: neutron number of particle
c k0              : index of incident particle
c aejec,parA      : mass number of particle
c numparx         : maximal particle number
c Qfactor         : Q-factor for neutron/proton distinction
c p               : particle number
c maxpar          : maximal particle number
c Qenum,Qdenom    : help variables
c ppi             : proton particle number
c hpi             : proton hole number
c npi             : proton exciton number
c pnu             : neutron particle number
c hnu             : neutron hole number
c nnu             : neutron exciton number
c factor1-7       : help variables
c nfac            : n!
c
      if (flag2comp) return
      ZoverA=real(Zinit)/real(Ainit)
      NoverA=real(Ninit)/real(Ainit)
      zproj=parZ(k0)
      nproj=parN(k0)
      do 110 type=1,6
        do 120 p=1,numparx
          Qfactor(type,p)=1.
  120   continue
        if (parskip(type)) goto 110
        zejec=parZ(type)
        nejec=parN(type)
        aejec=parA(type)
        do 130 p=1,maxpar
          if (p.lt.aejec) goto 130
          Qenum=0.
          Qdenom=0.
          do 140 ppi=zproj,p-nproj
            hpi=ppi-zproj
            npi=ppi+hpi
            pnu=p-ppi
            hnu=p-nproj-ppi
            nnu=pnu+hnu
            factor1=ZoverA**npi
            factor2=NoverA**nnu
            factor3=1./real(nfac(ppi)*nfac(pnu))
            factor4=real(nfac(hpi)*nfac(hnu))
            Qdenom=Qdenom+factor1*factor2*factor3/factor4
            if (ppi.lt.zejec.or.pnu.lt.nejec) goto 140
            factor5=ZoverA**(npi-zejec)
            factor6=NoverA**(nnu-nejec)
            factor7=1./real(nfac(ppi-zejec)*nfac(pnu-nejec))
            Qenum=Qenum+factor5*factor6*factor7/factor4
  140     continue
          if (Qdenom.gt.0.) Qfactor(type,p)=(Qenum/Qdenom)*
     +      real(nfac(p-aejec))/real(nfac(p))
  130   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
