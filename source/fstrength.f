      function fstrength(Zcomp,Ncomp,Efs,Egamma,irad,l)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 4, 2012
c | Task  : Gamma ray strength functions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      integer      Zcomp,Ncomp,irad,l,i,nen,nT,it,jt,itemp
      real         fstrength,Efs,Egamma,sgr1,egr1,ggr1,kgr1,egr2,ggr2,
     +             Egam2,e,Tnuc,ggredep0,ggredep,enum,denom,factor1,
     +             factor2,eb,ee,Eq(0:numgamqrpa),gamb,game,f1,f2,fb,fe,
     +             tpr1,epr1,gpr1,epr2,gpr2,Tb,Te
c
c ************************* Strength functions *************************
c
c fstrength: gamma-ray strength function
c Ex       : excitation energy
c Egamma   : gamma energy
c irad     : variable to indicate M(=0) or E(=1) radiation
c ngr      : number of GR
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c gamadjust: logical for energy-dependent gamma adjustment
c adjust   : subroutine for energy-dependent parameter adjustment
c factor1,.: multiplication factor
c sgr,sgr1 : strength of GR
c egr,egr1 : energy of GR
c ggr,ggr1 : width of GR
c kgr,kgr1 : constant for gamma-ray strength function
c egr2,ggr2: help variables
c Egam2    : help variable
c strength : strength function of Kopecky-Uhl (1) or Brink-Axel (2)
c
c Models for E1 gamma-ray strength function:
c
c 1. Kopecky-Uhl
c 2. Brink-Axel
c 3. Goriely HFBCS
c 4. Goriely HFB
c 5. Goriely Hybrid model
c
      fstrength=0.
      do 10 i=1,ngr(Zcomp,Ncomp,irad,l)
       if (gamadjust(Zcomp,Ncomp)) then
          key='sgr'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor1)
          key='sgradjust'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor2)
          sgr1=factor1*factor2*sgr(Zcomp,Ncomp,irad,l,i)*
     +      sgradjust(Zcomp,Ncomp,irad,l,i)
          key='ggr'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor1)
          key='ggradjust'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor2)
          ggr1=factor1*factor2*ggr(Zcomp,Ncomp,irad,l,i)*
     +      ggradjust(Zcomp,Ncomp,irad,l,i)
          key='egr'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor1)
          key='egradjust'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor2)
          egr1=factor1*factor2*egr(Zcomp,Ncomp,irad,l,i)*
     +      egradjust(Zcomp,Ncomp,irad,l,i)
        else
          sgr1=sgr(Zcomp,Ncomp,irad,l,i)
          egr1=egr(Zcomp,Ncomp,irad,l,i)
          ggr1=ggr(Zcomp,Ncomp,irad,l,i)
        endif
        kgr1=kgr(Zcomp,Ncomp,irad,l)
        egr2=egr1**2
        ggr2=ggr1**2
        Egam2=Egamma**2
c
c 1. Kopecky-Uhl generalized Lorentzian.
c
c Efs      : fast particle energy for gamma ray strength function
c Tnuc     : nuclear temperature
c k0       : index of incident particle
c Einc     : incident energy
c delta    : energy shift
c S        : separation energy per particle
c alev     : level density parameter
c ggredep0 : energy dependent damping width at zero gamma energy
c ggredep  : energy dependent damping width
c twopi    : 2.*pi
c enum     : enumerator of Lorentzian
c denom    : denominator of Lorentzian
c factor1-2: help variables
c qrpaexist: flag for existence of tabulated QRPA strength functions
c
        Tnuc=0.
        if (strength.eq.1.and.l.eq.1.and.irad.eq.1) then
          if (k0.gt.0.or.Egamma.ne.Einc) then
            e=min(Efs,20.)+S(Zcomp,Ncomp,k0)-delta(Zcomp,Ncomp,0)-Egamma
            if (e.gt.0..and.alev(Zcomp,Ncomp).gt.0.)
     +        Tnuc=sqrt(e/alev(Zcomp,Ncomp))
          endif
          ggredep0=ggr1*twopi**2*Tnuc**2/egr2
          ggredep=ggredep0+ggr1*Egam2/egr2
          enum=ggredep*Egamma
          denom=(Egam2-egr2)**2+Egam2*ggredep**2
          factor1=enum/denom
          factor2=0.7*ggredep0/(egr1**3)
          fstrength=fstrength+kgr1*sgr1*ggr1*(factor1+factor2)
        endif
c
c 2. Brink-Axel standard Lorentzian.
c
        if (strength.eq.2.or.((strength.eq.3.or.strength.eq.4).and.
     +    .not.qrpaexist(Zcomp,Ncomp)).or.l.ne.1.or.irad.ne.1) then
          enum=0.
          if (Egamma.gt.0.001) then
            enum=ggr2*Egamma**(3-2*l)
            denom=(Egam2-egr2)**2+Egam2*ggr2
            fstrength=fstrength+kgr1*sgr1*enum/denom
          endif
        endif
c
c 3+4. Tabulated QRPA strength functions
c
c locate    : subroutine to find value in ordered table
c numgamqrpa: number of energies for QRPA strength function
c eqrpa,Eq  : energy grid for QRPA strength function
c fe1qrpa   : tabulated QRPA strength function
c eb,ee,....: help variables
c
        if ((strength.eq.3.or.strength.eq.4).and.qrpaexist(Zcomp,Ncomp)
     +    .and.l.eq.1.and.irad.eq.1) then
          if ((k0.gt.0.or.Egamma.ne.Einc).and.nTqrpa.gt.1) then
            e=min(Efs,20.)+S(Zcomp,Ncomp,k0)-delta(Zcomp,Ncomp,0)-Egamma
            if (e.gt.0..and.alev(Zcomp,Ncomp).gt.0.)
     +        Tnuc=sqrt(e/alev(Zcomp,Ncomp))
            nT=nTqrpa
            do 110 it=1,nTqrpa
              if (Tqrpa(it).gt.Tnuc) then
                nT=it-1
                goto 120
              endif
  110       continue
  120       Tb=Tqrpa(nT)
            if (nT.lt.nTqrpa) then
              Te=Tqrpa(nT+1)
            else
              Te=Tb
            endif
            itemp=2
          else
            Tb=0.
            Te=0.
            itemp=1
            nT=1
          endif
          do 130 nen=0,numgamqrpa
            Eq(nen)=eqrpa(Zcomp,Ncomp,nen)
  130     continue
          do 140 it=1,itemp
            jt=nT
            if (it.eq.2) jt=nT+1
            if (jt.gt.nTqrpa) jt=nTqrpa
            if (Egamma.le.Eq(numgamqrpa)) then
              call locate(Eq,0,numgamqrpa,Egamma,nen)
              eb=Eq(nen)
              ee=Eq(nen+1)
              gamb=fe1qrpa(Zcomp,Ncomp,nen,jt)
              game=fe1qrpa(Zcomp,Ncomp,nen+1,jt)
              if (gamb.gt.0..and.game.gt.0.) then
                f1=log10(gamb)+(Egamma-eb)/(ee-eb)*
     +            (log10(game)-log10(gamb))
                f2=10.**f1
              else
                f2=gamb+(Egamma-eb)/(ee-eb)*(game-gamb)
              endif
            else
              eb=Eq(numgamqrpa-1)
              ee=Eq(numgamqrpa)
              gamb=fe1qrpa(Zcomp,Ncomp,numgamqrpa-1,jt)
              game=fe1qrpa(Zcomp,Ncomp,numgamqrpa,jt)
              if (gamb.gt.0..and.game.gt.0.) then
                f1=log10(gamb)+(Egamma-eb)/(ee-eb)*
     +            (log10(game)-log10(gamb))
                f2=10.**f1
              else
                f2=gamb+(Egamma-eb)/(ee-eb)*(game-gamb)
              endif
            endif
            if (it.eq.1) fb=f2
            if (it.eq.2) fe=f2
  140     continue
          if (nTqrpa.gt.1.and.Tb.ne.Te) then
            if (fb.gt.0..and.fe.gt.0.) then
              f1=log10(fb)+(Tnuc-Tb)/(Te-Tb)*
     +          (log10(fe)-log10(fb))
              f2=10.**f1
            else
              f2=fb+(Tnuc-Tb)/(Te-Tb)*(fe-fb)
            endif
          endif
          fstrength=f2
        endif
c
c 5. Goriely Hybrid model
c
        if (strength.eq.5.and.l.eq.1.and.irad.eq.1) then
          if (k0.gt.0.or.Egamma.ne.Einc) then
            e=min(Efs,20.)+S(Zcomp,Ncomp,1)-delta(Zcomp,Ncomp,0)-Egamma
            if (e.gt.0..and.alev(Zcomp,Ncomp).gt.0.)
     +        Tnuc=sqrt(e/alev(Zcomp,Ncomp))
          endif
          if (Egamma.gt.0.) then
            ggredep=0.7*ggr1*(Egamma/egr1+twopi**2*Tnuc**2/Egamma/egr1)
            enum=ggredep*Egamma
            denom=(Egam2-egr2)**2+Egam2*ggredep*ggr1
            factor1=enum/denom
            fstrength=fstrength+kgr1*sgr1*ggr1*factor1
          endif
        endif
   10 continue
c
c Inclusion of additional extra strength (Pygmy Resonance),
c only if explicitly specified in the input
c
      tpr1=tpr(Zcomp,Ncomp,irad,l,1)
      if (Egamma.gt.0.001.and.tpr1.gt.0.) then
        if (gamadjust(Zcomp,Ncomp)) then
          key='tpr'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor1)
          key='tpradjust'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor2)
          tpr1=factor1*factor2*tpr(Zcomp,Ncomp,irad,l,1)*
     +      tpradjust(Zcomp,Ncomp,irad,l,1)
          key='gpr'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor1)
          key='gpradjust'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor2)
          gpr1=factor1*factor2*gpr(Zcomp,Ncomp,irad,l,1)*
     +      gpradjust(Zcomp,Ncomp,irad,l,1)
          key='epr'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor1)
          key='epradjust'
          call adjust(Egamma,key,Zcomp,Ncomp,0,l,factor2)
          epr1=factor1*factor2*epr(Zcomp,Ncomp,irad,l,1)*
     +      epradjust(Zcomp,Ncomp,irad,l,1)
        else
          epr1=epr(Zcomp,Ncomp,irad,l,1)
          gpr1=gpr(Zcomp,Ncomp,irad,l,1)
        endif
        kgr1=kgr(Zcomp,Ncomp,irad,l)
        epr2=epr1**2
        gpr2=gpr1**2
        Egam2=Egamma**2
        enum=gpr2*Egamma**(3-2*l)
        denom=(Egam2-epr2)**2+Egam2*gpr2
        fstrength=fstrength+kgr1*tpr1*enum/denom
      endif
c
c Reduction of gamma-strength for isospin forbidden transitions into
c Z=N nuclei
c
c fiso: correction factor for isospin forbidden transitions
c
      fstrength=fstrength/fiso(Zcomp,Ncomp,k0)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
