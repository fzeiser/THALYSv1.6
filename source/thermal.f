      subroutine thermal
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 22, 2012
c | Task  : Estimate of thermal cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      Zcomp,Ncomp,type,nen,idc,nex,i1,i2
      real         Etherm,ctherm,xscapres,xscap1,ratio,ratiores,ratiop,
     +             ratioresp,xspres,xsp1,ratioalpha,ratioresalpha,
     +             xsalphares,xsalpha1,Ereslog,elog,ealog,Eratio,xsa,R,
     +             Rres,xsres,xsalog,xsreslog,xs
c
c *********************** Extrapolate cross sections *******************
c
c For non-threshold channels, the cross sections are extrapolated
c down to 1.e-5 eV. Capture values at thermal energies are used.
c For energies up to 1 eV, the 1/sqrt(E) law is used. Between 1 eV
c and the first energy at which TALYS performs the statistical model
c calculation, we use logarithmic interpolation.
c
c Zcomp        : charge number index for compound nucleus
c Ncomp        : neutron number index for compound nucleus
c Etherm       : thermal energy
c ctherm       : constant for 1/v capture cross section function
c xscapres     : capture cross section in resonance region
c E1v          : energy at end of 1/v region
c xscap1       : capture cross section at first incident energy
c xspopnuc     : population cross section per nucleus
c ratio        : ratio thermal/first energy for gamma capture
c xscaptherm   : thermal capture cross section
c xsptherm     : thermal (n,p) cross section
c xsalphatherm : thermal (n,a) cross section
c ratiores     : ratio start of resonance region/first energy
c ratiop       : ratio thermal/first energy for proton
c ratioresp    : ratio start of resonance region/first energy
c ratioalpha   : ratio thermal/first energy for proton
c ratioresalpha: ratio start of resonance region/first energy
c
      Zcomp=0
      Ncomp=0
      Etherm=2.53e-8
      ctherm=sqrt(Etherm)*xscaptherm
      xscapres=ctherm/sqrt(E1v)
      xscap1=xspopnuc(Zcomp,Ncomp)
      if (xscap1.gt.0.) then
        ratio=xscaptherm/xscap1
        ratiores=xscapres/xscap1
      else
        ratio=1.
        ratiores=1.
      endif
c
c Protons
c
      ratiop=ratio
      ratioresp=ratiores
      if (xsptherm.ne.0.) then
        ctherm=sqrt(Etherm)*xsptherm
        xspres=ctherm/sqrt(E1v)
        xsp1=xspopnuc(1,0)
        if (xsp1.gt.0.) then
          ratiop=xsptherm/xsp1
          ratioresp=xspres/xsp1
        else
          ratiop=1.
        endif
      endif
c
c Alpha particles
c
      ratioalpha=ratio
      ratioresalpha=ratiores
      if (xsalphatherm.ne.0.) then
        ctherm=sqrt(Etherm)*xsalphatherm
        xsalphares=ctherm/sqrt(E1v)
        xsalpha1=xspopnuc(2,2)
        if (xsalpha1.gt.0.) then
          ratioalpha=xsalphatherm/xsalpha1
          ratioresalpha=xsalphares/xsalpha1
        else
          ratioalpha=1.
        endif
      endif
c
c Determine cross sections on low-energy grid
c
c Ereslog      : logarithm of energy at start of resonance region
c numinclow    : number of incident energies below Elow
c elog,ealog   : help variables
c Eratio       : energy ratio
c eninc        : incident energy in MeV
c fxsnonel     : non-elastic cross section
c fxselastot   : total elastic cross section (shape + compound)
c fxstotinc    : total cross section (neutrons only) for incident
c                channel
c fxscompel    : compound elastic cross section
c fxselasinc   : total elastic cross section (neutrons only) for
c                incident channel
c fxsreacinc   : reaction cross section for incident channel
c fxscompnonel : total compound non-elastic cross section
c fxsdirdiscsum: total direct cross section
c fxspreeqsum  : total preequilibrium cross section summed over
c                particles
c
      Ereslog=log(E1v)
      do 110 nen=1,numinclow
        elog=log(eninc(nen))
        ealog=log(eninc(numinclow+1))
        Eratio=sqrt(Etherm)/sqrt(eninc(nen))
        fxsnonel(nen)=0.
        fxselastot(nen)=0.
        fxstotinc(nen)=0.
        fxscompel(nen)=0.
        fxselasinc(nen)=0.
        fxsreacinc(nen)=0.
        fxscompnonel(nen)=0.
        fxsdirdiscsum(nen)=0.
        fxspreeqsum(nen)=0.
c
c Exclusive channel cross sections
c
c flagchannels : flag for exclusive channels calculation
c idnum        : counter for exclusive channel
c fxschannel   : channel cross section
c fxsgamchannel: gamma channel cross section
c fxsgamdischan: discrete gamma channel cross section
c numlev       : maximum number of included discrete levels
c Ethrexc      : threshold incident energy for exclusive channel
c idchannel    : identifier for exclusive channel
c xseps        : limit for cross sections
c xsalog       : help variable
c xsreslog     : cross section at start of resonance region
c pol1         : subroutine for interpolation of first order
c fxsratio     : ratio of exclusive cross section over residual
c                production cross section (for exclusive gamma ray
c                intensities)
c Nlast        : last discrete level
c fxschaniso   : channel cross section per isomer
c fexclyield   : exclusive channel yield per isomer
c
        if (flagchannels) then
          do 120 idc=0,idnum
            fxschannel(nen,idc)=0.
            fxsgamchannel(nen,idc)=0.
            do 130 i1=1,numlev
              do 130 i2=0,i1
                fxsgamdischan(nen,idc,i1,i2)=0.
  130       continue
            if (eninc(nen).le.Ethrexcl(idc,0)) goto 120
            if (idchannel(idc).eq.100000) goto 120
            xsa=xschannel(idc)
            if (xsa.lt.xseps) goto 120
            R=ratio
            Rres=ratiores
            if (idchannel(idc).eq.010000) then
              R=ratiop
              Rres=ratioresp
            endif
            if (idchannel(idc).eq.000001) then
              R=ratioalpha
              Rres=ratioresalpha
            endif
            if (eninc(nen).gt.E1v) then
              xsres=xsa*Rres
              if (xsres.le.0..or.xsa.le.0.) goto 120
              xsalog=log(xsa)
              xsreslog=log(xsres)
              call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
              fxschannel(nen,idc)=exp(xs)
            else
              xs=xsa*R*Eratio
              fxschannel(nen,idc)=xs
            endif
            fxsnonel(nen)=fxsnonel(nen)+fxschannel(nen,idc)
            fxsratio(nen,idc)=xsratio(idc)
            do 140 nex=0,Nlast(0,0,0)
              fxschaniso(nen,idc,nex)=0.
              if (eninc(nen).le.Ethrexcl(idc,nex)) goto 140
              xsa=xschaniso(idc,nex)
              if (xsa.lt.xseps) goto 140
              if (eninc(nen).gt.E1v) then
                xsres=xsa*Rres
                if (xsres.le.0..or.xsa.le.0.) goto 140
                xsalog=log(xsa)
                xsreslog=log(xsres)
                call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                fxschaniso(nen,idc,nex)=exp(xs)
              else
                xs=xsa*R*Eratio
                fxschaniso(nen,idc,nex)=xs
              endif
              fexclyield(nen,idc,nex)=exclyield(idc,nex)
  140       continue
            xsa=xsgamchannel(idc)
            if (xsa.lt.xseps) goto 120
            if (eninc(nen).gt.E1v) then
              xsres=xsa*Rres
              if (xsres.le.0..or.xsa.le.0.) goto 120
              xsalog=log(xsa)
              xsreslog=log(xsres)
              call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
              fxsgamchannel(nen,idc)=exp(xs)
            else
              xs=xsa*R*Eratio
              fxsgamchannel(nen,idc)=xs
            endif
            do 150 i1=1,numlev
              do 150 i2=0,i1
                xsa=xsgamdischan(idc,i1,i2)
                if (xsa.lt.xseps) goto 150
                if (eninc(nen).gt.E1v) then
                  xsres=xsa*Rres
                  if (xsres.le.0..or.xsa.le.0.) goto 150
                  xsalog=log(xsa)
                  xsreslog=log(xsres)
                  call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                  fxsgamdischan(nen,idc,i1,i2)=exp(xs)
                else
                  xs=xsa*R*Eratio
                  fxsgamdischan(nen,idc,i1,i2)=xs
                endif
  150       continue
  120     continue
        endif
c
c Binary cross sections
c
c fxsbinary: cross section from initial compound to residual nucleus
c Ethresh  : threshold incident energy for residual nucleus
c parZ     : charge number of particle
c parN     : neutron number of particle
c k0       : index of incident particle
c
        do 210 type=0,6
          fxsbinary(nen,type)=0.
          if (eninc(nen).le.Ethresh(parZ(type),parN(type),0)) goto 210
          if (type.eq.k0) goto 210
          xsa=xsbinary(type)
          if (xsa.lt.xseps) goto 210
          R=ratio
          Rres=ratiores
          if (type.eq.2) then
            R=ratiop
            Rres=ratioresp
          endif
          if (type.eq.6) then
            R=ratioalpha
            Rres=ratioresalpha
          endif
          if (eninc(nen).gt.E1v) then
            xsres=xsa*Rres
            if (xsres.le.0..or.xsa.le.0.) goto 210
            xsalog=log(xsa)
            xsreslog=log(xsres)
            call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
            fxsbinary(nen,type)=exp(xs)
          else
            xs=xsa*R*Eratio
            fxsbinary(nen,type)=xs
          endif
  210   continue
c
c Residual production cross sections
c
c maxZ     : maximal number of protons away from the initial compound
c            nucleus
c maxN     : maximal number of neutrons away from the initial compound
c            nucleus
c fxspopnuc: population cross section per nucleus
c fxspopex : population cross section summed over spin and parity
c fxsbranch: branching ratio for isomeric cross section
c flagastro: flag for calculation of astrophysics reaction rate
c
        do 310 Zcomp=0,maxZ
          do 310 Ncomp=0,maxN
            fxspopnuc(nen,Zcomp,Ncomp)=0.
            if (eninc(nen).le.Ethresh(Zcomp,Ncomp,0)) goto 310
            if (Zcomp.eq.parZ(k0).and.Ncomp.eq.parN(k0)) goto 310
            xsa=xspopnuc(Zcomp,Ncomp)
            if (xsa.lt.xseps) goto 310
            R=ratio
            Rres=ratiores
            if (Zcomp.eq.1.and.Ncomp.eq.0) then
              R=ratiop
              Rres=ratioresp
            endif
            if (Zcomp.eq.2.and.Ncomp.eq.2) then
              R=ratioalpha
              Rres=ratioresalpha
            endif
            if (eninc(nen).gt.E1v) then
              xsres=xsa*Rres
              if (xsres.le.0..or.xsa.le.0.) goto 310
              xsalog=log(xsa)
              xsreslog=log(xsres)
              call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
              fxspopnuc(nen,Zcomp,Ncomp)=exp(xs)
            else
              xs=xsa*R*Eratio
              fxspopnuc(nen,Zcomp,Ncomp)=xs
            endif
            if (.not.flagchannels) fxsnonel(nen)=fxsnonel(nen)+
     +        fxspopnuc(nen,Zcomp,Ncomp)
            do 320 nex=0,Nlast(Zcomp,Ncomp,0)
              fxspopex(nen,Zcomp,Ncomp,nex)=0.
              if (eninc(nen).le.Ethresh(Zcomp,Ncomp,nex)) goto 320
              xsa=xspopex(Zcomp,Ncomp,nex)
              if (xsa.lt.xseps) goto 320
              if (eninc(nen).gt.E1v) then
                xsres=xsa*Rres
                if (xsres.le.0..or.xsa.le.0.) goto 320
                xsalog=log(xsa)
                xsreslog=log(xsres)
                call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                fxspopex(nen,Zcomp,Ncomp,nex)=exp(xs)
              else
                xs=xsa*R*Eratio
                fxspopex(nen,Zcomp,Ncomp,nex)=xs
              endif
              fxsbranch(nen,Zcomp,Ncomp,nex)=xsbranch(Zcomp,Ncomp,nex)
  320       continue
            if (flagastro) xsastro(Zcomp,Ncomp,nen)=
     +        fxspopnuc(nen,Zcomp,Ncomp)
  310   continue
c
c Reactions to discrete states
c
c fxsexclusive: exclusive single channel cross section
c fxsdisctot  : total cross section summed over discrete states
c fxsexclcont : exclusive single channel cross section for continuum
c fxsngn      : total (projectile,gamma-ejectile) cross section
c fxsdisc     : total cross section for discrete state
c fxsdirdisc  : direct cross section for discrete state
c fxscompdisc : compound cross section for discrete state
c
        do 410 type=0,6
          fxsexclusive(nen,type)=0.
          fxsdisctot(nen,type)=0.
          fxsexclcont(nen,type)=0.
          fxsngn(nen,type)=0.
          if (eninc(nen).le.Ethresh(parZ(type),parN(type),0)) goto 410
          if (type.eq.k0) goto 410
          R=ratio
          Rres=ratiores
          if (type.eq.2) then
            R=ratiop
            Rres=ratioresp
          endif
          if (type.eq.6) then
            R=ratioalpha
            Rres=ratioresalpha
          endif
          xsa=xsexclusive(type)
          if (xsa.ge.xseps) then
            if (eninc(nen).gt.E1v) then
              xsres=xsa*Rres
              if (xsres.gt.0..or.xsa.gt.0.) then
                xsalog=log(xsa)
                xsreslog=log(xsres)
                call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                fxsexclusive(nen,type)=exp(xs)
              endif
            else
              xs=xsa*R*Eratio
              fxsexclusive(nen,type)=xs
            endif
          endif
          xsa=xsdisctot(type)
          if (xsa.ge.xseps) then
            if (eninc(nen).gt.E1v) then
              xsres=xsa*Rres
              if (xsres.gt.0..or.xsa.gt.0.) then
                xsalog=log(xsa)
                xsreslog=log(xsres)
                call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                fxsdisctot(nen,type)=exp(xs)
              endif
            else
              xs=xsa*R*Eratio
              fxsdisctot(nen,type)=xs
            endif
          endif
          xsa=xsexclcont(type)
          if (xsa.ge.xseps) then
            if (eninc(nen).gt.E1v) then
              xsres=xsa*Rres
              if (xsres.gt.0..or.xsa.gt.0.) then
                xsalog=log(xsa)
                xsreslog=log(xsres)
                call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                fxsexclcont(nen,type)=exp(xs)
              endif
            else
              xs=xsa*R*Eratio
              fxsexclcont(nen,type)=xs
            endif
          endif
          xsa=xsngn(type)
          if (xsa.ge.xseps) then
            if (eninc(nen).gt.E1v) then
              xsres=xsa*Rres
              if (xsres.gt.0..or.xsa.gt.0.) then
                xsalog=log(xsa)
                xsreslog=log(xsres)
                call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                fxsngn(nen,type)=exp(xs)
              endif
            else
              xs=xsa*R*Eratio
              fxsngn(nen,type)=xs
            endif
          endif
          do 420 nex=0,Nlast(parZ(type),parN(type),0)
            fxsdisc(nen,type,nex)=0.
            fxsdirdisc(nen,type,nex)=0.
            fxscompdisc(nen,type,nex)=0.
            xsa=xsdisc(type,nex)
            if (xsa.ge.xseps) then
              if (eninc(nen).gt.E1v) then
                xsres=xsa*Rres
                if (xsres.gt.0..or.xsa.gt.0.) then
                  xsalog=log(xsa)
                  xsreslog=log(xsres)
                  call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                  fxsdisc(nen,type,nex)=exp(xs)
                endif
              else
                xs=xsa*R*Eratio
                fxsdisc(nen,type,nex)=xs
              endif
            endif
            xsa=xsdirdisc(type,nex)
            if (xsa.ge.xseps) then
              if (eninc(nen).gt.E1v) then
                xsres=xsa*Rres
                if (xsres.gt.0..or.xsa.gt.0.) then
                  xsalog=log(xsa)
                  xsreslog=log(xsres)
                  call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                  fxsdirdisc(nen,type,nex)=exp(xs)
                endif
              else
                xs=xsa*R*Eratio
                fxsdirdisc(nen,type,nex)=xs
              endif
            endif
            xsa=xscompdisc(type,nex)
            if (xsa.ge.xseps) then
              if (eninc(nen).gt.E1v) then
                xsres=xsa*Rres
                if (xsres.gt.0..or.xsa.gt.0.) then
                  xsalog=log(xsa)
                  xsreslog=log(xsres)
                  call pol1(Ereslog,ealog,xsreslog,xsalog,elog,xs)
                  fxscompdisc(nen,type,nex)=exp(xs)
                endif
              else
                xs=xsa*R*Eratio
                fxscompdisc(nen,type,nex)=xs
              endif
            endif
  420     continue
  410   continue
c
c Total cross sections
c
        fxselastot(nen)=xselastot
        fxstotinc(nen)=fxselastot(nen)+fxsnonel(nen)
        fxscompel(nen)=xscompel
        fxselasinc(nen)=xselasinc
        fxsreacinc(nen)=fxsnonel(nen)+fxscompel(nen)
        fxscompnonel(nen)=fxsnonel(nen)
        fxsdirdiscsum(nen)=xsdirdiscsum
        fxspreeqsum(nen)=xspreeqsum
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
