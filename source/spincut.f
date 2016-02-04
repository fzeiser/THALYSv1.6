      function spincut(Zix,Nix,ald,Eex,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : April 4, 2012
c | Task  : Spin cutoff factor
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      character*80 key
      integer      Zix,Nix,ibar,ldmod
      real         spincut,ald,Eex,factor,Rs,s2,scutconst,Em,aldm,
     +             ignatyuk,Umatch,s2m,s2d,Ed,U
c
c *********************** Spin cutoff parameter ************************
c
c spincut     : spin cutoff factor
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c ald         : level density parameter
c Eex         : excitation energy
c ibar        : fission barrier
c
c Models for spin cutoff factor:
c
c spincutmodel 1: s.c. = I0 * (a/alimit) * sqrt (U/a) =
c                        I0/alimit * sqrt (a.U)
c spincutmodel 2: s.c. = I0 * sqrt (U/a)
c                 where I0 = 0.01389 * A^(5/3)
c
c We interpolate between the discrete level spin cutoff factor and
c the value at the matching energy.
c
c 1. Below the matching energy
c
c adjust      : subroutine for energy-dependent parameter adjustment
c factor      : multiplication factor
c ldadjust    : logical for energy-dependent level density adjustment
c ldmodel     : level density model
c spincutmodel: model for spin cutoff factor for ground state
c scutconst   : constant for spin cutoff factor
c Rspincut    : adjustable constant (global) for spin cutoff factor
c s2adjust    : adjustable constant (Z,A,barrier-dependent) for spin
c               cutoff parameter
c Irigid0     : undeformed rigid body value of moment of inertia
c alimit      : asymptotic level density parameter
c Em          : matching energy
c Exmatch     : matching point for Ex
c S           : separation energy per particle
c Ucrit       : critical U
c pair        : pairing energy
c Pshift      : adjustable pairing shift
c aldm        : level density parameter at matching energy
c ignatyuk    : function for energy dependent level density parameter a
c Umatch      : Exmatch - pairing
c aldcrit     : critical level density parameter
c Tcrit       : critical temperature
c Ediscrete   : energy of middle of discrete level region
c s2m,s2d,Ed  : help variables
c scutoffdisc : spin cutoff factor for discrete level region
c flagcolldamp: flag for damping of collective effects in effective
c               level density (without explicit collective enhancement)
c               Only used for Bruyeres-le-Chatel (Pascal Romain) fission
c               model
c beta2       : deformation parameter
c
      key='rspincut'
      call adjust(Eex,key,0,0,0,0,factor)
      Rs=factor*Rspincut
      if (ldadjust(Zix,Nix)) then
        key='s2adjust'
        call adjust(Eex,key,Zix,Nix,0,ibar,factor)
        s2=factor*s2adjust(Zix,Nix,ibar)
      else
        s2=s2adjust(Zix,Nix,ibar)
      endif
      ldmod=ldmodel(Zix,Nix)
      spincut=1.
      if (spincutmodel.eq.1) then
        scutconst=Rs*s2*Irigid0(Zix,Nix)/alimit(Zix,Nix)
      else
        scutconst=Rs*s2*Irigid0(Zix,Nix)
      endif
      Em=Exmatch(Zix,Nix,ibar)
      if (ldmod.eq.2.or.ldmod.ge.4) Em=S(Zix,Nix,1)
      if (ldmod.eq.3)
     +  Em=Ucrit(Zix,Nix)-pair(Zix,Nix)-Pshift(Zix,Nix,ibar)
      if (Eex.le.Em) then
        aldm=ignatyuk(Zix,Nix,Em,ibar)
        Umatch=Em-pair(Zix,Nix)-Pshift(Zix,Nix,ibar)
        if (Umatch.gt.0.) then
          if (spincutmodel.eq.1) then
            if (ldmod.eq.3) then
              s2m=scutconst*aldcrit(Zix,Nix)*Tcrit(Zix,Nix)
            else
              s2m=scutconst*sqrt(aldm*Umatch)
            endif
          else
            if (ldmod.eq.3) then
              s2m=scutconst*Tcrit(Zix,Nix)
            else
              s2m=scutconst*sqrt(Umatch/aldm)
            endif
          endif
        else
          s2m=scutoffdisc(Zix,Nix,ibar)
        endif
        if (flagcolldamp.and.ibar.ne.0)
     +    s2m=(1.+beta2(Zix,Nix,ibar)/3.)*s2m
        s2d=scutoffdisc(Zix,Nix,ibar)
        Ed=Ediscrete(Zix,Nix,ibar)
        spincut=s2d
        if (Em.ne.Ed.and.Eex.gt.Ed)
     +    spincut=s2d+(Eex-Ed)/(Em-Ed)*(s2m-s2d)
      else
c
c 2. Above the matching energy
c
c U       : excitation energy minus pairing energy
c delta   : energy shift
c
        U=Eex-delta(Zix,Nix,ibar)
        if (U.gt.0.) then
          if (spincutmodel.eq.1) then
            spincut=scutconst*sqrt(ald*U)
          else
            spincut=scutconst*sqrt(U/ald)
          endif
        else
          spincut=scutoffdisc(Zix,Nix,ibar)
        endif
      endif
      if (flagcolldamp.and.ibar.ne.0)
     +  spincut=(1.+beta2(Zix,Nix,ibar)/3.)*spincut
      spincut=max(scutoffdisc(Zix,Nix,ibar),spincut)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
