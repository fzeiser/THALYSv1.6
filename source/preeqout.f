      subroutine preeqout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : March 28, 2010
c | Task  : Output of pre-equilibrium cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell
      integer Zix,Nix,h,k,nen,n,J,type,p
      real    Eex,gs,ignatyuk,phdens,gsp,gsn,damp,preeqpair,phdens2,
     +        nonpski
c
c ************************ Pre-equilibrium *****************************
c
c 1. Output of particle-hole state densities
c
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c surfwell  : flag for surface effects in finite well
c flag2comp : flag for two-component pre-equilibrium model
c h         : hole number
c Etotal    : total energy of compound system (target + projectile)
c Eex       : excitation energy
c gs,g      : single-particle level density parameter
c flaggshell: flag for energy dependence of single particle level
c             density parameter g
c ignatyuk  : function for energy dependent level density parameter a
c alev      : level density parameter
c preeqpair : pre-equilibrium pairing energy
c pairmodel : model for preequilibrium pairing energy
c phdens    : particle-hole state density
c Efermi    : depth of Fermi well
c gp,gsp    : single-particle proton level density parameter
c gn,gsn    : single-particle neutron level density parameter
c damp      : shell damping factor
c phdens2   : function for two-component particle-hole state density
c n         : exciton number
c maxexc    : maximal exciton number
c RnJ       : spin distribution for particle-hole states
c RnJsum    : (2J+1)*sum over spin distributions
c Esurf     : well depth for surface interaction
c
      Zix=0
      Nix=0
      surfwell=.false.
      write(*,'(/" ++++++++++ PARTIAL STATE DENSITIES ++++++++++")')
      if (.not.flag2comp) then
        write(*,'(/" Particle-hole state densities"/)')
        write(*,'("    Ex   P(n=3)    gs    ",8(i1,"p",i1,"h",6x)/)')
     +    ((h+k,h,k=0,1),h=1,4)
        do 10 nen=1,int(Etotal)
          Eex=real(nen)
          gs=g(0,0)
          if (flaggshell) gs=gs*ignatyuk(Zix,Nix,Eex,0)/alev(0,0)
          write(*,'(1x,3f7.3,1p,8e10.3)') Eex,
     +      preeqpair(Zix,Nix,3,Eex,pairmodel),gs,
     +      ((phdens(Zix,Nix,h+k,h,gs,Eex,Efermi,surfwell),k=0,1),h=1,4)
   10   continue
      else
        write(*,'(/" Particle-hole state densities",/)')
        write(*,'("    Ex   P(n=3)   gp     gn  ",26x,
     +    "Configuration p(p) h(p) p(n) h(n)")')
        write(*,'(28x,9(2x,4i2)/)') 1,1,0,0, 0,0,1,1, 1,1,1,0, 1,0,1,1,
     +    2,1,0,0, 0,0,2,1, 2,2,0,0, 0,0,2,2, 1,1,1,1
        do 20 nen=1,int(Etotal)
          Eex=real(nen)
          gsp=gp(0,0)
          gsn=gn(0,0)
          if (flaggshell) then
            damp=ignatyuk(Zix,Nix,Eex,0)/alev(0,0)
            gsp=gsp*damp
            gsn=gsn*damp
          endif
          write(*,'(1x,4f7.3,1p,9e10.3)') Eex,
     +      preeqpair(Zix,Nix,3,Eex,pairmodel),gsp,gsn,
     +      phdens2(Zix,Nix,1,1,0,0,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,0,0,1,1,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,1,1,1,0,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,1,0,1,1,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,2,1,0,0,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,0,0,2,1,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,2,2,0,0,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,0,0,2,2,gsp,gsn,Eex,Efermi,surfwell),
     +      phdens2(Zix,Nix,1,1,1,1,gsp,gsn,Eex,Efermi,surfwell)
   20   continue
      endif
      write(*,'(/" Particle-hole spin distributions"/)')
      write(*,'("   n    ",9(" J=",i2,"       ")," Sum"/)') (J,J=0,8)
      do 30 n=1,maxexc
        write(*,'(1x,i3,1p,10e12.4)') n,(RnJ(n,J),J=0,8),RnJsum(n)
   30 continue
      write(*,'(/" Effective well depth for surface interaction:",f12.5,
     +  " MeV")') Esurf
c
c 2. Output of pre-equilibrium cross sections
c
c preeqnorm   : preequilibrium normalization factor
c parskip     : logical to skip outgoing particle
c ebegin      : first energy point of energy grid
c eend        : last energy point of energy grid
c parname     : name of particle
c p           : particle number
c maxpar      : maximal particle number
c nonpski     : preequilibrium cross section without pickup etc.
c xspreeq     : preequilibrium cross section per particle type and
c               outgoing energy
c xspreeqps   : preequilibrium cross section per particle type and
c               outgoing energy for pickup and stripping
c xspreeqki   : preequilibrium cross section per particle type and
c               outgoing energy for knockout and inelastic
c xspreeqbu   : preequilibrium cross section per particle type and
c               outgoing energy for breakup
c egrid       : outgoing energy grid
c xsstep      : preequilibrium cross section per particle type, stage
c               and outgoing energy
c xspreeqtot  : preequilibrium cross section per particle type
c xspreeqtotps: preequilibrium cross section per particle type for
c               pickup and stripping
c xspreeqtotki: preequilibrium cross section per particle type for
c               knockout and inelastic
c xspreeqtotbu: preequilibrium cross section per particle type for
c               breakup
c xssteptot   : preequilibrium cross section per particle type and stage
c xspreeqsum  : total preequilibrium cross section summed over particles
c
      write(*,'(/" ++++++++++ TOTAL PRE-EQUILIBRIUM CROSS SECTIONS",
     +  " ++++++++++")')
      if (preeqnorm.ne.0.)
     +  write(*,'(/" Pre-equilibrium normalization factor: ",f8.5/)')
     +  preeqnorm
      do 110 type=0,6
        if (parskip(type)) goto 110
        if (ebegin(type).ge.eend(type)) goto 110
        write(*,'(/" Pre-equilibrium cross sections for ",a8/)')
     +    parname(type)
        write(*,'("    E     Total",6("       p=",i1),
     +    "     Total  Pickup/Strip Knockout Breakup",/)')
     +    (p,p=1,maxpar)
        do 120 nen=ebegin(type),eend(type)
          nonpski=xspreeq(type,nen)-xspreeqps(type,nen)-
     +      xspreeqki(type,nen)-xspreeqbu(type,nen)
          write(*,'(1x,f7.3,1p,11e10.3)') egrid(nen),xspreeq(type,nen),
     +      (xsstep(type,p,nen),p=1,maxpar),nonpski,
     +      xspreeqps(type,nen),xspreeqki(type,nen),xspreeqbu(type,nen)
  120   continue
        nonpski=xspreeqtot(type)-xspreeqtotps(type)-xspreeqtotki(type)-
     +    xspreeqtotbu(type)
        write(*,'(/8x,1p,11e10.3)') xspreeqtot(type),
     +    (xssteptot(type,p),p=1,maxpar),nonpski,xspreeqtotps(type),
     +    xspreeqtotki(type),xspreeqtotbu(type)
        write(*,'(/" Integrated:",f12.5/)') xspreeqtot(type)
  110 continue
      write(*,'(" Total pre-equilibrium cross section:",f12.5)')
     +  xspreeqsum
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
