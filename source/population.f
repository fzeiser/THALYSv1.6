      subroutine population
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 21, 2012
c | Task  : Processing of pre-equilibrium spectra into population bins
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Zix,Nix,NL,nexout,nen1,na1,nb1,nen2,na2,nb2,parity,J,
     +        pc,p,h,pcpi,ppi,hpi,pcnu,pnu,hnu,nex
      real    SS,Eout,Elow,Ea1,Eb1,Ehigh,Ea2,Eb2,xsa,xsb,xslow,xshigh,
     +        xs,xscheck(0:numpar),norm
c
c ****** Process pre-equilibrium spectra into population arrays ********
c
c flagomponly: flag to execute ONLY an optical model calculation
c parskip    : logical to skip outgoing particle
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c SS,S       : separation energy per particle
c Nlast,NL   : last discrete level
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c flagmulpre : flag for multiple pre-equilibrium calculation
c mulpreZN   : logical for multiple pre-equilibrium per nucleus
c maxex      : maximum excitation energy bin for compound nucleus
c Eout       : outgoing energy
c Etotal     : total energy of compound system (target + projectile)
c Ex         : excitation energy
c egrid      : outgoing energy grid
c Elow,Ehigh : help variables
c deltaEx    : excitation energy bin for population arrays
c locate     : subroutine to find value in ordered table
c na1,nb1    : help variables
c Ea1,Eb1    : help variables
c xsa,xsb,xs : help variables
c xspreeq    : preequilibrium cross section per particle type and
c              outgoing energy
c
c The pre-equilibrium cross sections have been calculated on the
c emission energy grid. They are interpolated on the excitation
c energy grids of the level populations (both for the total and the
c spin/parity-dependent cases) to enable futher decay of the residual
c nuclides.
c
      if (flagomponly) return
      do 10 type=0,6
        if (parskip(type)) goto 10
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        SS=S(0,0,type)
        NL=Nlast(Zix,Nix,0)
        if (ebegin(type).ge.eend(type)) goto 10
c
c We include multiple pre-equilibrium emission only for neutrons and
c protons.
c
        if (flagmulpre.and.(type.eq.1.or.type.eq.2))
     +    mulpreZN(Zix,Nix)=.true.
c
c Loop over excitation energies. Determine the emission energy that
c corresponds with the excitation energy.
c
        do 20 nexout=NL+1,maxex(Zix,Nix)
          Eout=Etotal-SS-Ex(Zix,Nix,nexout)
          if (Eout.lt.egrid(ebegin(type))) goto 20
          Elow=Eout-0.5*deltaEx(Zix,Nix,nexout)
          call locate(egrid,ebegin(type),eend(type),Elow,nen1)
          na1=nen1
          nb1=nen1+1
          Ea1=egrid(na1)
          Eb1=min(egrid(nb1),Etotal-SS)
          Ehigh=Eout+0.5*deltaEx(Zix,Nix,nexout)
          call locate(egrid,ebegin(type),eend(type),Ehigh,nen2)
          na2=nen2
          nb2=nen2+1
          Ea2=egrid(na2)
          Eb2=min(egrid(nb2),Etotal-SS)
          xsa=xspreeq(type,na1)
          xsb=xspreeq(type,nb1)
c
c Add contribution from giant resonances
c
c flaggiant: flag for collective contribution from giant resonances
c xsgr     : smoothed giant resonance cross section
c pol1     : subroutine for polynomial interpolation of first order
c xslow,...: help variable
c
          if (flaggiant) then
            xsa=xsa+xsgr(type,na1)
            xsb=xsb+xsgr(type,nb1)
          endif
          call pol1(Ea1,Eb1,xsa,xsb,Elow,xslow)
          xsa=xspreeq(type,na2)
          xsb=xspreeq(type,nb2)
          if (flaggiant) then
            xsa=xsa+xsgr(type,na2)
            xsb=xsb+xsgr(type,nb2)
          endif
c
c Determine interpolated value.
c
c preeqpopex: pre-equilibrium population cross section summed over
c             spin and parity
c
          call pol1(Ea2,Eb2,xsa,xsb,Ehigh,xshigh)
          xs=0.5*(xslow+xshigh)*deltaEx(Zix,Nix,nexout)
          if (xs.lt.1.e-30) xs=0.
          preeqpopex(Zix,Nix,nexout)=xs
c
c If the pre-equilibrium spin distribution is chosen, the spectrum
c is interpolated on the spin/parity dependent population.
c
c pespinmodel: model for pre-equilibrium spin distribution or compound
c              spin distribution for pre-equilibrium cross section
c parity     : parity
c J          : total angular momentum
c maxJph     : maximal spin for particle-hole states
c xspreeqJP  : preequilibrium cross section per particle type,
c              outgoing energy, spin and parity
c xsgrstate  : smoothed giant resonance cross section per state
c preeqpop   : pre-equilibrium population cross section
c
          if (pespinmodel.eq.3) then
            do 30 parity=-1,1,2
              do 30 J=0,maxJph
                xsa=xspreeqJP(type,na1,J,parity)
                xsb=xspreeqJP(type,nb1,J,parity)
                if (flaggiant.and.J.le.3) then
                  xsa=xsa+xsgrstate(type,J,1,na1)+
     +              xsgrstate(type,J,2,na1)
                  xsb=xsb+xsgrstate(type,J,1,nb1)+
     +              xsgrstate(type,J,2,nb1)
                endif
                call pol1(Ea1,Eb1,xsa,xsb,Elow,xslow)
                xsa=xspreeqJP(type,na2,J,parity)
                xsb=xspreeqJP(type,nb2,J,parity)
                if (flaggiant.and.J.le.3) then
                  xsa=xsa+xsgrstate(type,J,1,na2)+
     +              xsgrstate(type,J,2,na2)
                  xsb=xsb+xsgrstate(type,J,1,nb2)+
     +              xsgrstate(type,J,2,nb2)
                endif
                call pol1(Ea2,Eb2,xsa,xsb,Ehigh,xshigh)
                xs=0.5*(xslow+xshigh)*deltaEx(Zix,Nix,nexout)
                if (xs.lt.1.e-30) xs=0.
                preeqpop(Zix,Nix,nexout,J,parity)=xs
   30       continue
          endif
c
c A similar interpolation is done for multiple pre-equilibrium emission.
c
          if (mulpreZN(Zix,Nix)) then
c
c 1. One-component model
c
c flag2comp: flag for two-component pre-equilibrium model
c pc       : composite particle number
c p0       : initial particle number
c maxpar   : maximal particle number
c p        : particle number
c parA     : mass number of particle
c h        : hole number
c xsstep   : preequilibrium cross section per particle type, stage
c            and outgoing energy
c xspopph  : population cross section per particle-hole configuration
c
            if (.not.flag2comp) then
              do 40 pc=p0,maxpar
                p=pc-parA(type)
                h=pc-p0
                if (p.lt.0.or.h.lt.0) goto 40
                xsa=xsstep(type,pc,na1)
                xsb=xsstep(type,pc,nb1)
                call pol1(Ea1,Eb1,xsa,xsb,Elow,xslow)
                xsa=xsstep(type,pc,na2)
                xsb=xsstep(type,pc,nb2)
                call pol1(Ea2,Eb2,xsa,xsb,Ehigh,xshigh)
                xs=0.5*(xslow+xshigh)*deltaEx(Zix,Nix,nexout)
                if (xs.lt.1.e-30) xs=0.
                xspopph(Zix,Nix,nexout,p,h)=xs
   40         continue
            else
c
c 2. Two-component model
c
c pcpi    : composite proton particle number
c ppi0    : initial proton number
c ppi     : proton particle number
c parZ    : charge number of particle
c hpi     : proton hole number
c pcnu    : composite neutron particle number
c pnu0    : initial neutron number
c pnu     : neutron particle number
c parN    : neutron number of particle
c hnu     : neutron hole number
c xsstep2 : two-component preequilibrium cross section
c xspopph2: population cross section per two-component particle-hole
c            configuration
c
              do 50 pcpi=ppi0,maxpar
                ppi=pcpi-parZ(type)
                hpi=pcpi-ppi0
                if (ppi.lt.0.or.hpi.lt.0) goto 50
                do 60 pcnu=pnu0,maxpar
                  pnu=pcnu-parN(type)
                  hnu=pcnu-pnu0
                  if (pnu.lt.0.or.hnu.lt.0) goto 60
                  xsa=xsstep2(type,pcpi,pcnu,na1)
                  xsb=xsstep2(type,pcpi,pcnu,nb1)
                  call pol1(Ea1,Eb1,xsa,xsb,Elow,xslow)
                  xsa=xsstep2(type,pcpi,pcnu,na2)
                  xsb=xsstep2(type,pcpi,pcnu,nb2)
                  call pol1(Ea2,Eb2,xsa,xsb,Ehigh,xshigh)
                  xs=0.5*(xslow+xshigh)*deltaEx(Zix,Nix,nexout)
                  if (xs.lt.1.e-30) xs=0.
                  xspopph2(Zix,Nix,nexout,ppi,hpi,pnu,hnu)=xs
   60           continue
   50         continue
            endif
          endif
   20   continue
   10 continue
c
c ***************** Correct for interpolation errors *******************
c
c xscheck   : help variable to check total population
c xspreeqtot: preequilibrium cross section per particle type
c norm      : normalization factor
c xsgrtot   : total smoothed giant resonance cross section
c
c Normalization of the pre-equilibrium population cross section. Due
c to interpolation, the part of the continuum population that comes from
c pre-equilibrium is not exactly equal to the total pre-equilibrium
c cross section. The normalization is done over the whole excitation
c energy range.
c
      do 110 type=0,6
        if (parskip(type)) goto 110
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        NL=Nlast(Zix,Nix,0)
        xscheck(type)=0.
        do 120 nex=NL+1,maxex(Zix,Nix)
          xscheck(type)=xscheck(type)+preeqpopex(Zix,Nix,nex)
  120   continue
        norm=1.
        if (xscheck(type).ne.0.)
     +    norm=(xspreeqtot(type)+xsgrtot(type))/xscheck(type)
        do 130 nex=NL+1,maxex(Zix,Nix)
          preeqpopex(Zix,Nix,nex)=preeqpopex(Zix,Nix,nex)*norm
          if (pespinmodel.eq.3) then
            do 140 parity=-1,1,2
              do 140 J=0,maxJph
                preeqpop(Zix,Nix,nex,J,parity)=
     +            preeqpop(Zix,Nix,nex,J,parity)*norm
  140       continue
          endif
  130   continue
  110 continue
c
c ************************* Check of population ************************
c
c flagcheck: flag for output of numerical checks
c parname  : name of particle
c
c The difference that was already corrected above is illustrated by the
c following table.
c
      if (flagcheck) then
        write(*,'(/" ########## POPULATION CHECK ##########"/)')
        write(*,'(" Particle Pre-equilibrium Population"/)')
        do 210 type=0,6
          if (parskip(type)) goto 210
          xs=xspreeqtot(type)
          if (flaggiant) xs=xs+xsgrtot(type)
          write(*,'(1x,a8,2(f12.5))') parname(type),xs,xscheck(type)
  210   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
