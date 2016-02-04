      subroutine inverseread(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning, Stephane Hilaire, Eric Bauge and Pascal Romain
c | Date  : April 27, 2013
c | Task  : Read ECIS results for outgoing particles and energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist
      integer          Zcomp,Ncomp,type,nen,Zix,Nix,nJ,i,nS,k,lev,l,
     +                 ispin
      real             groundspin2,rj,factor,jres
      double precision xs,Tcoef,teps
c
c ************ Read total, reaction and elastic cross section **********
c
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c csfile    : file with inverse reaction cross sections
c parskip   : logical to skip outgoing particle
c ebegin    : first energy point of energy grid
c eendmax   : last energy point of energy grid for maximum incident
c             energy
c xs        : help variable
c xstot     : total cross section (neutrons only) for
c xsreac    : reaction cross section
c xsopt     : optical model reaction cross section
c xselas    : total elastic cross section (neutrons only)
c ecisstatus: status of ECIS file
c
      inquire (file=csfile,exist=lexist)
      if (.not.lexist) then
        write(*,'(" TALYS-error: The first calculation of a run",
     +    " should always be done with ecissave y and",
     +    " eciscalc y")')
        stop
      endif
      open (unit=3,status='unknown',file=csfile)
      do 10 type=1,6
        if (parskip(type)) goto 10
        do 20 nen=ebegin(type),eendmax(type)
          read(3,'()')
          if (type.eq.1) then
            read(3,*) xs
            xstot(type,nen)=max(real(xs),0.)
          endif
          read(3,*) xs
          xsreac(type,nen)=max(real(xs),0.)
          xsopt(type,nen)=xsreac(type,nen)
          if (type.eq.1) then
            read(3,*) xs
            xselas(type,nen)=max(real(xs),0.)
          endif
   20   continue
   10 continue
      close (unit=3,status=ecisstatus)
c
c ******************* Read transmission coefficients *******************
c
c transfile  : file with transmission coefficients
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c groundspin2: 2 * spin of ground state
c jdis       : spin of level
c nJ         : number of total J values for transmission coefficients
c rJ         : compound nucleus spin
c nS         : number of for transmission coefficients per J-value
c lev        : level number
c l          : orbital angular momentum
c jres       : j-value
c Tjl,Tcoef  : transmission coefficients as a function of
c              particle type, energy, spin and l-value
c numl       : maximal number of l-values in TALYS
c colltype   : type of collectivity (D, V or R)
c flagrot    : flag for use of rotational optical model per
c              outgoing particle, if available
c factor     : help variable
c ispin      : spin index
c parspin    : spin of particle
c
c The transmission coefficient Tjl has four indices: particle type,
c energy, spin and  l-value. For spin-1/2 particles, we use the
c indices -1 and 1 for the two spin values. For spin-1 particles,
c we use -1, 0 and 1 and for spin-0 particles we use 0 only.
c
c For rotational nuclei, the rotational transmission
c coefficients are transformed into into spherical equivalents.
c
      open (unit=7,status='unknown',file=transfile)
      do 110 type=1,6
        if (parskip(type)) goto 110
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        groundspin2=int(2.*jdis(Zix,Nix,0))
        do 120 nen=ebegin(type),eendmax(type)
          read(7,'(55x,i5)') nJ
          do 130 i=1,nJ
            read(7,'(f10.1,5x,i5)') rj,nS
            do 140 k=1,nS
              read(7,'(i3,i6,f9.1,e20.8)') lev,l,jres,Tcoef
              if (l.gt.numl) goto 140
              if (lev.eq.1) then
                if (colltype(Zix,Nix).ne.'S'.and.flagrot(type)) then
                  factor=(2.*rj+1.)/(2.*jres+1.)/(groundspin2+1.)
                else
                  factor=1.
                endif
                if (parspin(type).eq.0.5) then
                  ispin=int(2.*(jres-real(l)))
                else
                  ispin=int(jres-real(l))
                endif
                Tjl(type,nen,ispin,l)=Tjl(type,nen,ispin,l)+
     +            factor*max(real(Tcoef),0.)
              endif
  140       continue
  130     continue
  120   continue
  110 continue
c
c ************** Processing of transmission coefficients ***************
c
c Transmission coefficients averaged over spin and determination of
c maximal l-value. ECIS stops its output of transmission coefficients
c somewhat too early. For the highest l values the transmission
c coefficient for (l+spin) is not written in the output. Since these
c are small numbers we put them equal to the value for (l-spin).
c
c eend      : last energy point of energy grid
c Tl        : transmission coefficients as a function of
c             particle type, energy and l-value (averaged over spin)
c translimit: limit for transmission coefficient
c teps      : help variable
c transeps  : limit for transmission coefficient
c lmax      : maximal l-value for transmission coefficients
c
      do 210 type=1,6
        if (parskip(type)) goto 210
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
        if (type.ne.3.and.type.ne.6) then
          do 220 nen=ebegin(type),eend(type)
            do 230 l=0,numl
              if (Tjl(type,nen,-1,l).ne.0.and.Tjl(type,nen,1,l).eq.0)
     +          Tjl(type,nen,1,l)=Tjl(type,nen,-1,l)
              if (Tjl(type,nen,-1,l).eq.0.and.Tjl(type,nen,1,l).ne.0.
     +          and.l.gt.0) Tjl(type,nen,-1,l)=Tjl(type,nen,1,l)
              Tl(type,nen,l)=((l+1)*Tjl(type,nen,1,l)+
     +          l*Tjl(type,nen,-1,l))/(2*l+1)
              teps=Tl(type,nen,0)*translimit/(2*l+1)
              teps=max(teps,transeps)
              if (Tjl(type,nen,-1,l).lt.teps.and.Tjl(type,nen,1,l).
     +          lt.teps) then
                lmax(type,nen)=l-1
                goto 220
              endif
              lmax(type,nen)=l
  230       continue
  220     continue
        endif
c
c 2. Spin 1 particles: Deuterons
c
        if (type.eq.3) then
          do 240 nen=ebegin(type),eend(type)
            do 250 l=0,numl
              if (Tjl(type,nen,-1,l).ne.0.and.Tjl(type,nen,0,l).eq.0)
     +          Tjl(type,nen,0,l)=Tjl(type,nen,-1,l)
              if (Tjl(type,nen,-1,l).ne.0.and.Tjl(type,nen,1,l).eq.0)
     +          Tjl(type,nen,1,l)=Tjl(type,nen,-1,l)
              if (Tjl(type,nen,-1,l).eq.0.and.Tjl(type,nen,1,l).ne.0.
     +          and.l.gt.0)
     +          Tjl(type,nen,-1,l)=Tjl(type,nen,1,l)
              Tl(type,nen,l)=((2*l+3)*Tjl(type,nen,1,l)+
     +          (2*l+1)*Tjl(type,nen,0,l)+
     +          (2*l-1)*Tjl(type,nen,-1,l))/(3*(2*l+1))
              teps=Tl(type,nen,0)*translimit/(2*l+1)
              teps=max(teps,transeps)
              if (Tjl(type,nen,-1,l).lt.teps.and.Tjl(type,nen,0,l).
     +          lt.teps.and.Tjl(type,nen,1,l).lt.teps) then
                lmax(type,nen)=l-1
                goto 240
              endif
              lmax(type,nen)=l
  250       continue
  240     continue
        endif
c
c 3. Spin 0 particles: Alpha-particles
c
        if (type.eq.6) then
          do 260 nen=ebegin(type),eend(type)
            do 270 l=0,numl
              Tl(type,nen,l)=Tjl(type,nen,0,l)
              teps=Tl(type,nen,0)*translimit/(2*l+1)
              teps=max(teps,transeps)
              if (Tl(type,nen,l).lt.teps) then
                lmax(type,nen)=l-1
                goto 260
              endif
              lmax(type,nen)=l
  270       continue
  260     continue
        endif
  210 continue
      close (unit=7,status=ecisstatus)
      open (unit=10,status='unknown',file='ecis.invin')
      close (unit=10,status=ecisstatus)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
