      subroutine eciscompound
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 15, 2012
c | Task  : Create ECIS input file for compound cross section
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nex,Zix,Nix,kopt,i
      real    eopt
c
c *********************** Write standard input *************************
c
c title        : title of ECIS input file
c ecis1,ecis2  : 100 input flags ('T' or 'F') for ECIS
c ncoll        : number of nuclear states
c njmax        : maximal number of j-values in ECIS
c iterm        : number of iterations
c npp          : number of optical potentials
c rmatch       : matching radius
c nsp1         : number of uncoupled states and continua
c nsp2         : number of uncoupled states with angular distribution
c ncont        : number of continua
c Einc         : incident energy in MeV
c targetspin   : spin of target
c tarparity    : parity of target
c spin         : spin of incident particle
c projmass     : mass of projectile
c resmass      : mass of target nucleus
c prodZ        : product of charges of projectile and target nucleus
c jcomp        : spin of level
c pcomp        : parity of level
c elevelcomp   : energy of level
c spincomp     : spin of incident particle
c ejeccomp     : mass of projectile
c masscomp     : mass of nucleus with indices (Z,N)
c prodZcomp    : product of charges
c Zix          : charge number index for residual nucleus
c parZ         : charge number of particle
c typecomp     : particle type
c Nix          : neutron number index for residual nucleus
c parN         : neutron number of particle
c eopt         : incident energy
c specmass     : specific mass for target nucleus and all particles
c kopt         : optical model index for fast particle
c optical      : subroutine for determination of optical potential
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c disp         : flag for dispersive optical model
c angbeg       : first angle
c anginc       : angle increment
c angend       : last angle
c bz1          : elastic enhancement factor
c parinclude   : logical to include outgoing particle
c tgo          : slow s-wave neutron gamma width/spacing
c S            : separation energy per particle
c egr          : energy of GR
c ggr          : width of GR
c parskip      : logical to skip outgoing particle
c Zinit        : charge number of initial compound nucleus
c aldcomp      : level density parameter with indices (Z,N)
c Umcomp       : matching point for U (excitation energy - pairing
c                energy)
c tempcomp     : nuclear temperature
c E0comp       : constant of temperature formula
c Excomp       : Matching Ex
c
      write(1,'(a72)') title
      write(1,'(a50)') ecis1
      write(1,'(a50)') ecis2
      write(1,'(4i5)') ncoll,njmax,iterm,npp
      write(1,'(10x,f10.5,10x,3("    1.e-10"))') rmatch
      write(1,'()')
      write(1,'(2i5,10x,i5)') nsp1,nsp2,ncont
      if (Einc.ge.0.01) then
        write(1,'(f5.2,2i2,a1,5f10.5)') targetspin,0,1,
     +    tarparity,Einc,spin,projmass,resmass,prodZ
      else
        write(1,'(f5.2,2i2,a1,1p,e10.3,0p,4f10.5)') targetspin,0,1,
     +    tarparity,Einc,spin,projmass,resmass,prodZ
      endif
      do 10 nex=1,nsp1
        write(1,'(f5.2,2i2,a1,5f10.5)') jcomp(nex),0,nex+1,pcomp(nex),
     +  elevelcomp(nex),spincomp(nex),ejeccomp(nex),masscomp(nex),
     +  prodZcomp(nex)
   10 continue
      do 20 nex=0,nsp1
        Zix=parZ(typecomp(nex))
        Nix=parN(typecomp(nex))
        eopt=Einc-real(elevelcomp(nex)/specmass(Zix,Nix,typecomp(nex)))
        eopt=max(eopt,0.001)
        kopt=typecomp(nex)
        call optical(Zix,Nix,kopt,eopt)
        write(1,'(3f10.5)') v,rv,av
        write(1,'(3f10.5)') w,rw,aw
        write(1,'(3f10.5)') vd,rvd,avd
        write(1,'(3f10.5)') wd,rwd,awd
        write(1,'(3f10.5)') vso,rvso,avso
        write(1,'(3f10.5)') wso,rwso,awso
        write(1,'(3f10.5)') rc,0.,0.
        write(1,'(3f10.5)') 0.,0.,0.
   20 continue
      write(1,'(3f10.5)') angbeg,anginc,angend
      write(1,'(f10.5)') bz1
      if (parinclude(0)) write(1,'(1p,e10.3,0p,4f10.5)')
     +  tgo,S(0,0,1),0.,egr(0,0,1,1,1),ggr(0,0,1,1,1)
      do 30 nex=0,ncont
        if (parskip(0).and.nex.eq.0) goto 30
        write(1,'(1p,7e10.3)') real(Zinit),aldcomp(nex),
     +    Umcomp(nex),tempcomp(nex),0.,E0comp(nex),Excomp(nex)
   30 continue
      write(1,'(3i5)') 1,1,0
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      if (disp(Zix,Nix,k0)) then
        do 40 i=1,npp
          write(1,'(10x,2i5)') 2,2
          write(1,'(10x,f10.5,40x,f10.5)') ef(Zix,Nix,k0),w2(Zix,Nix,k0)
          write(1,'(20x,2f10.5)') d3(Zix,Nix,k0),d2(Zix,Nix,k0)
   40   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
