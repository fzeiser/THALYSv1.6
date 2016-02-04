      subroutine ecisdwbamac(itype,type)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 17, 2009
c | Task  : Create ECIS input file for macroscopic DWBA calculation for
c |         MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 ch
      integer     itype,type,Zix,Nix,Z,J,i,Zixi,Nixi,kopt
      real        eopt
c
c **************************** Write input *****************************
c
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c disp       : flag for dispersive optical model
c title      : title of ECIS input file
c ecis1,ecis2: 100 input flags ('T' or 'F') for ECIS
c ncoll      : number of nuclear states
c njmax      : maximal number of j-values in ECIS
c iterm      : number of iterations
c npp        : number of optical potentials
c rmatch     : matching radius
c Emsdin     : incident MSD energy
c parspin    : spin of particle
c parmass    : mass of particle in a.m.u.
c nucmass    : mass of nucleus
c parZ       : charge number of particle
c parN       : neutron number of particle
c Zinit      : charge number of initial compound nucleus
c maxJmsd    : maximal spin for MSD calculation
c ch         : help variable
c Exmsd      : excitation energy for MSD energy grid
c betamsd    : deformation parameter
c
      Zix=Zindex(0,0,type)
      Nix=Nindex(0,0,type)
      Z=ZZ(0,0,type)
      if (disp(Zix,Nix,itype)) then
        ecis1(10:10)='T'
      else
        ecis1(10:10)='F'
      endif
      write(9,'(a72)') title
      write(9,'(a50)') ecis1
      write(9,'(a50)') ecis2
      write(9,'(4i5)') ncoll,njmax,iterm,npp
      write(9,'(10x,f10.5,10x,3("    1.e-10"))') rmatch
      write(9,'(f5.2,2i2,a1,5f10.5)') 0.,1,2,'+',Emsdin,
     +  parspin(itype),parmass(itype),nucmass(parZ(itype),parN(itype)),
     +  real(Zinit-parZ(itype))*parZ(itype)
      do 10 J=0,maxJmsd
        if (mod(J,2).eq.0) then
          ch='+'
        else
          ch='-'
        endif
        write(9,'(f5.2,2i2,a1,5f10.5)') real(J),0,3,ch,Exmsd,
     +    parspin(type),parmass(type),nucmass(Zix,Nix),
     +    real(Z*parZ(type))
        write(9,'(2i5)') 1,J+1
   10 continue
      do 20 J=0,maxJmsd
        write(9,'(i5,5x,f10.5)') J,betamsd
   20 continue
      do 30 i=1,3
c
c Transition potential and potential for incident channel
c
c Zixi,Nixi: help variables
c kopt     : index for fast particle
c eopt     : incident energy (help variable)
c Emsdout  : outgoing MSD energy
c
        if (i.le.2) then
          Zixi=parZ(itype)
          Nixi=parN(itype)
          kopt=itype
          if (i.eq.1) then
            eopt=0.5*(Emsdin+Emsdout)
          else
            eopt=Emsdin
          endif
        endif
c
c Potential for outgoing channel
c
c optical      : subroutine for determination of optical potential
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c angbeg       : first angle
c anginc       : angle increment
c angend       : last angle
c
        if (i.eq.3) then
          Zixi=parZ(type)
          Nixi=parN(type)
          kopt=type
          eopt=Emsdout
        endif
        call optical(Zixi,Nixi,kopt,eopt)
        write(9,'(3f10.5)') v,rv,av
        write(9,'(3f10.5)') w,rw,aw
        write(9,'(3f10.5)') vd,rvd,avd
        write(9,'(3f10.5)') wd,rwd,awd
        write(9,'(3f10.5)') vso,rvso,avso
        write(9,'(3f10.5)') wso,rwso,awso
        write(9,'(3f10.5)') rc,0.,0.
        write(9,'(3f10.5)') 0.,0.,0.
   30 continue
      write(9,'(3f10.5)') angbeg,anginc,angend
      if (disp(Zixi,Nixi,itype)) then
        do 40 i=1,3
          write(9,'(10x,2i5)') 2,2
          write(9,'(10x,f10.5,40x,f10.5)') ef(Zixi,Nixi,itype),
     +      w2(Zixi,Nixi,itype)
          write(9,'(20x,2f10.5)') d3(Zixi,Nixi,itype),
     +      d2(Zixi,Nixi,itype)
   40   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
