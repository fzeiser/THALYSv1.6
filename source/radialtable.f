      subroutine radialtable(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Eric Bauge
c | Date  : May 11, 2011
c | Task  : Tabulated radial matter densities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*4  radchar
      character*90 radfile
      integer      Zix,Nix,Z,A,ia,i,nradjlm,j
      real         h,dr,rp0,rp1,rn0,rn1,ap,an,expo
c
c ****************** Read radial matter densities **********************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c AA,A       : mass number of residual nucleus
c radialfile : radial matter density file
c radchar    : help variable
c radialmodel: model for radial matter densities (JLM OMP only)
c radfile    : radial matter density file
c path       : directory containing structure files to be read
c lenpath    : length of pathname
c ia         : mass number from radial matter density table
c nradjlm    : number of radial points
c h          : radial step size
c numjlm     : maximum number of radial points
c radjlm     : radial points for JLM potential
c rhojlmp    : density for protons
c rhojlmn    : density for neutrons
c jlmexist   : flag for existence of tabulated radial matter density
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      if (radialfile(Zix)(1:1).ne.' ') then
        radfile=radialfile(Zix)
      else
        radchar='z   '
        write(radchar(2:4),'(i3.3)') Z
        if (radialmodel.eq.1) then
          radfile=path(1:lenpath)//'optical/jlm/goriely/'//radchar
        else
          radfile=path(1:lenpath)//'optical/jlm/hilaire/'//radchar
        endif
        inquire (file=radfile,exist=lexist)
        if (.not.lexist) return
      endif
      open (unit=2,status='old',file=radfile)
   10 read(2,'(4x,i4,i5,f7.3)',end=60) ia,nradjlm,h
      if (A.ne.ia) then
        do 20 i=1,nradjlm
          read(2,*)
   20   continue
        goto 10
      endif
c
c Hilaire file starts at 0 fm, which we skip.
c
      if (radialmodel.eq.2) then
        read(2,*)
        nradjlm=nradjlm-1
      endif
      do 30 i=1,min(nradjlm,numjlm)
        read(2,*) radjlm(Zix,Nix,i),
     +    (rhojlmp(Zix,Nix,i,j),j=1,5),(rhojlmn(Zix,Nix,i,j),j=1,5)
   30 continue
c
c Extrapolate for untabulated values
c
      if (nradjlm.lt.numjlm) then
        do 40 j=1,5
          rp0=rhojlmp(Zix,Nix,nradjlm,j)
          rp1=rhojlmp(Zix,Nix,nradjlm-1,j)
          rn0=rhojlmn(Zix,Nix,nradjlm,j)
          rn1=rhojlmn(Zix,Nix,nradjlm-1,j)
          ap=0.
          if (rp0.gt.0..and.rp1.gt.0.) ap=-log(rp0/rp1)/h
          an=0.
          if (rn0.gt.0..and.rn1.gt.0.) an=-log(rn0/rn1)/h
          do 50 i=nradjlm+1,numjlm
             dr=h*(i-nradjlm)
             if (j.eq.1) radjlm(Zix,Nix,i)=radjlm(Zix,Nix,nradjlm)+dr
             expo=ap*dr
             if (abs(expo).le.80.) rhojlmp(Zix,Nix,i,j)=rp0*exp(-expo)
             expo=an*dr
             if (abs(expo).le.80.) rhojlmn(Zix,Nix,i,j)=rn0*exp(-expo)
   50     continue
   40   continue
      endif
c
c Set JLM flags
c
c flagjlm : flag for using semi-microscopic JLM OMP
c alphaomp: alpha optical model (1=normal, 2= McFadden-Satchler,
c           3-5= folding potential)

      if (flagjlm) then
        jlmexist(Zix,Nix,1)=.true.
        jlmexist(Zix,Nix,2)=.true.
      endif
      if (alphaomp.ge.3) jlmexist(Zix,Nix,6)=.true.
   60 close (unit=2)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
