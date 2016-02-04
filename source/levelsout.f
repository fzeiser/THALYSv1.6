      subroutine levelsout(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 3, 2009
c | Task  : Output of discrete levels
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,Z,N,A,type,i,k
c
c *************************** Discrete levels **************************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c NN,N       : neutron number of residual nucleus
c AA,A       : mass number of residual nucleus
c nuc        : symbol of nucleus
c nucmass    : mass of nucleus
c parskip    : logical to skip outgoing particle
c parname    : name of particle
c S          : separation energy per particle
c nlev       : number of levels for nucleus
c edis       : energy of level
c jdis       : spin of level
c cparity    : parity of level (character)
c parlev     : parity of level
c tau        : lifetime of state in seconds
c jassign    : flag for assignment of spin
c passign    : flag for assignment of parity
c ENSDF      : string from original ENSDF discrete level file
c nbranch    : number of branching levels
c branchratio: gamma-ray branching ratio to level
c bassign    : flag for assignment of branching ratio
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      write(*,'(/" NUCLEAR STRUCTURE INFORMATION FOR Z=",i3," N=",i3,
     +  " (",i3,a2,") ")') Z,N,A,nuc(Z)
      write(*,'(/" Mass in a.m.u.     : ",f10.6)') nucmass(Zix,Nix)
      write(*,'(/" Separation energies:")')
      write(*,'(/" Particle        S         "/)')
      do 10 type=1,6
        if (parskip(type)) goto 10
        write(*,'(1x,a8,f12.5)') parname(type),S(Zix,Nix,type)
   10 continue
      write(*,'(/" Discrete levels of Z=",i3," N=",i3," (",i3,a2,") ")')
     +   Z,N,A,nuc(Z)
      write(*,'(/" Number  Energy Spin Parity  Branching  ",
     +  "Ratio (%) Lifetime(sec) Assignment        ENSDF"/)')
      do 20 i=0,nlev(Zix,Nix)
        if (tau(Zix,Nix,i).ne.0.) then
          write(*,'(1x,i3,4x,f7.4,1x,f4.1,3x,a1,24(" "),2x,1p,e10.3,
     +      7x,2a1,a18)') i,edis(Zix,Nix,i),jdis(Zix,Nix,i),
     +      cparity(parlev(Zix,Nix,i)),tau(Zix,Nix,i),
     +      jassign(Zix,Nix,i),passign(Zix,Nix,i),ENSDF(Zix,Nix,i)
        else
          write(*,'(1x,i3,4x,f7.4,1x,f4.1,3x,a1,36(" "),7x,2a1,a18)') i,
     +      edis(Zix,Nix,i),jdis(Zix,Nix,i),cparity(parlev(Zix,Nix,i)),
     +      jassign(Zix,Nix,i),passign(Zix,Nix,i),ENSDF(Zix,Nix,i)
        endif
        do 30 k=1,nbranch(Zix,Nix,i)
          write(*,'(31x,"--->",i3,2x,f8.4,18x,a1)')
     +      branchlevel(Zix,Nix,i,k),branchratio(Zix,Nix,i,k)*100.,
     +      bassign(Zix,Nix,i,k)
   30   continue
   20 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
