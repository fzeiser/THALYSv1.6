      subroutine weakcoupling(Zix,Nix,type)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 5, 2004
c | Task  : Weak coupling model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,type,Z,N,A,Zcore,Ncore,Zixcore,Nixcore,i,i2,
     +        middle,j,i3,i4,i5
      real    spinbeg,spinend,factor
c
c ************************* Weak coupling model ************************
c
c Zix  : charge number index for residual nucleus
c Nix  : neutron number index for residual nucleus
c ZZ,Z : charge number of residual nucleus
c NN,N : neutron number of residual nucleus
c AA,A : mass number of residual nucleus
c Zcore: charge number of core nucleus
c core : even-even core for weakcoupling (-1 or 1)
c Ncore: neutron number of core nucleus
c Zinit: charge number of initial compound nucleus
c Ninit: neutron number of initial compound nucleus
c
c The even-even core is determined by subtracting or adding the extra
c nucleon.
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      if (mod(Z,2).eq.1) then
        Zcore=Z+core
        Ncore=N
      else
        Zcore=Z
        Ncore=N+core
      endif
      Zix=abs(Zinit-Zcore)
      Nix=abs(Ninit-Ncore)
      Z=Zcore
c
c The levels and deformation parameters for the core nucleus are
c retrieved.
c
c strucexist: flag to determine whether structure info for this nucleus
c             exists
c levels    : subroutine for discrete levels
c deformpar : subroutine for deformation parameters
c Zixcore   : charge number index for core nucleus
c Nixcore   : neutron number index for core nucleus
c parZ      : charge number of particle
c parN      : neutron number of particle
c
      if (.not.strucexist(Zix,Nix)) call levels(Zix,Nix)
      call deformpar(Zix,Nix)
      Zixcore=Zix
      Nixcore=Nix
      Zix=parZ(type)
      Nix=parN(type)
c
c The levels of the core nucleus are distributed over the odd nucleus
c by angular momentum splitting.
c
c 10: Loop over core states
c
c numlev2   : maximum number of levels
c k0        : index of incident particle
c leveltype : type of level (rotational (R) or vibrational (V))
c deform    : deformation parameter
c edis      : energy of level
c
      do 10 i=0,numlev2
        if (i.eq.0.and.type.eq.k0) goto 10
        if (leveltype(Zix,Nix,i).eq.'R'.or.leveltype(Zix,Nix,i).eq.'V')
     +    goto 10
        if (deform(Zixcore,Nixcore,i).ne.0.) then
c
c Determine corresponding energy region in odd nucleus
c
          middle=0
          do 20 i2=1,numlev2
            if (edis(Zix,Nix,i2).gt.edis(Zixcore,Nixcore,i)) then
              middle=i2+1
              goto 30
            endif
   20     continue
c
c 40: Loop over possible split levels
c
c spinbeg   : begin of possible spin values
c spinend   : end of possible spin values
c jdis      : spin of level
c parlev    : parity of level
c
   30     spinbeg=abs(jdis(Zix,Nix,0)-jdis(Zixcore,Nixcore,i))
          spinend=jdis(Zix,Nix,0)+jdis(Zixcore,Nixcore,i)
          do 40 j=int(spinbeg),int(spinend)
c
c Search for odd-spin level that corresponds with j. Start from the
c level found in loop 20, and work outwards until level is found.
c Then assign parameters to this level.
c
c nlevmax2: maximum number of levels
c deftype : deformation length (D) or parameter (B)
c jcore   : spin of level of core nucleus
c pcore   : parity of level of core nucleus
c factor  : help variable
c
            do 50 i3=0,numlev2
              do 60 i4=1,-1,-2
                i5=middle+i4*i3
                if (i3.eq.0.and.i4.eq.-1) goto 50
                if (i5.lt.1.or.i5.gt.nlevmax2(Zix,Nix)) goto 60
                if (deform(Zix,Nix,i5).ne.0.) goto 60
                if (leveltype(Zix,Nix,i5).ne.'D') goto 60
                if (int(jdis(Zix,Nix,i5)).eq.j) then
                  jcore(Zix,Nix,i5)=jdis(Zixcore,Nixcore,i)
                  pcore(Zix,Nix,i5)=parlev(Zixcore,Nixcore,i)
                  factor=sqrt((2.*(j+0.5)+1.)/
     +              ((2.*jdis(Zixcore,Nixcore,i)+1.)*
     +              (2.*jdis(Zix,Nix,0)+1.)))
                  deform(Zix,Nix,i5)=factor*deform(Zixcore,Nixcore,i)
                  if (deftype(Zix,Nix).eq.'D'.and.
     +              deftype(Zixcore,Nixcore).eq.'B') deform(Zix,Nix,i5)=
     +              deform(Zix,Nix,i5)*1.24*(A**onethird)
                  if (deftype(Zix,Nix).eq.'B'.and.
     +              deftype(Zixcore,Nixcore).eq.'D') deform(Zix,Nix,i5)=
     +              deform(Zix,Nix,i5)/(1.24*(A**onethird))
                  goto 40
                endif
   60         continue
   50       continue
   40     continue
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
