      subroutine rotclass2(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Pascal Romain
c | Date  : March 27, 2012
c | Task  : Build rotational bands on class2 states
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,nbi,i,itstot,its,jts,pa1,pa2
      real    Ecut,jstart,jstep,Erk10,rj,Erot,Eband,en1,en2,rj1,rj2
c
c ******************* Rotational bands building ************************
c
c Calculation of maximum energy of class2 states
c
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c nfisbar   : number of fission barriers
c Emaxclass2: maximum energy for class2 states
c fbaradjust: adjustable factors for fission parameters
c             (default 1.)
c fbarrier  : height of fission barriers
c fecont    : start of continuum energy
c widthc2   : width of class2 states
c
      if (nfisbar(Zix,Nix).eq.2) then
        Emaxclass2(Zix,Nix,1)=fbaradjust(Zix,Nix,1)*fbarrier(Zix,Nix,1)+
     +    fecont(Zix,Nix,1)+0.5*widthc2(Zix,Nix,1)
      endif
      if (nfisbar(Zix,Nix).eq.3) then
        Emaxclass2(Zix,Nix,1)=fbaradjust(Zix,Nix,1)*fbarrier(Zix,Nix,1)+
     +    fecont(Zix,Nix,1)+0.5*widthc2(Zix,Nix,1)
        Emaxclass2(Zix,Nix,2)=fbaradjust(Zix,Nix,2)*fbarrier(Zix,Nix,2)+
     +    fecont(Zix,Nix,2)+0.5*widthc2(Zix,Nix,2)
      endif
c
c Rotational bands construction
c
c nclass2     : number of sets of class2 states
c nfisc2rot   : number of rotational class2 states per set
c Ecut        : help variable
c nfisc2hb    : number of head band class2 states per set
c jstart,jstep: help variables
c Erk10       : energy shift between k=0- and k=1-
c jfisc2hb    : spin of head band class2 states
c pfisc2hb    : parity of head band class2 states
c rj,itstot   : help variables
c Erot        : rotational energy
c minertc2    : moment of inertia for class2 states set
c Eband       : help variable
c numrot      : number of rotational states
c efisc2hb    : energy of head band class2 states
c efisc2rot   : energy of rotational class2 states
c jfisc2rot   : spin of rotational class2 states
c pfisc2rot   : parity of rotational class2 states
c
      do 10 nbi=1,nclass2(Zix,Nix)
        nfisc2rot(Zix,Nix,nbi)=0
        Ecut=Emaxclass2(Zix,Nix,nbi)
        do 20 i=1,nfisc2hb(Zix,Nix,nbi)
          jstart=jfisc2hb(Zix,Nix,nbi,i)
          jstep=1.
          Erk10=0.
          if (jfisc2hb(Zix,Nix,nbi,i).eq.0) jstep=2.
          if ((jfisc2hb(Zix,Nix,nbi,i).eq.0).and.
     +      (pfisc2hb(Zix,Nix,nbi,i).eq.-1)) then
            jstart=1.
            Erk10=1./minertc2(Zix,Nix,nbi)
          endif
          rj=jstart-jstep
   30     rj=rj+jstep
          Erot=(rj*(rj+1.)-jstart*(jstart+1.))/
     +      (2.*minertc2(Zix,Nix,nbi))
          Eband=efisc2hb(Zix,Nix,nbi,i)+Erot+Erk10
          if (Eband.gt.Ecut) goto 20
          nfisc2rot(Zix,Nix,nbi)=nfisc2rot(Zix,Nix,nbi)+1
          itstot=nfisc2rot(Zix,Nix,nbi)
          if (itstot.gt.numrot) then
            write(*,'(" TALYS-error: Due to the number of",
     +        " class 2 states, the number of rotational states",
     +        " exceeds ",i3)') numrot
            write(*,'(" numrot in talys.cmb should be increased")')
            stop
          endif
          efisc2rot(Zix,Nix,nbi,itstot)=Eband
          jfisc2rot(Zix,Nix,nbi,itstot)=rj
          pfisc2rot(Zix,Nix,nbi,itstot)=pfisc2hb(Zix,Nix,nbi,i)
          goto 30
   20   continue
   10 continue
c
c ****** Reorganize class2 states by increasing excitation energy ******
c
c en1,en2,rj1,rj2,pa1,pa2: help variables
c
      do 110 nbi=1,nclass2(Zix,Nix)
        do 120 its=1,nfisc2rot(Zix,Nix,nbi)-1
          do 130 jts=its+1,nfisc2rot(Zix,Nix,nbi)
            en1=efisc2rot(Zix,Nix,nbi,its)
            rj1=jfisc2rot(Zix,Nix,nbi,its)
            pa1=pfisc2rot(Zix,Nix,nbi,its)
            en2=efisc2rot(Zix,Nix,nbi,jts)
            rj2=jfisc2rot(Zix,Nix,nbi,jts)
            pa2=pfisc2rot(Zix,Nix,nbi,jts)
            if (en1.gt.en2) then
              efisc2rot(Zix,Nix,nbi,its)=en2
              jfisc2rot(Zix,Nix,nbi,its)=rj2
              pfisc2rot(Zix,Nix,nbi,its)=pa2
              efisc2rot(Zix,Nix,nbi,jts)=en1
              jfisc2rot(Zix,Nix,nbi,jts)=rj1
              pfisc2rot(Zix,Nix,nbi,jts)=pa1
            endif
  130     continue
  120   continue
  110 continue
      if (nfisbar(Zix,Nix).eq.2) Emaxclass2(Zix,Nix,1)=
     +  Emaxclass2(Zix,Nix,1)-0.5*widthc2(Zix,Nix,1)
      if (nfisbar(Zix,Nix).eq.3) then
        Emaxclass2(Zix,Nix,1)=Emaxclass2(Zix,Nix,1)-
     +    0.5*widthc2(Zix,Nix,1)
        Emaxclass2(Zix,Nix,2)=Emaxclass2(Zix,Nix,2)-
     +    0.5*widthc2(Zix,Nix,2)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
