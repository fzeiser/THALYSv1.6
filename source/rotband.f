      subroutine rotband(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 10, 2007
c | Task  : Build rotational bands on transition states
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,nbi,i,itstot,its,jts,pa1,pa2
      real    jstart,jstep,Erk10,rj,Erot,Eband,en1,en2,rj1,rj2
c
c ******************* Rotational bands building ************************
c
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c nfisbar     : number of fission barrier parameters
c nfistrrot   : number of rotational transition states for barrier
c nfistrhb    : number of head band transition states for barrier
c jstart,jstep: help variables
c jfistrhb    : spin of head band transition states
c pfistrhb    : parity of head band transition states
c Erk10       : energy shift between k=0- and k=1-
c rj,itstot   : help variables
c Erot        : rotational energy
c minertia    : moment of inertia of fission barrier
c Eband       : help variable
c efistrhb    : energy of head band transition states
c fecont      : start of continuum energy
c numrot      : number of rotational states
c efistrrot   : energy of rotational transition states
c jfistrrot   : spin of rotational transition states
c pfistrrot   : parity of rotational transition states
c
      do 10 nbi=1,nfisbar(Zix,Nix)
        nfistrrot(Zix,Nix,nbi)=0
        do 20 i=1,nfistrhb(Zix,Nix,nbi)
          jstart=jfistrhb(Zix,Nix,nbi,i)
          jstep=1.
          Erk10=0.
          if (jfistrhb(Zix,Nix,nbi,i).eq.0) jstep=2.
          if ((jfistrhb(Zix,Nix,nbi,i).eq.0).and.
     +      (pfistrhb(Zix,Nix,nbi,i).eq.-1)) then
            jstart=1.
            Erk10=1./minertia(Zix,Nix,nbi)
          endif
          rj=jstart-jstep
   30     rj=rj+jstep
          Erot=(rj*(rj+1.)-jstart*(jstart+1.))/
     +      (2.*minertia(Zix,Nix,nbi))
          Eband=efistrhb(Zix,Nix,nbi,i)+Erot+Erk10
          if (Eband.gt.fecont(Zix,Nix,nbi)) goto 20
          nfistrrot(Zix,Nix,nbi)=nfistrrot(Zix,Nix,nbi)+1
          itstot=nfistrrot(Zix,Nix,nbi)
          if (itstot.gt.numrot) then
            write(*,'(" TALYS-error: Due to the number of",
     +        " head band states, the number of rotational states",
     +        " exceeds ",i3)') numrot
            write(*,'(" numrot in talys.cmb should be increased")')
            stop
          endif
          efistrrot(Zix,Nix,nbi,itstot)=Eband
          jfistrrot(Zix,Nix,nbi,itstot)=rj
          pfistrrot(Zix,Nix,nbi,itstot)=pfistrhb(Zix,Nix,nbi,i)
          goto 30
   20   continue
   10 continue
c
c *** Reorganize transition states by increasing excitation energy *****
c
c en1,en2,rj1,rj2,pa1,pa2: help variables
c
      do 110 nbi=1,nfisbar(Zix,Nix)
        do 120 its=1,nfistrrot(Zix,Nix,nbi)-1
          do 130 jts=its+1,nfistrrot(Zix,Nix,nbi)
            en1=efistrrot(Zix,Nix,nbi,its)
            rj1=jfistrrot(Zix,Nix,nbi,its)
            pa1=pfistrrot(Zix,Nix,nbi,its)
            en2=efistrrot(Zix,Nix,nbi,jts)
            rj2=jfistrrot(Zix,Nix,nbi,jts)
            pa2=pfistrrot(Zix,Nix,nbi,jts)
            if (en1.gt.en2) then
              efistrrot(Zix,Nix,nbi,its)=en2
              jfistrrot(Zix,Nix,nbi,its)=rj2
              pfistrrot(Zix,Nix,nbi,its)=pa2
              efistrrot(Zix,Nix,nbi,jts)=en1
              jfistrrot(Zix,Nix,nbi,jts)=rj1
              pfistrrot(Zix,Nix,nbi,jts)=pa1
            endif
  130     continue
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
