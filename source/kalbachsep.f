      subroutine kalbachsep
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 10, 2004
c | Task  : Separation energies for Kalbach systematics
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Z,A
      real    Ac,Ab,Zc,Zb,Nc,Nb,Ib
c
c ************** Spherical Myers-Swiatecki parameters ******************
c
c ZZ,Z    : charge number of residual nucleus
c AA,A    : mass number of residual nucleus
c Ainit   : mass number of initial compound nucleus
c Zinit   : charge number of initial compound nucleus
c Ac,.... : help variables
c Ib      : energy required for breakup
c parsym  : symbol of particle
c Smyers  : Myers-Swiatecki separation energy
c onethird: 1/3
c twothird: 2/3
c
c We use the 'old' Myers-Swiatecki parameters as used by Kalbach
c in Phys. Rev. C37, 2350, (1988).
c
      do 10 type=1,6
        Z=ZZ(0,0,type)
        A=AA(0,0,type)
        Ac=real(Ainit)
        Ab=real(A)
        Zc=real(Zinit)
        Zb=real(Z)
        Nc=Ac-Zc
        Nb=Ab-Zb
        Ib=0.
        if (parsym(type).eq.'d') Ib=2.225
        if (parsym(type).eq.'t') Ib=8.482
        if (parsym(type).eq.'h') Ib=7.718
        if (parsym(type).eq.'a') Ib=28.296
        Smyers(type)=15.68*(Ac-Ab)-28.07*((Nc-Zc)**2/Ac-(Nb-Zb)**2/Ab)-
     +    18.56*(Ac**twothird-Ab**twothird)+
     +    33.22*((Nc-Zc)**2/Ac**(2.*twothird)-(Nb-Zb)**2/Ab**
     +    (2.*twothird))-0.717*(Zc**2/(Ac**onethird)-
     +    Zb**2/(Ab**onethird))+1.211*(Zc**2/Ac-Zb**2/Ab)-Ib
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
