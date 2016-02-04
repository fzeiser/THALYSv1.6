      subroutine excitonout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 5, 2006
c | Task  : Output of exciton model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer p,h,n,type
      real    lambdaplus
c
c ************************ Exciton model *******************************
c
      write(*,'(/" ++++++++++ EXCITON MODEL ++++++++++")')
c
c 1. Output of matrix element
c
c Ecomp     : total energy of composite system
c M2constant: constant for matrix element in exciton model
c p         : particle number
c p0        : initial particle number
c maxpar    : maximal particle number
c h         : hole number
c n         : exciton number
c matrix    : subroutine for matrix element for exciton model
c Ainit     : mass number of initial compound nucleus
c M2        : square of matrix element
c
      write(*,'(/" 1. Matrix element for E= ",f7.3/)') Ecomp
      write(*,'(" Constant for matrix element: ",f7.3/)') M2constant
      write(*,'("  p h       M2"/)')
      do 10 p=p0,maxpar
        h=p-p0
        n=p+h
        call matrix(Ainit,n)
        write(*,'(1x,2i2,2x,1p,e12.5)') p,h,M2
   10 continue
c
c 2. Output of Q-factors
c
c parname: name of particle
c Qfactor: Q-factor for neutron/proton distinction
c
      write(*,'(/" 2. Q-factors"/)')
      write(*,'("  p h  ",6(a8,1x),/)') (parname(type),type=1,6)
      do 20 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,7f9.5)') p,h,(Qfactor(type,p),type=1,6)
   20 continue
c
c 3. Output of emission rates or escape widths
c
c wemispart: emission rate per particle and exciton number
c wemistot : total emission rate per exciton number
c hbar     : Planck's constant / 2.pi in MeV.s
c
      write(*,'(/" 3. Emission rates or escape widths"/)')
      write(*,'(" A. Emission rates ( /sec)"/)')
      write(*,'("  p h",3x,7(a8,4x),"Total"/)') (parname(type),type=0,6)
      do 30 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,1p,8e12.5)') p,h,
     +    (wemispart(type,p,h),type=0,6),wemistot(p,h)
   30 continue
      write(*,'(/" B. Escape widths (MeV)"/)')
      write(*,'("  p h",2x,7(a8,4x),"Total"/)') (parname(type),type=0,6)
      do 40 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,1p,8e12.5)') p,h,
     +    (wemispart(type,p,h)*hbar,type=0,6),wemistot(p,h)*hbar
   40 continue
c
c 4. Output of transition rates or damping widths and total widths
c
c lambdaplus: transition rate for n --> n+2
c
      write(*,'(/" 4. Internal transition rates or damping widths,",
     +  " total widths"/)')
      write(*,'(" A. Internal transition rates ( /sec)"/)')
      write(*,'("  p h    lambdaplus"/)')
      do 50 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,1p,e15.5)') p,h,lambdaplus(0,0,p,h)
   50 continue
      write(*,'(/" B. Damping widths (MeV)"/)')
      write(*,'("  p h     gammaplus"/)')
      do 60 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,1p,e15.5)') p,h,lambdaplus(0,0,p,h)*hbar
   60 continue
      write(*,'(/" C. Total widths (MeV)"/)')
      write(*,'("  p h     gammatot"/)')
      do 70 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,1p,e15.5)') p,h,
     +    (lambdaplus(0,0,p,h)+wemistot(p,h))*hbar
   70 continue
c
c 5. Output of depletion factors
c
c depletion: depletion factor at each stage
c
      write(*,'(/" 5. Depletion factors"/)')
      write(*,'("  p h  depletion"/)')
      do 80 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,f10.5)') p,h,depletion(p,h)
   80 continue
c
c 6. Output of lifetimes of exciton states
c
c tauexc: mean lifetime
c
      write(*,'(/" 6. Lifetimes")')
      write(*,'(/"  p h   mean lifetime"/)')
      do 90 p=p0,maxpar
        h=p-p0
        write(*,'(1x,2i2,1p,e15.5)') p,h,tauexc(p,h)
   90 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
