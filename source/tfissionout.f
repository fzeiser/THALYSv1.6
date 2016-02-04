      subroutine tfissionout(Zcomp,Ncomp,nex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 30, 2013
c | Task  : Output of fission transmission coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,nex,Z,N,A,odd,J,J2
c
c *************** Output of fission transmission coefficients **********
c
c Zcomp: charge number index for compound nucleus
c Ncomp: neutron number index for compound nucleus
c nex  : excitation energy bin
c ZZ,Z : charge number of residual nucleus
c NN,N : neutron number of residual nucleus
c AA,A : mass number of residual nucleus
c nuc  : symbol of nucleus
c Exinc: excitation energy of entrance bin
c odd  : odd (1) or even (0) nucleus
c maxJ : maximal J-value
c J2   : 2 * J
c tfis : fission transmission coefficients
c
      Z=ZZ(Zcomp,Ncomp,0)
      N=NN(Zcomp,Ncomp,0)
      A=AA(Zcomp,Ncomp,0)
      write(*,'(/" Fission transmission coefficients for Z=",i3,
     +  " N=",i3," (",i3,a2,") and an excitation energy of ",f7.3,
     +  " MeV"/)') Z,N,A,nuc(Z),Exinc
      write(*,'("   J      T(J,-)      T(J,+)"/)')
      odd=mod(A,2)
      do 10 J=0,maxJ(Zcomp,Ncomp,nex)
        J2=2*J+odd
        write(*,'(1x,f4.1,2x,1p,2e12.5)') 0.5*J2,tfis(J,-1),tfis(J,1)
   10 continue
      write(*,*)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
