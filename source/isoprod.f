      subroutine isoprod
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 17, 2012
c | Task  : Calculate isotope production
c +---------------------------------------------------------------------
c
c ********** Calculation and output of production yields ***************
c
c prodinitial: subroutine for Initialization of isotope production info
c prodres    : subroutine for residual production cross sections
c prodrates  : subroutine to calculate reaction rates
c prodyield  : subroutine to calculate production yields
c prodout    : subroutine for output of isotope production
c
      call prodinitial
      call prodres
      call prodrates
      call prodyield
      call prodout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
