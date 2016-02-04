      function func1(a,b,c,t1,t2,dab)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : August 22, 2004
c | Task  : Calculation of the x,x1,x2 terms for GOE
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      real             func1,a,b,c,dab
      double precision t1,t2,f1,f2

c ****************** Calculation of the x,x1,x2 terms ******************
c
      f1=dab*(1.-t1)*((a/(1.+t1*a)+b/(1.+t1*b)+2*c/(1.-t1*c))**2)
      f2=(1.+dab)*(a*(1.+a)/((1.+t1*a)*(1.+t2*a))+
     +  b*(1.+b)/((1.+t1*b)*(1.+t2*b))+2*c*(1.-c)/((1.-t1*c)*(1.-t2*c)))
      func1=real(f1+f2)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
