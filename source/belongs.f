      function belongs(x,x1,x2)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 9, 2004
c | Task  : Test if x belongs to the interval [xlow,xup] or [xup,xlow]
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c x     : value to test
c xl1   : 1st limit
c xl2   : 2nd limit
c ratio : if 0<ratio<1 then x is between xl1 and xl2
c
      logical  belongs
      real     x,x1,x2,x1x,x1x2,norm12,norm1x,epsilon,pscal,dnorm
c
c **********************************************************************
c
      belongs=.true.
      x1x=x-x1
      x1x2=x2-x1
      norm12=abs(x1x2)
      norm1x=abs(x1x)
      epsilon=norm12/1.0e14
      pscal=x1x*x1x2
      dnorm=norm1x-norm12
c
c test if vector (x-x1,0) and vector (x2-x1,0) have the same direction
c
      if (pscal.lt.-epsilon) belongs=.false.
c
c test if norm of vector (x-x1,0) lower than norm of vector (x2-x1,0)
c
      if (dnorm.gt.epsilon) belongs=.false.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
