      function invect(x,y,x1,y1,x2,y2)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 9, 2004
c | Task  : Test if (x,y) belongs to the segment defined by the points
c |         (x1,y1),(x2,y2)
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c x,y    : coordinates of the point to test
c x1,y1  : coordinates of the 1st point of the segment
c x2,y2  : coordinates of the 2nd point of the segment
c x1x,y1x: coordinates of the vector (x-x1,y-y1)
c x12,y12: coordinates of the vector (x2-x1,y2-y1)
c pscal  : scalar product of the vector (x-x1,y-y1) with (x2-x1,y2-y1)
c norm1  : norm of vector (x-x1,y-y1)
c norm2  : norm of vector (x2-x1,y2-y1)
c
c ********************************* Method *****************************
c
c We test if the point (x,y) belongs to the segment by testing if the
c scalar product of the vector (x-x1,y-y1) with the vector (x2-x1,y2-y1)
c is positive and if the norm of the vector (x-x1,y-y1) is lower than
c the norm of the vector (x2-x1,y2-y1)
c
      logical invect
      real    x,y,x1,y1,x2,y2,x1x,y1x,x12,y12,norm1,norm2,epsilon,pscal,
     +        dnorm
c
c *************************** Initialisations **************************
c
      invect=.false.
      x1x=x-x1
      y1x=y-y1
      x12=x2-x1
      y12=y2-y1
c
c norms and scalar product
c
      norm1=x1x*x1x+y1x*y1x
      norm2=x12*x12+y12*y12
      epsilon=norm2/1.0e14
      pscal=x12*x1x+y12*y1x
c
c test if the two vectors have the same direction (if not return)
c
      if (pscal.lt.-epsilon) return
c
c test if norm of vector (x-x1,y-y1) lower than norm of segment
c
      dnorm=norm1-norm2
      if (dnorm.gt.epsilon) return
      invect=.true.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
