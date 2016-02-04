      function intri(x,y,x1,y1,x2,y2,x3,y3)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 9, 2004
c | Task  : Test if (x,y) belongs to the triangle defined by the points
c |         (x1,y1),(x2,y2),(x3,y3)
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c x,y   : coordinates of the point to test
c x1,y1 : coordinates of the 1st summit of the triangle
c x2,y2 : coordinates of the 2nd summit of the triangle
c x3,y3 : coordinates of the 3rd summit of the triangle
c xg,yg : center of gravity of the triangle (x1,y1),(x2,y2),(x3,y3)
c sign1 : help variable in talys.cmb
c signg1: help variable
c sign2 : help variable in talys.cmb
c signg2: help variable
c sign3 : help variable in talys.cmb
c signg3: help variable
c
      logical intri
      real    x,y,x1,y1,x2,y2,x3,y3,xg,yg,signg1,signg2,signg3,sideline,
     +        sign1,sign2,sign3
c
c ********************************* Method *****************************
c
c We test if the center of gravity (xg,yg) of the three summits and the
c point to test are on the same side of each of the three segments
c defining the triangle
c
      intri=.false.
c
c Center of gravity
c
      xg=(x1+x2+x3)/3.0
      yg=(y1+y2+y3)/3.0
c
c Sideline for each point
c
c first segment (x1,y1) --> (x2,y2)
c
      sign1=sideline(x,y,x1,y1,x2,y2)
      signg1=sideline(xg,yg,x1,y1,x2,y2)
c
c second segment (x2,y2) --> (x3,y3)
c
      sign2=sideline(x,y,x2,y2,x3,y3)
      signg2=sideline(xg,yg,x2,y2,x3,y3)
c
c third segment (x3,y3) --> (x1,y1)
c
      sign3=sideline(x,y,x3,y3,x1,y1)
      signg3=sideline(xg,yg,x3,y3,x1,y1)
c
c Final test
c
      if ((sign1.eq.signg1).or.(sign1.eq.0.0)) then
        if ((sign2.eq.signg2).or.(sign2.eq.0.0)) then
          if ((sign3.eq.signg3).or.(sign3.eq.0.0)) then
            intri=.true.
          endif
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
