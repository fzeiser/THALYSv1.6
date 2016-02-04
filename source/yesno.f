      function yesno(flag)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 7, 2004
c | Task  : Assign y or n to logical value
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      character*1 yesno
      logical     flag
c
c *********************** y or n assignment ****************************
c
c flag : distribution
c yesno: y or n function
c
      if (flag) then
        yesno='y'
      else
        yesno='n'
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
