      function surface(type,elab)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 7, 2004
c | Task  : Surface effects in exciton model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type
      real    surface,elab
c
c ************************* Kalbach formula ****************************
c
c surface : well depth for first hole
c elab    : incident energy
c Ainit   : mass number of initial compound nucleus
c onethird: 1/3
c
      if (type.eq.1) then
        surface=12.+26.*elab**4/(elab**4+(245./(Ainit**onethird))**4)
      else
        surface=22.+16.*elab**4/(elab**4+(450./(Ainit**onethird))**4)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
