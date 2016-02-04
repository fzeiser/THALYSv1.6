      subroutine dwbaint
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : Interpolate DWBA cross sections for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nen1,J,iang,nen2
c
c ******************* Interpolate DWBA cross sections ******************
c
c msdbins2  : number of energy points for MSD calculation
c maxJmsd   : maximal spin for MSD calculation
c xsdwin    : DWBA cross section as a function of incident energy,
c             outgoing energy and angular momentum
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c xsdw      : DWBA angular distribution as a function of incident
c             energy, outgoing energy, angular momentum and angle
c
      do 10 nen1=0,msdbins2,2
        do 20 J=0,maxJmsd
          xsdwin(nen1,nen1,J,0)=xsdwin(nen1,nen1+2,J,0)
          xsdwin(nen1+1,nen1+1,J,0)=xsdwin(nen1,nen1+2,J,0)
          if (flagddx) then
            do 30 iang=0,nanglecont
              xsdw(nen1+1,nen1+1,J,iang,0)=xsdw(nen1,nen1+2,J,iang,0)
   30       continue
          endif
   20   continue
        do 40 nen2=nen1+1,msdbins2-1,2
          do 40 J=0,maxJmsd
            xsdwin(nen1,nen2,J,0)=0.5*(xsdwin(nen1,nen2-1,J,0)+
     +        xsdwin(nen1,nen2+1,J,0))
            if (flagddx) then
              do 50 iang=0,nanglecont
                xsdw(nen1,nen2,J,iang,0)=0.5*(xsdw(nen1,nen2-1,J,iang,0)
     +            +xsdw(nen1,nen2+1,J,iang,0))
   50         continue
            endif
   40   continue
   10 continue
      do 60 nen1=1,msdbins2-1,2
        do 60 nen2=nen1+1,msdbins2,2
          do 60 J=0,maxJmsd
            xsdwin(nen1,nen2,J,0)=0.5*(xsdwin(nen1-1,nen2,J,0)+
     +        xsdwin(nen1+1,nen2,J,0))
            if (flagddx) then
              do 70 iang=0,nanglecont
                xsdw(nen1,nen2,J,iang,0)=0.5*(xsdw(nen1-1,nen2,J,iang,0)
     +            +xsdw(nen1+1,nen2,J,iang,0))
   70         continue
            endif
   60 continue
      do 80 nen1=1,msdbins2-1,2
        do 80 nen2=nen1+2,msdbins2-1,2
          do 80 J=0,maxJmsd
            xsdwin(nen1,nen2,J,0)=0.5*(xsdwin(nen1-1,nen2-1,J,0)+
     +        xsdwin(nen1+1,nen2+1,J,0))
            if (flagddx) then
              do 90 iang=0,nanglecont
                xsdw(nen1,nen2,J,iang,0)=0.5*
     +            (xsdw(nen1-1,nen2-1,J,iang,0)+
     +            xsdw(nen1+1,nen2+1,J,iang,0))
   90         continue
            endif
   80 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
