      subroutine interangle
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 10, 2004
c | Task  : Intermediate angles by addition theorem for MSD model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer iangout,iangint,iphi,iang
      real    angstep,angout,angint,cosine,sine,phi,cosgam,gam
c
c ************************ Addition theorem ****************************
c
c angstep       : angle step in degrees
c nanglecont    : number of angles for continuum
c iangout       : outgoing angle index
c angout        : outgoing angle
c deg2rad       : conversion factor for degrees to radians
c anglecont     : angle in degrees for continuum
c iangint       : intermediate angle index
c angint        : intermediate angle
c cosine,sine   : help variables
c twopi         : 2.*pi
c rad2deg       : conversion factor for radians to degrees
c phi,cosgam,gam: help variables
c nangleint     : number of possibilities to link intermediate angle to
c                 final angle
c
c This is necessary to transform the ingoing angle of the second step
c and link it with the outgoing angle of the first step.
c
      angstep=180./nanglecont
      do 10 iangout=0,nanglecont
        angout=anglecont(iangout)*deg2rad
        do 20 iangint=0,nanglecont
          angint=anglecont(iangint)*deg2rad
          cosine=cos(angout)*cos(angint)
          sine=sin(angout)*sin(angint)
          do 30 iang=0,numangcont
            nangleint(iangout,iangint,iang)=0
   30     continue
          do 40 iphi=0,2*nanglecont-1
            phi=real(iphi)*0.5/nanglecont*twopi
            cosgam=cosine+cos(phi)*sine
            if (abs(cosgam).gt.1.0) cosgam=1.0
            gam=acos(cosgam)*rad2deg
            iang=int((gam+angstep/2.)/angstep)
            nangleint(iangout,iangint,iang)=
     +        nangleint(iangout,iangint,iang)+1
   40     continue
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
