      subroutine dwbaout(itype,type,nen1,nen2)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 7, 2006
c | Task  : Output of ECIS results for DWBA for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 ofor1,ofor2
      integer      itype,type,nen1,nen2,J,iang
c
c ********************** Write DWBA cross sections *********************
c
c ofor1,ofor2: help variable
c maxJmsd    : maximal spin for MSD calculation
c parsym     : symbol of particle
c Emsdin     : incident MSD energy
c Exmsd      : excitation energy for MSD energy grid
c Emsdout    : outgoing MSD energy
c flagddx    : flag for output of double-differential cross sections
c nanglecont : number of angles for continuum
c anglecont  : angle in degrees for continuum
c xsdw       : DWBA angular distribution as a function of incident
c              energy, outgoing energy, angular momentum and angle
c xsdwin     : DWBA cross section as a function of incident energy,
c              outgoing energy and angular momentum
c
      ofor1='(/,"  Angle",3x,n(" J =",i2,7x))'
      ofor2='(/,10x,n(" J =",i2,7x))'
      write(ofor1(17:17),'(i1)') min(maxJmsd+1,8)
      write(ofor2(8:8),'(i1)') min(maxJmsd+1,8)
      write(*,'(/," (",a1,",",a1,") DWBA cross sections ",
     +  "(mb/Sr) for E-in=",f6.2," MeV, Ex=",f6.2," MeV, E-out=",f6.2,
     +  " MeV")') parsym(itype),parsym(type),Emsdin,Exmsd,Emsdout
      if (flagddx) then
        write(*,fmt=ofor1) (J,J=0,min(maxJmsd,7))
        write(*,'()')
        do 10 iang=0,nanglecont
          write(*,'(1x,f5.1,1p,8e13.5)') anglecont(iang),
     +      (xsdw(nen1,nen2,J,iang,0),J=0,min(maxJmsd,7))
   10   continue
      endif
      write(*,fmt=ofor2) (J,J=0,min(maxJmsd,7))
      write(*,'(/," Angle",1p,8e13.5)') (xsdwin(nen1,nen2,J,0),
     +  J=0,min(maxJmsd,7))
      write(*,'(" integr.")')
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
