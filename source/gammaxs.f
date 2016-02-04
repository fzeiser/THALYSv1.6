      function gammaxs(Zcomp,Ncomp,Egamma)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 5, 2009
c | Task  : Gamma ray cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,irad,l
      real    gammaxs,Egamma,xsgdr,fstrength,levinger,freedeut,fpauli,
     +        xsqd
c
c ************** Calculate photo-absorption cross section **************
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c Einc     : incident energy in MeV
c Egamma   : gamma energy
c xsgdr    : photo-absorption cross section from GDR part
c fstrength: gamma ray strength function
c kgr      : constant for gamma-ray strength function
c
c 1. GDR part
c
      xsgdr=0.
      do 10 irad=0,1
        do 10 l=1,gammax
          xsgdr=xsgdr+fstrength(Zcomp,Ncomp,Einc,Egamma,irad,l)/
     +      kgr(Zcomp,Ncomp,irad,l)*Egamma
   10 continue
c
c 2. QD part
c
c levinger: Levinger parameter
c freedeut: free deuteron cross section
c fpauli  : Pauli-blocking function of Chadwick
c xsqd    : photo-absorption cross section from QD part
c Ntarget : neutron number of target nucleus
c Ztarget : charge number of target nucleus
c Atarget : mass number of target nucleus
c
      levinger=6.5
      if (Egamma.gt.2.224) then
        freedeut=61.2/(Egamma**3)*(Egamma-2.224)**1.5
      else
        freedeut=0.
      endif
      if (Egamma.le.140.) then
        if (Egamma.ge.20.) then
          fpauli=8.3714e-2-9.8343e-3*Egamma+4.1222e-4*Egamma*Egamma-
     +      3.4762e-6*(Egamma**3)+9.3537e-9*(Egamma**4)
        else
          fpauli=exp(-73.3/Egamma)
        endif
      else
        fpauli=exp(-24.2348/Egamma)
      endif
      xsqd=levinger*Ntarget*Ztarget/real(Atarget)*freedeut*fpauli
c
c Total absorption cross section
c
      gammaxs=xsgdr+xsqd
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
