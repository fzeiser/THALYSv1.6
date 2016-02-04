      subroutine onecontinuumA(itype,type)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 1, 2008
c | Task  : Unnormalized one-step direct cross sections for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer itype,type,h,p,Zix,Nix,nen1,nen2,J,iang
      real    gs,QQ,ignatyuk,rJ,omegaJ(0:numJmsd),omega,total,xs
c
c ************* Calculate continuum one-step cross sections ************
c
c h         : hole number
c p         : particle number
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c g,gs      : single-particle level density parameter
c QQ        : Q-value
c S         : separation energy per particle
c msdbins2  : number of energy points for MSD calculation
c Emsdin    : incident MSD energy
c specmass  : specific mass for target nucleus
c Emsd      : MSD energy grid
c Emsdout   : outgoing MSD energy
c Exmsd     : excitation energy for MSD energy grid
c flaggshell: flag for energy dependence of single particle level
c             density parameter g
c ignatyuk  : function for energy dependent level density parameter a
c alev      : level density parameter
c maxJmsd   : maximal spin for MSD calculation
c rJ        : help variable
c omegaJ    : help variable
c omega     : particle-hole state density
c total,xs  : help variables
c xsdwin    : DWBA cross section as a function of incident energy,
c             outgoing energy and angular momentum
c xscont1   : continuum one-step direct cross section for MSD
c             (unnormalized)
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c xsdw      : DWBA angular distribution as a function of incident
c             energy, outgoing energy, angular momentum and angle
c xscontad1 : continuum one-step direct angular distribution for MSD
c             (unnormalized)
c
      h=1
      p=1
      Zix=Zindex(0,0,type)
      Nix=Nindex(0,0,type)
      gs=g(Zix,Nix)
      QQ=S(0,0,itype)-S(0,0,type)
      do 10 nen1=0,msdbins2
        Emsdin=specmass(Zix,Nix,itype)*Emsd(nen1)
        do 20 nen2=nen1,msdbins2
          Emsdout=Emsd(nen2)
          Exmsd=Emsdin-Emsdout+QQ
          if (Exmsd.lt.0..or.(Exmsd+0.1).ge.Emsdin) goto 20
          if (flaggshell) gs=g(Zix,Nix)*ignatyuk(Zix,Nix,Exmsd,0)/
     +      alev(Zix,Nix)
          do 30 J=0,maxJmsd
            rJ=real(J)
            omegaJ(J)=omega(Zix,Nix,p,h,gs,Exmsd,rJ)
   30     continue
          total=0.
          do 40 J=0,maxJmsd
            xs=xsdwin(nen1,nen2,J,0)
            total=total+omegaJ(J)*xs
   40     continue
          xscont1(itype,type,nen1,nen2)=total
          if (flagddx) then
            do 50 iang=0,nanglecont
              total=0.
              do 60 J=0,maxJmsd
                xs=xsdw(nen1,nen2,J,iang,0)
                total=total+omegaJ(J)*xs
   60         continue
              xscontad1(itype,type,nen1,nen2,iang)=total
   50       continue
          endif
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
