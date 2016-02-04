      subroutine opticaldeut(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 11, 2013
c | Task  : Other optical potential for deuterons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,i,Z,A,N
      real    eopt,A13,mu,summu,beta,asym
c
c Alternative options for deuteron OMP
c
c deuteronomp=1 --> Watanabe potential
c deuteronomp=2 --> Daehnick potential
c deuteronomp=3 --> Bojowald potential
c deuteronomp=4 --> Han potential
c deuteronomp=5 --> Haixia An potential
c
c Overwrite some of the previous values.
c
c ompadjust    : subroutine for local optical model parameter adjustment
c eopt         : incident energy
c ZZ,Z         : charge number of residual nucleus
c AA,A         : mass number of residual nucleus
c N            : neutron number of residual nucleus
c onethird     : 1/3
c asym         : asymmetry parameter
c magic        : magic numbers
c beta,mu      : help variables
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c
      call ompadjust(eopt,3)
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      N=A-Z
      A13=real(A)**onethird
      asym=(N-Z)/real(A)
c
c deuteronomp=2 --> Global OMP of W.W. Daehnick, J.D. Childs, Z. Vrcelj,
c Phys. Rev C21, 2253 (1980): relativistic case
c
      if (deuteronomp.eq.2) then
        v=88.0+0.88*Z/A13-0.283*eopt
        rv=1.17
        av=0.717+0.0012*eopt
        beta=-(0.01*eopt)**2
        w=(12+0.031*eopt)*(1.-exp(beta))
        rw=1.376-0.01*sqrt(eopt)
        summu=0.
        do i=1,8
          mu=(0.5*(magic(i)-N))**2
          summu=summu+exp(-mu)
        enddo
        aw=0.52+0.07*A13-0.04*summu
        vd=0.
        rvd=rw
        avd=aw
        wd=(12.+0.031*eopt)*exp(beta)
        rwd=rvd
        awd=avd
        vso=5.0
        rvso=1.04
        avso=0.60
        wso=0.37*A13-0.03*eopt
        rwso=0.80
        awso=0.25
        rc=1.30
      endif
c
c deuteronomp=3 --> Global OMP of J. Bojowald et al,
c Phys. Rev. C38, 1153 (1988).
c
      if (deuteronomp.eq.3) then
        v=81.32-0.24*eopt+1.43*Z/A13
        rv=1.18
        av=0.636+0.035*A13
        w=max(0.,0.132*(eopt-45.))
        rw=1.27
        aw=0.768+0.021*A13
        vd=0.
        rvd=rw
        avd=aw
        wd=max(0.,7.80+1.04*A13-0.712*w)
        rwd=rvd
        awd=avd
        vso=6.
        rvso=0.78+0.038*A13
        avso=rvso
        wso=0.
        rwso=rvso
        awso=avso
        rc=1.30
      endif
c
c deuteronomp=4 --> Global OMP of Y. Han, Y. Shi and Q. Shen,
c Phys. Rev. C74, 044615 (2006).
c
      if (deuteronomp.eq.4) then
        v=82.18-0.148*eopt-0.000886*eopt*eopt-34.811*asym+1.058*Z/A13
        rv=1.174
        av=0.809
        w=max(0.,-4.916+0.0555*eopt+4.42e-5*eopt*eopt+35.*asym)
        rw=1.563
        aw=0.700+0.045*A13
        vd=0.
        rvd=1.328
        avd=0.465+0.045*A13
        wd=20.968-0.0794*eopt-43.398*asym
        rwd=rvd
        awd=avd
        vso=3.703
        rvso=1.234
        avso=0.813
        wso=-0.206
        rwso=rvso
        awso=avso
        rc=1.698
      endif
c
c deuteronomp=5 --> Global OMP of Haixia An and Chonghai Cai,
c Phys. Rev. C73, 054605 (2006).
c
      if (deuteronomp.eq.5) then
        v=91.85-0.249*eopt-0.000116*eopt*eopt+0.642*Z/A13
        rv=1.152-0.00776/A13
        av=0.719+0.0126*A13
        w=1.104+0.0622*eopt
        rw=1.305+0.0997/A13
        aw=0.855-0.100*A13
        vd=0.
        rvd=1.334+0.152/A13
        avd=0.531+0.062*A13
        wd=10.83-0.0306*eopt
        rwd=rvd
        awd=avd
        vso=3.557
        rvso=0.972
        avso=1.011
        wso=0.
        rwso=rvso
        awso=avso
        rc=1.303
      endif
c
c Possible adjustment of parameters
c
c Fv1,....     : adjustable factors for OMP (default 1.)
      v=Fv1*v
      rv=Frv*rv
      av=Fav*av
      w=Fw1*w
      rw=Frw*rw
      aw=Faw*aw
      vd=Fd1*vd
      rvd=Frvd*rvd
      avd=Favd*avd
      wd=Fd1*wd
      rwd=Frwd*rwd
      awd=Fawd*awd
      vso=Fvso1*vso
      rvso=Frvso*rvso
      avso=Favso*avso
      wso=Fwso1*wso
      rwso=Frwso*rwso
      awso=Fawso*awso
      rc=Frc*rc
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
