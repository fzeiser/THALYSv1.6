      subroutine inversenorm(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 19, 2009
c | Task  : Normalization and extrapolation of reaction cross sections
c |         and transmission coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,type,Z,A,nen,l,ispin
      real    xsprev,enuc,xs,tripathi,xstripathi(numen),norm
c
c ************ Normalization with semi-empirical results ***************
c
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c parskip   : logical to skip outgoing particle
c flagsys   : flag for reaction cross section from systematics
c ZZ,Z      : charge number of residual nucleus
c AA,A      : mass number of residual nucleus
c xsprev    : help variable
c threshnorm: normalization factor at trheshold
c ebegin    : first energy point of energy grid
c eendmax   : last energy point of energy grid for maximum incident
c             energy
c xsopt     : optical model reaction cross section
c enuc      : incident energy in MeV per nucleon
c egrid     : outgoing energy grid
c parA      : mass number of particle
c tripathi  : function for semi-empirical reaction cross section of
c             Tripathi et al.
c xstripathi: help variable
c parZ      : charge number of particle
c norm      : normalization factor
c xsreac    : reaction cross section
c xselas    : total elastic cross section (neutrons only)
c lcc,l     : orbital angular momentum
c numl      : maximal number of l-values in TALYS
c Tl        : transmission coefficients as a function of particle type,
c             energy and l-value (averaged over spin)
c ispin     : spin index
c Tjl       : transmission coefficients as a function of particle type,
c             energy, spin and l-value
c
c The normalization is only performed if the option for semi-empirical
c reaction cross sections is enabled.
c The semi-empirical results have a too sharp cutoff at low energies.
c Therefore, for the lowest energies the optical model results are
c renormalized with the ratio at the threshold.
c
      do 10 type=1,6
        if (parskip(type)) goto 10
        if (.not.flagsys(type)) goto 10
        Z=ZZ(Zcomp,Ncomp,type)
        A=AA(Zcomp,Ncomp,type)
        xsprev=0.
        do 20 nen=ebegin(type),eendmax(type)
          if (xsopt(type,nen).eq.0.) goto 20
          enuc=egrid(nen)/parA(type)
          xs=tripathi(parZ(type),parA(type),Z,A,enuc)
          xstripathi(nen)=xs
          if (xsprev.eq.0.and.xs.ne.0.) threshnorm(type)=
     +      xs/xsopt(type,nen)
          xsprev=xs
   20   continue
        do 30 nen=ebegin(type),eendmax(type)
          if (xsopt(type,nen).eq.0.) goto 30
          xs=xstripathi(nen)
          if (xs.eq.0.) then
            norm=threshnorm(type)
            xs=xsopt(type,nen)*threshnorm(type)
          else
            norm=xs/xsopt(type,nen)
          endif
          xsreac(type,nen)=xs
          if (type.eq.1) xselas(type,nen)=xselas(type,nen)+
     +      xsopt(type,nen)-xs
          do 40 l=0,numl
            Tl(type,nen,l)=Tl(type,nen,l)*norm
            do 50 ispin=-1,1
              Tjl(type,nen,ispin,l)=Tjl(type,nen,ispin,l)*norm
   50       continue
   40     continue
   30   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
