      subroutine incidentnorm
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 11, 2004
c | Task  : Normalization of reaction cross sections and transmission
c |         coefficients for incident channel
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,ispin
      real    enuc,xs,tripathi,norm
c
c ************ Normalization with semi-empirical results ***************
c
c flagsys   : flag for reaction cross section from systematics
c k0        : index of incident particle
c xsoptinc  : optical model reaction cross section for incident channel
c enuc      : incident energy in MeV per nucleon
c Einc      : incident energy in MeV
c parA      : mass number of particle
c xs        : help variable
c tripathi  : function for semi-empirical reaction cross section of
c             Tripathi et al.
c parZ      : charge number of particle
c Ztarget   : charge number of target nucleus
c Atarget   : mass number of target nucleus
c norm      : normalization factor
c threshnorm: normalization factor at trheshold
c xsreacinc : reaction cross section for incident channel
c xselasinc : total elastic cross section (neutrons only) for incident
c             channel
c l         : orbital angular momentum
c numl      : maximal number of l-values in TALYS
c Tlinc     : transmission coefficients as a function of l for the
c             incident channel, averaged over spin
c ispin     : spin index
c Tjlinc    : transmission coefficients as a function of spin and l for
c             the incident channel
c
c The normalization is only performed if the option for semi-empirical
c reaction cross sections is enabled.
c
      if (.not.flagsys(k0)) return
      if (xsoptinc.eq.0.) return
      enuc=Einc/parA(k0)
      xs=tripathi(parZ(k0),parA(k0),Ztarget,Atarget,enuc)
      if (xs.eq.0.) then
        norm=threshnorm(k0)
        xs=xsoptinc*threshnorm(k0)
      else
        norm=xs/xsoptinc
      endif
      xsreacinc=xs
      if (k0.eq.1) xselasinc=xselasinc+xsoptinc-xs
      do 10 l=0,numl
        Tlinc(l)=Tlinc(l)*norm
        do 20 ispin=-1,1
          Tjlinc(ispin,l)=Tjlinc(ispin,l)*norm
   20   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
