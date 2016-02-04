      subroutine incidentgamma
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 5, 2009
c | Task  : Incident photons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer irad,l
      real    gammaxs,fstrength,Tgamma
c
c **** Photo-absorption cross section and transmission coefficients ****
c
c lmaxinc  : maximal l-value for transmission coefficients for
c            incident channel
c gammax   : number of l-values for gamma multipolarity
c xsreacinc: reaction cross section for incident channel
c irad     : variable to indicate M(=0) or E(=1) radiation
c gammaxs  : function for gamma ray cross sections
c Einc     : incident energy in MeV
c Tgamma   : gamma transmission coefficient
c twopi    : 2.*pi
c fstrength: gamma ray strength function
c Tjlinc   : transmission coefficients as a function of radiation type
c            and l for the incident channel
c
c Note that we use the first index of Tjlinc for the radiation type,
c instead of the particle spin index.
c
      lmaxinc=gammax
      xsreacinc=gammaxs(0,0,Einc)
      do 10 irad=0,1
        do 10 l=1,gammax
          Tgamma=twopi*(Einc**(2*l+1))*fstrength(0,0,Einc,Einc,irad,l)
          Tjlinc(irad,l)=Tgamma
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
