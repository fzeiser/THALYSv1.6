      subroutine widthfluc(ielas)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 22, 2004
c | Task  : Width fluctuation corrections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ielas
c
c ******************* Choice of width fluctuation model ****************
c
c ielas          : designator for elastic channel
c tnumo,tnumi    : counters for width fluctuation calculation
c numtrans       : number of transmission coefficients
c numhill        : maximum number of Hill-Wheeler points
c wmode          : designator for width fluctuation model
c moldauer       : subroutine for Moldauer width fluctuation correction
c tnum           : total number of transmission coefficients
c denomhf        : denominator for compound nucleus formula
c nmold,xmold    : variables for Gauss-Legendre integration
c tjlav,transjl  : array for width fluctuation calculation
c freedom        : number of degrees of freedom
c prodwidth      : product of widths
c Wab            : width fluctuation factor
c hrtw           : subroutine for HRTW width fluctuation correction
c sumhrtw,.......: variables for HRTW calculation
c goe            : subroutine for GOE triple integral width fluctuation
c                  correction
c agoe1,.........: variables for GOE triple integral calculation
c
c Width fluctuation models:
c wmode= 0: Hauser-Feshbach (no width fluctuations)
c wmode= 1: Moldauer
c wmode= 2: HRTW
c wmode= 3: GOE
c
      if (tnumo.gt.numtrans-numhill-2) return
      if (wmode.eq.1) call moldauer(tnum,tnumi,tnumo,denomhf,nmold,
     +  xmold,tjlav,freedom,prodwidth,Wab,numtrans,ielas)
      if (wmode.eq.2) call hrtw(transjl,tnumi,tnumo,denomhf,sumhrtw,
     +  vhrtw,whrtw,Wab,numtrans,ielas)
      if (wmode.eq.3) call goe(transjl,tnum,tnumi,tnumo,denomhf,ngoep,
     +  ngoes,ngoet,tjlav,agoe1,agoe2,agoe3,agoe4,agoe5,agoe6,agoe7,
     +  agoe8,sgoe1,sgoe2,sgoe3,sgoe4,sgoe5,Wab,numtrans,ielas)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
