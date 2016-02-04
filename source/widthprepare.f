      subroutine widthprepare
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 21, 2004
c | Task  : Preparation of width fluctuation corrections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ******************* Choice of width fluctuation model ****************
c
c wmode          : designator for width fluctuation model
c molprepare     : subroutine for preparation of Moldauer width
c                  fluctuation correction
c transjl,tjlav  : array for width fluctuation calculation
c tnum,tnuminc   : counter for width fluctuation calculation
c denomhf        : denominator for compound nucleus formula
c nmold,xmold,...: variables for Gauss-Legendre integration
c freedom        : number of degrees of freedom
c prodwidth      : product of widths
c numtrans       : number of transmission coefficients
c hrtwprepare    : subroutine for preparation of HRTW width fluctuation
c                  correction
c sumhrtw,.......: variables for HRTW calculation
c goeprepare     : subroutine for preparation of GOE triple integral
c                  width fluctuation correction
c agoe1,.........: variables for GOE triple integral calculation
c
c Width fluctuation models:
c wmode= 0: Hauser-Feshbach (no width fluctuations)
c wmode= 1: Moldauer
c wmode= 2: HRTW
c wmode= 3: GOE
c
c First, all width fluctuation variables that only depend on J and P and
c not on the other angular momentum quantum numbers are calculated.
c
      if (wmode.eq.1) call molprepare(transjl,tnum,denomhf,nmold,xmold,
     +  wmold,tjlav,freedom,prodwidth,numtrans,tnuminc)
      if (wmode.eq.2) call hrtwprepare(transjl,tnum,denomhf,tjlav,
     +  sumhrtw,vhrtw,whrtw,numtrans,tnuminc)
      if (wmode.eq.3) call goeprepare(transjl,tnum,denomhf,ngoep,
     +  ngoes,ngoet,xgoep,xgoes,xgoet,wgoep,wgoes,wgoet,tjlav,agoe1,
     +  agoe2,agoe3,agoe4,agoe5,agoe6,agoe7,agoe8,sgoe1,sgoe2,sgoe3,
     +  sgoe4,sgoe5,numtrans,tnuminc)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
