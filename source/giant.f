      subroutine giant
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 12, 2004
c | Task  : Giant resonance contribution
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,i,nen,iang,Zix,Nix
      real    wscale,grwidth,fac1,fac2,sumgauss,edist,gauss(0:numen),
     +        weight,diswidth
c
c ************* Smearing of giant resonances into spectra **************
c
c wscale         : scaling factor for giant resonance to Gaussian width
c l              : l-value
c betagr         : deformation parameter for giant resonance
c eoutgr         : outgoing energy
c k0             : index of incident particle
c eninccm        : center-of-mass incident energy in MeV
c Egrcoll        : energy of giant resonance
c Ggrcoll,grwidth: width of giant resonance
c fac1,fac2,edist: help variables
c sqrttwopi      : sqrt(2.*pi)
c sumgauss       : sum over Gaussians (for normalization)
c ebegin         : first energy point of energy grid
c eend           : last energy point of energy grid
c gauss          : Gaussian contribution
c edist          : help variable
c egrid          : outgoing energy grid
c weight         : Gaussian weight
c xsgrstate      : smoothed giant resonance cross section per state
c xsgrcoll       : giant resonance cross section
c xsgr           : smoothed giant resonance cross section
c flagddx        : flag for output of double-differential cross sections
c nanglecont     : number of angles for continuum
c xsgrad         : smoothed giant resonance angular distribution
c grcollad       : giant resonance angular distribution
c xsgrtot        : total smoothed giant resonance cross section
c xsgrsum        : sum over giant resonance cross sections
c
      wscale=0.42
      do 10 l=0,3
        do 20 i=1,2
          if (betagr(l,i).eq.0.) goto 20
          eoutgr(k0,l,i)=eninccm-Egrcoll(l,i)
          grwidth=Ggrcoll(l,i)*wscale
          fac1=1./(grwidth*sqrttwopi)
          fac2=1./(2.*grwidth**2)
          sumgauss=0.
          do 30 nen=ebegin(k0),eend(k0)
            gauss(nen)=0.
            edist=abs(eoutgr(k0,l,i)-egrid(nen))
            if (edist.gt.5.*grwidth) goto 30
            gauss(nen)=fac1*exp(-(edist**2)*fac2)
            sumgauss=sumgauss+gauss(nen)
   30     continue
          if (sumgauss.ne.0.) then
            do 40 nen=ebegin(k0),eend(k0)
              weight=gauss(nen)/sumgauss
              xsgrstate(k0,l,i,nen)=weight*xsgrcoll(k0,l,i)
              xsgr(k0,nen)=xsgr(k0,nen)+xsgrstate(k0,l,i,nen)
              if (flagddx) then
                do 50 iang=0,nanglecont
                  xsgrad(k0,nen,iang)=xsgrad(k0,nen,iang)+
     +              weight*grcollad(k0,l,i,iang)
   50           continue
              endif
   40       continue
          endif
          xsgrtot(k0)=xsgrtot(k0)+xsgrcoll(k0,l,i)
   20   continue
   10 continue
      xsgrsum=xsgrtot(k0)
c
c *********** Other collective contributions to the continuum **********
c
c xscollconttot: total collective cross section in the continuum
c Zindex,Zix   : charge number index for residual nucleus
c Nindex,Nix   : neutron number index for residual nucleus
c Nlast        : last discrete level
c numlev2      : maximum number of levels
c deform       : deformation parameter
c eoutdis      : outgoing energy of discrete state reaction
c diswidth     : width of discrete level peak
c elwidth      : width of elastic peak in MeV
c xscollcont   : collective cross section in the continuum
c xsdirdisc    : direct cross section for discrete state
c deltaE       : energy bin around outgoing energies
c collcontad   : collective angular distribution in the continuum
c directad     : direct angular distribution
c
      if (xscollconttot.eq.0.) return
      xsgrtot(k0)=xsgrtot(k0)+xscollconttot
      xsgrsum=xsgrsum+xscollconttot
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      do 110 i=Nlast(Zix,Nix,0)+1,numlev2
        if (deform(Zix,Nix,i).eq.0.) goto 110
        if (eoutdis(k0,i).le.0.) goto 110
        diswidth=elwidth*(eoutdis(k0,i)/eninccm)**1.5
        fac1=1./(diswidth*sqrttwopi)
        fac2=1./(2.*diswidth**2)
        sumgauss=0.
        do 120 nen=ebegin(k0),eend(k0)
          gauss(nen)=0.
          edist=abs(egrid(nen)-eoutdis(k0,i))
          if (edist.gt.5.*diswidth) goto 120
          gauss(nen)=fac1*exp(-(edist**2)*fac2)
          sumgauss=sumgauss+gauss(nen)
  120   continue
        if (sumgauss.ne.0.) then
          do 130 nen=ebegin(k0),eend(k0)
            weight=gauss(nen)/sumgauss/deltaE(nen)
            xscollcont(nen)=xscollcont(nen)+weight*xsdirdisc(k0,i)
            if (flagddx) then
              do 140 iang=0,nanglecont
                collcontad(nen,iang)=collcontad(nen,iang)+
     +            weight*directad(k0,i,iang)
  140         continue
            endif
  130     continue
        endif
  110 continue
c
c Add collective contribution to giant resonance results
c
      do 210 nen=ebegin(k0),eend(k0)
        xsgr(k0,nen)=xsgr(k0,nen)+xscollcont(nen)
        if (flagddx) then
          do 220 iang=0,nanglecont
            xsgrad(k0,nen,iang)=xsgrad(k0,nen,iang)+collcontad(nen,iang)
  220     continue
        endif
  210 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
