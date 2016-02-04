      subroutine spectra
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 5, 2011
c | Task  : Creation of spectra
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,NL,nend,nen,iang,fine,nhigh,nenout,i,iangdisc
      real    Elast,Ares,convfac1,convfac2,convfac3,Ehigh,Eout,Ea,Eb,
     +        factor,diswidth,fac1,fac2,gauss
c
c ********** Add smoothed discrete cross sections to spectra ***********
c
c xsdisc            : total cross section for discrete state
c k0                : index of incident particle
c Ltarget           : excited level of target
c xselastot         : total elastic cross section (shape + compound)
c parskip           : logical to skip outgoing particle
c xsparticle        : total particle production cross section
c Nlast,NL          : last discrete level
c parZ              : charge number of particle
c parN              : neutron number of particle
c flagadd           : flag for addition of discrete states to spectra
c flagnatural       : flag for calculation of natural element
c Elast             : last outgoing energy for smoothed discrete
c                     contribution
c Einc              : incident energy in MeV
c elwidth           : width of elastic peak in MeV
c eoutdis           : outgoing energy of discrete state reaction
c locate            : subroutine to find value in ordered table
c egrid,espec       : outgoing energy grid
c ebegin            : first energy point of energy grid
c nend              : help variable
c eend,eendout      : last energy point of energy grid
c flagEchannel      : flag for channel energy for emission spectrum
c Ares              : mass number of residual nucleus
c Ainit             : mass number of initial compound nucleus
c parA              : mass number of particle
c convfac1,.....    : conversion factor for reference system
c xsdiscout         : total smoothed cross section for discrete state
c xsgr              : smoothed giant resonance cross section
c xspreeq,xspreeqout: preequilibrium cross section per particle type
c                     and outgoing energy
c xsmpreeq          : multiple pre-equilibrium emission spectrum
c xsmpreeqout       : multiple pre-equilibrium emission spectrum
c xscomp,xscompout  : compound emission spectrum
c flagddx           : flag for output of double-differential cross
c                     sections
c nanglecont        : number of angles for continuum
c xsdiscoutad       : smoothed angular distribution for discrete state
c xsgrad            : smoothed giant resonance angular distribution
c xspreeqoutad      : preequilibrium angular distribution per particle
c                     type
c xsmpreeqoutad     : multiple preequilibrium angular distribution
c xscompoutad       : compound emission angular distribution
c xscompad          : compound emission angular distribution
c flagendf          : flag for information for ENDF-6 file
c fine              : refinement factor for high-energy tail of spectrum
c Ehigh,nhigh       : help variables
c segment           : number of segments to divide emission energy grid
c numendisc         : number of discrete outgoing energies
c Eout              : outgoing energy
c Ea,Eb,factor      : help variables
c flagaddel         : flag for addition of elastic peak to spectra
c diswidth          : width of discrete level peak
c eninccm           : center-of-mass incident energy in MeV
c fac1,fac2         : help variables
c sqrttwopi         : sqrt(2.*pi)
c gauss             : Gaussian contribution
c iangdisc          : counter for discrete angle
c nangle            : number of angles
c directad          : direct angular distribution
c discad            : discrete state angular distribution
c
      xsdisc(k0,Ltarget)=xselastot
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (xsparticle(type).eq.0.) goto 10
        NL=Nlast(parZ(type),parN(type),0)
        if (flagadd.or.type.ge.2) then
          if (flagnatural) then
            Elast=Einc-6.-elwidth
          else
            Elast=eoutdis(type,NL)-elwidth
          endif
          Elast=max(Elast,0.)
          call locate(egrid,ebegin(type),eend(type),Elast,nend)
        else
          nend=eend(type)
          eendout(type)=eend(type)
        endif
        convfac1=1.
        convfac2=0.
        convfac3=0.
        if (flagEchannel) then
          Ares=real(Ainit-parA(type))
          convfac1=Ares/Ainit
          convfac2=real(parA(type))/Ainit*real(parA(k0))/Ainit*Einc
          if (convfac1.gt.0..and.convfac2.gt.0.)
     +      convfac3=2.*sqrt(convfac1*convfac2)
        endif
        do 20 nen=ebegin(type),nend
          espec(type,nen)=convfac1*egrid(nen)+convfac2+
     +      convfac3*sqrt(egrid(nen))
          xsdiscout(type,nen)=xsgr(type,nen)
          xspreeqout(type,nen)=xspreeq(type,nen)
          xsmpreeqout(type,nen)=xsmpreeq(type,nen)
          xscompout(type,nen)=xscomp(type,nen)
          if (flagddx) then
            do 30 iang=0,nanglecont
              xsdiscoutad(type,nen,iang)=xsgrad(type,nen,iang)
              xspreeqoutad(type,nen,iang)=xspreeqad(type,nen,iang)
              xsmpreeqoutad(type,nen,iang)=xsmpreeqad(type,nen,iang)
              xscompoutad(type,nen,iang)=xscompad(type,nen,iang)
   30       continue
          endif
   20   continue
        if (.not.(flagadd.or.type.ge.2)) goto 100
        if (flagendf) then
          fine=2
        else
          fine=10
        endif
        Ehigh=egrid(eend(type))-egrid(nend)
        nhigh=int(fine*segment*Ehigh)
        nhigh=min(nhigh,numendisc-8*fine)
        eendout(type)=nend+nhigh+8*fine
        do 40 nenout=nend+1,eendout(type)
          espec(type,nenout)=egrid(nend)+(nenout-nend)/
     +      (segment*real(fine))
          Eout=espec(type,nenout)
          espec(type,nenout)=convfac1*Eout+convfac2+convfac3*sqrt(Eout)
          call locate(egrid,nend,eend(type),Eout,nen)
          Ea=Eout-egrid(nen)
          Eb=egrid(nen+1)-egrid(nen)
          if (Ea.lt.Eb) then
            factor=Ea/Eb
          else
            factor=1.
          endif
          xsdiscout(type,nenout)=xsgr(type,nen)+
     +      factor*(xsgr(type,nen+1)-xsgr(type,nen))
          xspreeqout(type,nenout)=xspreeq(type,nen)+
     +      factor*(xspreeq(type,nen+1)-xspreeq(type,nen))
          xsmpreeqout(type,nenout)=xsmpreeq(type,nen)+
     +      factor*(xsmpreeq(type,nen+1)-xsmpreeq(type,nen))
          xscompout(type,nenout)=xscomp(type,nen)+
     +      factor*(xscomp(type,nen+1)-xscomp(type,nen))
          if (flagddx) then
            do 50 iang=0,nanglecont
              xsdiscoutad(type,nenout,iang)=xsgrad(type,nen,iang)+
     +          factor*(xsgrad(type,nen+1,iang)-xsgrad(type,nen,iang))
              xspreeqoutad(type,nenout,iang)=xspreeqad(type,nen,iang)+
     +          factor*(xspreeqad(type,nen+1,iang)-
     +          xspreeqad(type,nen,iang))
              xsmpreeqoutad(type,nenout,iang)=xsmpreeqad(type,nen,iang)+
     +          factor*(xsmpreeqad(type,nen+1,iang)-
     +          xsmpreeqad(type,nen,iang))
              xscompoutad(type,nenout,iang)=xscompad(type,nen,iang)+
     +          factor*(xscompad(type,nen+1,iang)-
     +          xscompad(type,nen,iang))
   50       continue
          endif
          do 60 i=0,NL
            if (i.eq.Ltarget.and.type.eq.k0.and.k0.gt.1) goto 60
            if (i.eq.Ltarget.and..not.flagaddel) goto 60
            if (eoutdis(type,i).le.0.) goto 60
            diswidth=elwidth*(eoutdis(type,i)/eninccm)**1.5
            fac1=1./(diswidth*sqrttwopi)
            fac2=1./(2.*diswidth**2)
            gauss=fac1*exp(-(Eout-eoutdis(type,i))**2*fac2)
            xsdiscout(type,nenout)=xsdiscout(type,nenout)+
     +        gauss*xsdisc(type,i)
            if (flagddx) then
              do 70 iang=0,nanglecont
                iangdisc=real(iang*nangle)/nanglecont
                  xsdiscoutad(type,nenout,iang)=
     +            xsdiscoutad(type,nenout,iang)+
     +            gauss*discad(type,i,iangdisc)
   70         continue
            endif
   60     continue
   40   continue
c
c ************************ Create total spectra ************************
c
c xssumout  : cross section summed over mechanisms
c preeqratio: pre-equilibrium ratio
c xssumoutad: angular distribution summed over mechanisms
c
  100   do 110 nen=ebegin(type),eendout(type)
          xssumout(type,nen)=xspreeqout(type,nen)+xsmpreeqout(type,nen)+
     +      xscompout(type,nen)
          if (xssumout(type,nen).ne.0.) then
            preeqratio(type,nen)=(xspreeqout(type,nen)+
     +        xsmpreeqout(type,nen))/xssumout(type,nen)
          else
            preeqratio(type,nen)=0.
          endif
          xssumout(type,nen)=xssumout(type,nen)+xsdiscout(type,nen)
          if (flagddx) then
            do 120 iang=0,nanglecont
              xssumoutad(type,nen,iang)=xsdiscoutad(type,nen,iang)+
     +          xspreeqoutad(type,nen,iang)+
     +          xsmpreeqoutad(type,nen,iang)+xscompoutad(type,nen,iang)
  120         continue
          endif
  110   continue
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
