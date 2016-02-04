      subroutine comprecoil(Zcomp,Ncomp,nex,type,nexout,nenbeg,nenend)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 12, 2004
c | Task  : Recoils from compound decay
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,nex,type,nexout,nenbeg,nenend,nen,Zix,Nix
      real    Eout,ang,xsgridad(numen,0:numangcont)
      integer   iymax,iymaxp1
      integer   ix,iy
      integer   ierec,iarec,ierecbeg,ierecend
      real      compmass,ratio
      real      fracCM,fracCMloc(0:numenrec)
      real      vcmloc(0:numenrec)
      real      dEejcm
      real      pejcm1,pejcm2
      integer   iang
      real      fluxCM
      real      kalbach
      real      labs1,labs2,labsurf
      real      scovej1(1:numen2,0:2*numangcont+1)
      real      scovej2(1:numen2,0:2*numangcont+1)
      real      scovrec1(0:numenrec,0:2*numangrec+1)
      real      scovrec2(0:numenrec,0:2*numangrec+1)
      integer   xlim1(2),ylim1(2),xlim2(2),ylim2(2)
      real      ddxLAB
      integer   iymod,iys
      real      sumddxr
      real      ddxrecadd,ddxejadd
      real      surfbin,fluxadd
      real      vcmav
c
c *** Calculate recoil and light particle spectra in the LAB frame *****
c
c We are decaying from a compound nucleus bin with index nex and
c excitation energies characteristics
c        Exinc=Ex(Zcomp,Ncomp,nex) : middle of the bin
c        Ex0plus=Exinc+0.5*dExinc  : upper limit of the bin
c        Ex0min=Exinc-0.5*dExinc   : lower limit of the bin
c with a given ejectile (i.e. we are still within the type loop) in a
c residual nucleus bin nexout characterised by an excitation energy
c Exout=Ex(Zix,Nix,nexout) and the bin limits Ex1min
c (lower boundary of residual bin) and Ex1plus
c (upper boundary of residual bin). The corresponding ejectile emission
c energies have been deduced and are given by the egrid(ieject) values
c were ieject is between nenbeg and nenend
c emax (maximal emission energy) and emin (minimal emission energy)
c and Eout (emission energy)
c
c The bin that decays is also described by a distribution of kinetic
c energies and angles in the LAB frame given by the combination of
c the elements of the array ddxrec(Zcomp,Ncomp,nex,ierec,iarec) which
c gives the spectrum as function of the recoil energies ierec and lab
c angle iarec and the elements of the array
c ddxrectot(Zcomp,Ncomp,nex) which gives the integrated cross
c section distributed over the energy-angular bin. In other words,
c ddxrectot(Zcomp,Ncomp,nex) is the sum of the elements of
c ddxrectot(Zcomp,Ncomp,nex,ierec,iarec) times the attached lab
c surface element areareclab(Zcomp,Ncomp,ierec,iarec).
c Therefore, for each decaying bin, the fraction of the decay which
c corresponds to a nucleus moving with kinetic energy
c E=Erec(Zcomp,Ncomp,ierec) and the angular direction
c in the LAB is the ratio of ddxrec(Zcomp,Ncomp,nex,ierec,iarec)*
c areareclab(Zcomp,Ncomp,ierec,iarec) with
c ddxrectot(Zcomp,Ncomp,nex,ierec,iarec)
c
c Loop over the compound nucleus kinetic energies
c
c vcm      : mean center of mass velocity for a bin
c vcmmin   : minimum center of mass velocity for a bin
c vcm      : maximum center of mass velocity for a bin
c
      do 10 nen=nenbeg,nenend
        Eout=egrid(nen)
        do 20 iang=0,nanglecont
          ang=0.5*(angcontmin(iang)+angcontmax(iang))*deg2rad
          xsgridad(nen,iang)=compspect(nen)/fourpi+
     +      preeqspect(nen)*kalbach(type,Einc,Eout,ang)
   20   continue
   10 continue
      Zix=Zindex(Zcomp,Ncomp,type)
      Nix=Nindex(Zcomp,Ncomp,type)
      ejectmass=parmass(type)*amu
      recoilmass=nucmass(Zix,Nix)*amu
      if (type.ne.0) then
        ratio=recoilmass/(ejectmass*(ejectmass+recoilmass))
      endif
      compmass=nucmass(Zcomp,Ncomp)*amu
      sumddxr=0.
c
c loop over recoil angles and velocities to calculate the fraction of
c the bin population with a given angle-velocity direction
c
      angcm=0.
      vcmav=0.
      do 110 ierec=0,maxenrec
        vcmloc(ierec)=sqrt(2.*Erec(Zcomp,Ncomp,ierec)/compmass)
        fracCMloc(ierec)=0.
        do 120 iarec=0,nanglerec
c
c Calculate fraction of the decay having this CM velocity
c characteristics (i.e. fracCMloc)
c
          if (ddxrectot(Zcomp,Ncomp,nex).gt.0.)
     +      fracCMloc(ierec)=fracCMloc(ierec)+
     +      ddxrec(Zcomp,Ncomp,nex,ierec,iarec)*
     +      areareclab(Zcomp,Ncomp,ierec,iarec)/
     +      ddxrectot(Zcomp,Ncomp,nex)
  120   continue
        if (flagrecoilav) vcmav=vcmav+fracCMloc(ierec)*vcmloc(ierec)
  110 continue
      if (flagrecoilav.and.type.gt.0) then
        ierecbeg=1
        ierecend=1
        vcmloc(1)=vcmav
        fracCMloc(1)=1.
      else
        ierecbeg=0
        ierecend=maxenrec
      endif
c
c Determine CM outgoing energies corresponding to the decay
c from nex to nexout
c
c Eejcm1     : Lower limit of CM ejectile energy bin
c Eejcm2     : Upper limit of CM ejectile energy bin
c dEejcm     : Width of CM ejectile energy bin
c pejcm1     : Impulsion corresponding to Eejcm1
c pejcm2     : Impulsion corresponding to Eejcm2
c
        do 210 ierec=ierecbeg,ierecend
          vcm=vcmloc(ierec)
          fracCM=fracCMloc(ierec)
          do 220 nen=nenbeg,nenend
            if (compspect(nen)+preeqspect(nen).eq.0.) go to 220
            Eejcm1=Ebottom(nen)
            if (Eejcm1.eq.0.) Eejcm1=0.5*egrid(1)
            Eejcm2=Etop(nen)
            dEejcm=Eejcm2-Eejcm1
            if (type.ne.0) then
               pejcm1=ejectmass*sqrt(2*ratio*Eejcm1)
               pejcm2=ejectmass*sqrt(2*ratio*Eejcm2)
               vejcm1=pejcm1/ejectmass
               vejcm2=pejcm2/ejectmass
             else
               pejcm1=Eejcm1
               pejcm2=Eejcm2
            endif
c
c Deduce recoil excitation energies and velocities corresponding
c to the ejectile energies (not perfect but should be good enough)
c
c vreccm1    : Recoil velocity corresponding to Eejcm1
c vreccm2    : Recoil velocity corresponding to Eejcm2
c
            vreccm1=pejcm1/recoilmass
            vreccm2=pejcm2/recoilmass
c
c loop over ejectile angles in the CM frame
c
c angejcm1   : Lower value of ejectile CM angle (in degrees)
c angejcm2   : Upper value of ejectile CM angle (in degrees)
c angejcmr   : Middle value of ejectile CM angle (in radians)
c
            do 230 iang=0,nanglecont
              cosejcm1=cosangcontmin(iang)
              cosejcm2=cosangcontmax(iang)
              sinejcm1=sinangcontmin(iang)
              sinejcm2=sinangcontmax(iang)
c
c Total flux that must be spread in the LAB frame
c
c fluxCM     : Total flux that must be distributed in the LAB
c
              fluxCM=xsgridad(nen,iang)*dEejcm*dcosangcont(iang)*fracCM
              if (fluxCM.eq.0.) go to 230
c
c EJECTILE and RECOIL TREATMENT
c
c
c The ejectile as well as recoil LAB angles and energies corresponding
c to the CM points are deduced
c
c flaglabddx : flag for calculation of DDX in LAB system
c Eejlab11   : LAB ejectile energy corresp. to (Eejcm1,angejcm1)
c cosejlab11 : LAB ejectile angle cosine corresp. to (Eejcm1,angejcm1)
c Eejlab12   : LAB ejectile energy corresp. to (Eejcm1,angejcm2)
c cosejlab12 : LAB ejectile angle cosine corresp. to  (Eejcm1,angejcm2)
c Eejlab21   : LAB ejectile energy corresp. to (Eejcm2,angejcm1)
c cosejlab21 : LAB ejectile angle cosine corresp. to (Eejcm2,angejcm1)
c Eejlab22   : LAB ejectile energy corresp. to (Eejcm2,angejcm2)
c cosejlab22 : LAB ejectile angle cosine corresp. to (Eejcm2,angejcm2)
c Ereclab11  : LAB recoil energy corresp. to (Eejcm1,angejcm1)
c cosreclab11: LAB recoil angle cosine corresp. to (Eejcm1,angejcm1)
c Ereclab12  : LAB recoil energy corresp. to (Eejcm1,angejcm2)
c cosreclab12: LAB recoil angle cosine corresp. to (Eejcm1,angejcm2)
c Ereclab21  : LAB recoil energy corresp. to (Eejcm2,angejcm1)
c cosreclab21: LAB recoil angle cosine corresp. to (Eejcm2,angejcm1)
c Ereclab22  : LAB recoil energy corresp. to (Eejcm2,angejcm2)
c cosreclab22: LAB recoil angle cosine corresp. to (Eejcm2,angejcm2)
c
              call cm2lab(type)
c
c calculate lab bins occupation for the ejectile in the LAB
c We assume the image in the LAB of a CM triangle is a triangle
c This is theoretically wrong but should be a good approximation
c
c labsurf1   : total surface covered in the LAB by the ejectile image
c              in the LAB of the first CM triangle (Eejcm1,angejcm1),
c              (Eejcm1,angejcm2),(Eejcm2,angejcm1)
c xlim1      : limits in the LAB ejectile energy grid
c ylim1      : limits in the LAB ejectile angular grid
c scovej1    : surface covered for each LAB energy-angle bin
c labsurf2   : total surface covered in the LAB by the ejectile image
c              in the LAB of the second CM triangle (Eejcm2,angejcm1),
c              (Eejcm2,angejcm2),(Eejcm1,angejcm2)
c xlim2      : limits in the LAB ejectile energy grid
c ylim2      : limits in the LAB ejectile angular grid
c scovej2    : surface covered for each LAB energy-angle bin
c labsurf    : Total covered surface in the LAB
c
c ddxLAB     : Eejectile double differential cross section for the LAB
c              energy-angular bin corresponding to binddx
c
              if (flaglabddx) then
               call labsurface(Zcomp,Ncomp,type,1,
     +           Eejlab11,cosejlab11,sinejlab11,
     +           Eejlab12,cosejlab12,sinejlab12,
     +           Eejlab21,cosejlab21,sinejlab21,
     +           labs1,xlim1,ylim1,scovej1,scovrec1)
               call labsurface(Zcomp,Ncomp,type,1,
     +           Eejlab21,cosejlab21,sinejlab21,
     +           Eejlab22,cosejlab22,sinejlab22,
     +           Eejlab12,cosejlab12,sinejlab12,
     +           labs2,xlim2,ylim2,scovej2,scovrec2)
               labsurf=labs1+labs2
               if (labsurf.eq.0.) goto 230
               ddxLAB=fluxCM/labsurf
c
c store ejectile flux in lab ddx array
c
c ixl,iyl    : help variables
c ddxejlab   : LAB ddx array
c areaejlab  : Total surface of LAB ddx bins
c
c Deduce the real ejectile angles by coupling with
c compound nucleus moving directions in the LAB
c
                iymax=2*nanglecont+1
                iymaxp1=iymax+1
                do ix=xlim1(1),xlim1(2)
                  do iy=ylim1(1),ylim1(2)
                    iymod=mod(iy,iymaxp1)
                    surfbin=areaejlab(type,ix,iymod)
                    if (surfbin.ne.0.) then
                      fluxadd=ddxLAB*scovej1(ix,iymod)
                      ddxejadd=fluxadd/surfbin
                      if (iymod.le.nanglecont) then
                        iys=iymod
                      else
                        iys=iymax-iymod
                      endif
                      ddxejlab(type,ix,iys)=ddxejlab(type,ix,iys)+
     +                  ddxejadd
                    endif
                  enddo
                enddo
                do ix=xlim2(1),xlim2(2)
                  do iy=ylim2(1),ylim2(2)
                    iymod=mod(iy,iymaxp1)
                    surfbin=areaejlab(type,ix,iymod)
                    if (surfbin.ne.0.) then
                      fluxadd=ddxLAB*scovej2(ix,iymod)
                      ddxejadd=fluxadd/surfbin
                      if (iymod.le.nanglecont) then
                        iys=iymod
                      else
                        iys=iymax-iymod
                      endif
                      ddxejlab(type,ix,iys)=ddxejlab(type,ix,iys)+
     +                  ddxejadd
                    endif
                  enddo
                enddo
              endif
c
c calculate lab bins occupation for the recoil in the LAB
c We assume the image in the LAB of a CM triangle is a triangle
c This is theoretically wrong but should be a good approximation
c
c labsurf1   : total surface covered in the LAB by the recoil image
c              in the LAB of the first CM triangle (Eejcm1,angejcm1),
c              (Eejcm1,angejcm2),(Eejcm2,angejcm1)
c xlim1      : limits in the LAB ejectile energy grid
c ylim1      : limits in the LAB ejectile angular grid
c scovrec1   : surface covered for each LAB energy-angle bin
c labsurf2   : total surface covered in the LAB by the image in
c              the LAB of the second CM triangle (Eejcm2,angejcm1),
c              (Eejcm2,angejcm2),(Eejcm1,angejcm2)
c xlim2      : limits in the LAB ejectile energy grid
c ylim2      : limits in the LAB ejectile angular grid
c scovrec2   : surface covered for each LAB energy-angle bin
c labsurf    : Total covered surface in the LAB
c
c ddxLAB     : Double differential cross section for the LAB
c              energy-angular bin corresponding to binddx
c
              call labsurface(Zcomp,Ncomp,type,0,
     +          Ereclab11,cosreclab11,sinreclab11,
     +          Ereclab12,cosreclab12,sinreclab12,
     +          Ereclab21,cosreclab21,sinreclab21,
     +          labs1,xlim1,ylim1,scovej1,scovrec1)
              call labsurface(Zcomp,Ncomp,type,0,
     +          Ereclab21,cosreclab21,sinreclab21,
     +          Ereclab22,cosreclab22,sinreclab22,
     +          Ereclab12,cosreclab12,sinreclab12,
     +          labs2,xlim2,ylim2,scovej2,scovrec2)
              labsurf=labs1+labs2
              if (labsurf.eq.0.) goto 230
              ddxLAB=fluxCM/labsurf
c
c store recoil flux in lab ddx array
c
c ixl,iyl    : help variables
c ddxrec     : LAB ddx array
c areareclab : Total surface of LAB ddx bins
c
c Deduce the real ejectile angles by coupling with
c compound nucleus moving directions in the LAB
c
              iymax=2*nanglerec+1
              iymaxp1=iymax+1
              do ix=xlim1(1),xlim1(2)
                do iy=ylim1(1),ylim1(2)
                  iymod=mod(iy,iymaxp1)
                  surfbin=areareclab(Zix,Nix,ix,iymod)
                  if (surfbin.ne.0.) then
                    fluxadd=ddxLAB*scovrec1(ix,iymod)
                    ddxrecadd=fluxadd/surfbin
                    sumddxr=sumddxr+fluxadd
                    if (iymod.le.nanglerec) then
                      iys=iymod
                    else
                      iys=iymax-iymod
                    endif
                    ddxrec(Zix,Nix,nexout,ix,iys)=
     +                ddxrec(Zix,Nix,nexout,ix,iys)+ddxrecadd
                  endif
                enddo
              enddo
              do ix=xlim2(1),xlim2(2)
                do iy=ylim2(1),ylim2(2)
                  iymod=mod(iy,iymaxp1)
                  surfbin=areareclab(Zix,Nix,ix,iymod)
                  if (surfbin.ne.0.) then
                    fluxadd=ddxLAB*scovrec2(ix,iymod)
                    ddxrecadd=fluxadd/surfbin
                    sumddxr=sumddxr+fluxadd
                    if (iymod.le.nanglerec) then
                      iys=iymod
                    else
                      iys=iymax-iymod
                    endif
                    ddxrec(Zix,Nix,nexout,ix,iys)=
     +                ddxrec(Zix,Nix,nexout,ix,iys)+ddxrecadd
                  endif
                enddo
              enddo
  230       continue
  220     continue
  210 continue
      ddxrectot(Zix,Nix,nexout)=ddxrectot(Zix,Nix,nexout)+sumddxr
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
