      subroutine binaryrecoil
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 12, 2004
c | Task  : Recoil for binary reaction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Zix,Nix
      integer nen
      real    SS
c
c ************************** Local variables ***************************
c
c ejectmass    : Ejectile mass
c recoilmass   : Recoil mass
c compmass     : Composite system mass
c ratio        : Help variable defining the "classical" relation between
c                the CM (kinetic+recoil) energy and the CM ejectile
c                velocity
c vcm          : CM velocity with respect to the LAB
c angcm        : CM angle with respect to the LAB
c sumejbinlab  : Integrated ejectile binary spectrum in the LAB
c sumrecbinlab : Integrated recoil binary spectrum in the LAB
c dEejcm       : CM width of ejectile emission energy bin
c Eejcm1       : CM lower limit of emission energy bin (in talys.cmb)
c Eejcm2       : CM upper limit of emission energy bin (in talys.cmb)
c pejcm1       : Impulsion corresponding to Eejcm1
c pejcm2       : Impulsion corresponding to Eejcm2
c Exrec1       : Recoil Excitation energy corresponding to Eejcm1
c Exrec        : Recoil Excitation energy corresponding to Eejcm
c Exrec2       : Recoil Excitation energy corresponding to Eejcm
c vreccm1      : Recoil velocity corresponding to Eejcm1
c vreccm2      : Recoil velocity corresponding to Eejcm2
c iex          : Loop over excitation energy bins counter
c iex1         : Recoil excitation energy index corresponding to Eejcm1
c iex2         : Recoil excitation energy index corresponding to Eejcm1
c iexmin       : minimum of iex1 and iex2
c iexmax       : maximum of iex1 and iex2
c numbinrec    : number of excitation energy bin covered
c ddxCM        : Double differential cross section for a given CM
c                energy-angular bin
c fluxCM       : Total flux that must be distributed in the LAB frame
c irenorm      : index indicating if lab array renormalisation is needed
c Emaxi        : help variable
c dEori,dEnew  : help variables
c flaglabddx   : flag for calculation of DDX in LAB system
c renorm       : renormalisation factor
c Eejlab11     : LAB energy deduced from cm2lab subroutine
c cosejlab11   : LAB angle cosine deduced from cm2lab subroutine
c Eejlab12     : LAB energy deduced from cm2lab subroutine
c cosejlab12   : LAB angle cosine deduced from cm2lab subroutine
c Eejlab21     : LAB energy deduced from cm2lab subroutine
c cosejlab21   : LAB angle cosine deduced from cm2lab subroutine
c Eejlab22     : LAB energy deduced from cm2lab subroutine
c cosejlab22   : LAB angle cosine deduced from cm2lab subroutine
c Ereclab11    : LAB energy deduced from cm2lab subroutine
c cosreclab11  : LAB angle cosine deduced from cm2lab subroutine
c Ereclab12    : LAB energy deduced from cm2lab subroutine
c cosreclab12  : LAB angle cosine deduced from cm2lab subroutine
c Ereclab21    : LAB energy deduced from cm2lab subroutine
c cosreclab21  : LAB angle cosine deduced from cm2lab subroutine
c Ereclab22    : LAB energy deduced from cm2lab subroutine
c cosreclab22  : LAB angle cosine deduced from cm2lab subroutine
c labsurf1     : total surface covered in the LAB by the ejectile image
c                in the LAB of the first CM triangle (E1,ang1),
c                (E1,ang2),(E2,ang1)
c xlim1        : limits in the LAB ejectile energy grid
c ylim1        : limits in the LAB ejectile angular grid
c scovej1      : surface covered for each ejectile LAB energy-angle bin
c scovrec1     : surface covered for each recoil  LAB energy-angle bin
c labsurf2     : total surface covered in the LAB by the ejectile image
c                in the LAB of the second CM triangle (E2,ang1),
c                (E2,ang2),(E1,ang2)
c xlim2        : limits in the LAB ejectile energy grid
c ylim2        : limits in the LAB ejectile angular grid
c scovej2      : surface covered for each ejectile LAB energy-angle bin
c scovrec2     : surface covered for each recoil  LAB energy-angle bin
c labsurf      : Total covered surface in the LAB
c ddxLAB       : Ejectile double differential cross section for the LAB
c                energy-angular bin (fluxCM/labsurf)
c ix,iy        : help variables
c ddxejlab     : LAB ddx array
c areaejlab    : Total surface of LAB ddx bins
c ddxrecadd    : help variable
c ddxrec       : LAB ddx array
c areareclab   : Total surface of LAB ddx bins
c ddxrectot    : Total ddx integral
c

      real      compmass,ratio
      real      sumejbinlab(0:6),sumrecbinlab(0:6)
      real      dEejcm
      real      pejcm1,pejcm2
      real      Exrec1,Exrec2
      integer   iex
      integer   iex1,iex2
      integer   iexmin,iexmax,numbinrec
      integer   iang
      real      ddxCM
      real      fluxCM
      integer   irenorm
      real      Emaxi
      real      dEori,dEnew,renorm
      real      labsurf1,labsurf2,labsurf
      real      scovej1(numen2,0:2*numangcont+1)
      real      scovej2(numen2,0:2*numangcont+1)
      real      scovrec1(0:numenrec,0:2*numangrec+1)
      real      scovrec2(0:numenrec,0:2*numangrec+1)
      integer   xlim1(2),ylim1(2),xlim2(2),ylim2(2)
      real      ddxLAB
      integer   iymax,iymaxp1
      integer   ix,iy,iymod,iys
      real      ddxrecadd,surfbin,fluxadd
      real      sum,sumenrec,sumenej
      integer   iarec2
c
c parskip  : logical to skip outgoing particle
c flagspec : flag for output of spectra
c ebegin   : first energy point of energy grid
c eend     : last energy point of energy grid
c
      if (flagspec) then
        do 10 type=0,6
          if (parskip(type)) goto 10
c
c *** Calculate recoil and light particle spectra in the LAB frame *****
c
c Initialisations (angcm=0 since first compound system decay)
c
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          SS=S(0,0,type)
          ejectmass=parmass(type)*amu
          recoilmass=nucmass(Zix,Nix)*amu
          compmass=nucmass(0,0)*amu
          if (type.ne.0) then
            ratio=recoilmass/(ejectmass*(ejectmass+recoilmass))
          endif
          vcm=sqrt(2.*Erecinit/compmass)
          angcm=0.
          sumejbinlab(type)=0.
          sumrecbinlab(type)=0.
c
c Loop over ejectile CM outgoing energy bins (i.e. egrid(nen))
c corresponding to the light particle spectrum component in the CM frame
c as calculated in subroutine binaryspectra.f which have been stored in
c xsbinemis(type,nen) and xsbinemisad(type,nen,iang)
c
          do 20 nen=ebegin(type),eend(type)
            if (xsbinemis(type,nen).eq.0.) goto 20
            Eejcm1=Ebottom(nen)
            if (Eejcm1.eq.0.) Eejcm1=0.5*egrid(1)
            if (Eejcm1.gt.Etotal-S(0,0,type)) goto 20
            Eejcm2=Etop(nen)
            dEejcm=Eejcm2-Eejcm1
            if (type.ne.0) then
              pejcm1=ejectmass*sqrt(2.*ratio*Eejcm1)
              pejcm2=ejectmass*sqrt(2.*ratio*Eejcm2)
              vejcm1=pejcm1/ejectmass
              vejcm2=pejcm2/ejectmass
            else
              pejcm1=Eejcm1
              pejcm2=Eejcm2
            endif
c
c Deduce recoil excitation energies and velocities and find the
c excitation energy bin in which Exrec1 and Exrec2 are located
c
            Exrec1=Exinc-SS-Eejcm1
            Exrec2=Exinc-SS-Eejcm2
            vreccm1=pejcm1/recoilmass
            vreccm2=pejcm2/recoilmass
            iex1=0
            iex2=0
            do iex=0,maxex(Zix,Nix)
              if (Exrec1.ge.Ex(Zix,Nix,iex)) iex1=iex
              if (Exrec2.ge.Ex(Zix,Nix,iex)) iex2=iex
            enddo
            iexmin=min(iex1,iex2)
            iexmax=max(iex1,iex2)
            numbinrec=iexmax-iexmin+1
c
c a specific treatment should be done we cover more than one bin
c right now we do not go into the details
c
c
c loop over ejectile angles in the CM frame
c
            do 30 iang=0,nanglecont
              cosejcm1=cosangcontmin(iang)
              cosejcm2=cosangcontmax(iang)
              sinejcm1=sinangcontmin(iang)
              sinejcm2=sinangcontmax(iang)
c
c Total flux that must be spread in the LAB frame
c
              ddxCM=xsbinemisad(type,nen,iang)
              fluxCM=ddxCM*dEejcm*dcosangcont(iang)
              if (fluxCM.eq.0.) goto 30
c
c EJECTILE and RECOIL TREATMENT
c
c Ejectile and recoil LAB angles and energies deduced from CM points
c
              call cm2lab(type)
c
c We check if ejectile energies are greater than the maximum value
c If it is so we accordingly modify the maximum ejectile energy and
c renormalise the ejectile area array as well as the ejectile spectrum
c array in the lab
c
              if (flaglabddx) then
                irenorm=0
                Emaxi=Eejlabmax(type,iejlab(type))
                dEori=Emaxi-Eejlabmin(type,iejlab(type))
                if (Eejlab11.gt.Emaxi) then
                  Emaxi=Eejlab11
                  irenorm=1
                endif
                if (Eejlab12.gt.Emaxi) then
                  Emaxi=Eejlab12
                  irenorm=1
                endif
                if (Eejlab21.gt.Emaxi) then
                  Emaxi=Eejlab21
                  irenorm=1
                endif
                if (Eejlab22.gt.Emaxi) then
                  Emaxi=Eejlab22
                  irenorm=1
                endif
                if (irenorm.eq.1) then
                  Eejlabmax(type,iejlab(type))=Emaxi
                  dEnew=Emaxi-Eejlabmin(type,iejlab(type))
                  dEejlab(type,iejlab(type))=dEnew
                  renorm=dEori/dEnew
                  do iarec2=0,nanglecont
                    areaejlab(type,iejlab(type),iarec2)=
     +              areaejlab(type,iejlab(type),iarec2)/renorm
                    ddxejlab(type,iejlab(type),iarec2)=
     +              ddxejlab(type,iejlab(type),iarec2)*renorm
                  enddo
                endif
c
c calculate lab bins occupation for the ejectile in the LAB
c We assume the image in the LAB of a CM triangle is a triangle
c This is theoretically wrong but should be a good approximation
c
                call labsurface(0,0,type,1,
     +                          Eejlab11,cosejlab11,sinejlab11,
     +                          Eejlab12,cosejlab12,sinejlab12,
     +                          Eejlab21,cosejlab21,sinejlab21,
     +                          labsurf1,xlim1,ylim1,scovej1,scovrec1)
                call labsurface(0,0,type,1,
     +                          Eejlab21,cosejlab21,sinejlab21,
     +                          Eejlab22,cosejlab22,sinejlab22,
     +                          Eejlab12,cosejlab12,sinejlab12,
     +                          labsurf2,xlim2,ylim2,scovej2,scovrec2)
                labsurf=labsurf1+labsurf2
                if (labsurf.eq.0.) goto 30
                ddxLAB=fluxCM/labsurf
c
c store ejectile flux in lab ddx array
c
                iymax=2*nanglecont+1
                iymaxp1=iymax+1
                do ix=xlim1(1),xlim1(2)
                  do iy=ylim1(1),ylim1(2)
                    iymod=mod(iy,iymaxp1)
                    if (iymod.le.nanglecont) then
                        iys=iymod
                      else
                        iys=iymax-iymod
                    endif
                    if (areaejlab(type,ix,iymod).gt.0.)
     +                ddxejlab(type,ix,iys)=ddxejlab(type,ix,iys)+
     +                ddxLAB/areaejlab(type,ix,iymod)*scovej1(ix,iymod)
                  enddo
                enddo
                do ix=xlim2(1),xlim2(2)
                  do iy=ylim2(1),ylim2(2)
                    iymod=mod(iy,iymaxp1)
                    if (iymod.le.nanglecont) then
                        iys=iymod
                      else
                        iys=iymax-iymod
                    endif
                    if (areaejlab(type,ix,iymod).gt.0.)
     +                ddxejlab(type,ix,iys)=ddxejlab(type,ix,iys)+
     +                ddxLAB/areaejlab(type,ix,iymod)*scovej2(ix,iymod)
                  enddo
                enddo
              endif
c
c We check if Recoil energies are greater than the maximum
c value calculated. If it is so we accordingly modify the
c maximum recoil energy, the recoil area array as well as
c the recoil spectrum array in the lab
c
              irenorm=0
              Emaxi=Erecmax(Zix,Nix,maxenrec)
              dEori=Emaxi-Erecmin(Zix,Nix,maxenrec)
              if (Ereclab11.gt.Emaxi) then
                Emaxi=Ereclab11
                irenorm=1
              endif
              if (Ereclab12.gt.Emaxi) then
                Emaxi=Ereclab12
                irenorm=1
              endif
              if (Ereclab21.gt.Emaxi) then
                Emaxi=Ereclab21
                irenorm=1
              endif
              if (Ereclab22.gt.Emaxi) then
                Emaxi=Ereclab22
                irenorm=1
              endif
              if (irenorm.eq.1) then
                Erecmax(Zix,Nix,maxenrec)=Emaxi
                dEnew=Emaxi-Erecmin(Zix,Nix,maxenrec)
                renorm=dEori/dEnew
                do iarec2=0,nanglerec
                  areareclab(Zix,Nix,maxenrec,iarec2)=
     +            areareclab(Zix,Nix,maxenrec,iarec2)/renorm
                  do iex=0,maxex(Zix,Nix)
                    ddxrec(Zix,Nix,iex,maxenrec,iarec2)=
     +              ddxrec(Zix,Nix,iex,maxenrec,iarec2)*renorm
                  enddo
                enddo
              endif
c
c calculate lab bins occupation for the recoil in the LAB
c We assume the image in the LAB of a CM triangle is a triangle
c This is theoretically wrong but should be a good approximation.
c
              call labsurface(0,0,type,0,
     +                        Ereclab11,cosreclab11,sinreclab11,
     +                        Ereclab12,cosreclab12,sinreclab12,
     +                        Ereclab21,cosreclab21,sinreclab21,
     +                        labsurf1,xlim1,ylim1,scovej1,scovrec1)
              call labsurface(0,0,type,0,
     +                        Ereclab21,cosreclab21,sinreclab21,
     +                        Ereclab22,cosreclab22,sinreclab22,
     +                        Ereclab12,cosreclab12,sinreclab12,
     +                        labsurf2,xlim2,ylim2,scovej2,scovrec2)
              labsurf=labsurf1+labsurf2
              if (labsurf.eq.0.) goto 30
              ddxLAB=fluxCM/labsurf/numbinrec
c
c store recoil flux in lab ddx array
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
                    if (iymod.le.nanglerec) then
                      iys=iymod
                    else
                      iys=iymax-iymod
                    endif
                    do iex=iexmin,iexmax
                      ddxrec(Zix,Nix,iex,ix,iys)=
     +                ddxrec(Zix,Nix,iex,ix,iys)+ddxrecadd
                      ddxrectot(Zix,Nix,iex)=ddxrectot(Zix,Nix,iex)+
     +                  fluxadd

                    enddo
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
                    if (iymod.le.nanglerec) then
                      iys=iymod
                    else
                      iys=iymax-iymod
                    endif
                    do iex=iexmin,iexmax
                      ddxrec(Zix,Nix,iex,ix,iys)=
     +                ddxrec(Zix,Nix,iex,ix,iys)+ddxrecadd
                      ddxrectot(Zix,Nix,iex)=ddxrectot(Zix,Nix,iex)+
     +                  fluxadd
                    enddo
                  endif
                enddo
              enddo
  30       continue
  20     continue
c
c Integrate the recoil and ejectile spectra in the LAB
c
         if (flaglabddx) then
           sumenej=0.
           do nen=1,iejlab(type)
             sum=0.
             do iang=0,nanglecont
               sum=sum+ddxejlab(type,nen,iang)*dcosangcont(iang)
             enddo
             if (nen.le.iejlab(type)) then
               sumenej=sumenej+sum*twopi*dEejlab(type,nen)
             endif
           enddo
           sumejbinlab(type)=sumenej
         endif
         sumenrec=0.
         do iex=0,maxex(Zix,Nix)
           sum=0.
           do ix=0,maxenrec
             do iy=0,nanglerec
               sum=sum+ddxrec(Zix,Nix,iex,ix,iy)*
     +                 areareclab(Zix,Nix,ix,iy)
             enddo
           enddo
           sumenrec=sumenrec+sum*twopi
         enddo
         if ((Zix.eq.0).and.(Nix.eq.0)) then
           sumrecbinlab(type)=sumenrec-xstotinc
         else
           sumrecbinlab(type)=sumenrec
         endif
  10   continue
      endif
c
c Renormalise the recoil and ejectile spectra
c
      do 110 type=0,6
        if (parskip(type)) goto 110
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        if (flaglabddx) then
          if (sumejbinlab(type).le.1.e-14) then
            renorm=1.
          else
            renorm=binemissum(type)/sumejbinlab(type)
          endif
          do nen=1,iejlab(type)
            do iang=0,nanglecont
              ddxejlab(type,nen,iang)=ddxejlab(type,nen,iang)*renorm
            enddo
          enddo
        endif
        if (sumrecbinlab(type).le.1.e-14) then
          renorm=1.
        else
          renorm=binemissum(type)/sumrecbinlab(type)
        endif
        do iex=0,maxex(Zix,Nix)
          if ((type.eq.0).and.(iex.eq.maxex(0,0))) then
            ddxrectot(Zix,Nix,iex)=ddxrectot(Zix,Nix,iex)
          else
            ddxrectot(Zix,Nix,iex)=ddxrectot(Zix,Nix,iex)*renorm
          endif
          do ix=0,maxenrec
            do iy=0,nanglerec
              if ((type.eq.0).and.(iex.eq.maxex(0,0)).and.
     +          (ix.eq.irecinit).and.(iy.eq.0)) then
                ddxrec(Zix,Nix,iex,ix,iy)=ddxrec(Zix,Nix,iex,ix,iy)
              else
                ddxrec(Zix,Nix,iex,ix,iy)=ddxrec(Zix,Nix,iex,ix,iy)*
     +          renorm
              endif
            enddo
          enddo
        enddo
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
