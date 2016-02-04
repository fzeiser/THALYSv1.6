      subroutine angdisrecoil
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 12, 2004
c | Task  : Recoil angular distributions for discrete states
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,type,Zix,Nix,nex,iang
      real    SS
c
c recoil local variables
c
c ejectmass  : Ejectile mass
c recoilmass : Recoil mass
c compmass   : Composite system mass
c vcm        : Compound nucleus velocity
c ddxCM      : Double differential cross section for a given CM
c              energy-angular bin
c fluxCM     : Total flux that must be distributed in the LAB
c Eejlab11   : LAB ejectile energy corresponding to (E1,ang1)
c cosejlab11 : LAB ejectile angle cosine corresponding to (E1,ang1)
c Eejlab12   : LAB ejectile energy corresponding to (E1,ang2)
c cosejlab12 : LAB ejectile angle cosine corresponding to (E1,ang2)
c Eejlab21   : LAB ejectile energy corresponding to (E2,ang1)
c cosejlab21 : LAB ejectile angle cosine corresponding to (E2,ang1)
c Eejlab22   : LAB ejectile energy corresponding to (E2,ang2)
c cosejlab22 : LAB ejectile angle cosine corresponding to (E2,ang2)
c Ereclab11  : LAB recoil energy corresponding to (E1,ang1)
c cosreclab11: LAB recoil angle cosine corresponding to (E1,ang1)
c Ereclab12  : LAB recoil energy corresponding to (E1,ang2)
c cosreclab12: LAB recoil angle cosine corresponding to (E1,ang2)
c Ereclab21  : LAB recoil energy corresponding to (E2,ang1)
c cosreclab21: LAB recoil angle cosine corresponding to (E2,ang1)
c Ereclab22  : LAB recoil energy corresponding to (E2,ang2)
c cosreclab22: LAB recoil angle cosine corresponding to (E2,ang2)
c labsurf1   : total surface covered in the LAB by the ejectile image
c              in the LAB of the first CM triangle (Eejcm1,angejcm1),
c              (Eejcm1,angejcm2),(Eejcm2,angejcm1)
c xlim1      : limits in the LAB ejectile energy grid
c ylim1      : limits in the LAB ejectile angular grid
c flaglabddx : flag for calculation of DDX in LAB system
c scovej1    : surface covered for each LAB energy-angle bin
c labsurf2   : total surface covered in the LAB by the ejectile image
c              in the LAB of the second CM triangle (Eejcm2,angejcm1),
c              (Eejcm2,angejcm2),(Eejcm1,angejcm2)
c xlim2      : limits in the LAB ejectile energy grid
c ylim2      : limits in the LAB ejectile angular grid
c scovej2    : surface covered for each LAB energy-angle bin
c labsurf    : Total covered surface in the LAB
c
c ix,iy      : help variables
c ddxejlab   : LAB ddx array
c areaejlab  : Total surface of LAB ddx bins
c scovrec1   : surface covered for each LAB energy-angle bin
c scovrec2   : surface covered for each LAB energy-angle bin
c
c ddxLAB     : Double differential cross section for the LAB
c              energy-angular bin corresponding to binddx
c sumddxrec  : help variable
c ddxrecadd  : help variable
c ddxrec     : LAB ddx array
c areareclab : Total surface of LAB ddx bins
c ddxrectot  : Total ddx integral
c
      real      compmass,ratio
      real      xsCMlev0,xsCMlev
      real      Exrec1,Exrec2
      real      pejcm1,pejcm2
      real      sumang,inorm
      real      angnorm
      real      ddxCM,fluxCM
      integer   irenorm
      real      Emaxi
      real      dEori,dEnew,renorm
      real      labsurf,labsurf1,labsurf2
      dimension scovej1(1:numen2,0:2*numangcont+1)
      dimension scovej2(1:numen2,0:2*numangcont+1)
      dimension scovrec1(0:numenrec,0:2*numangrec+1)
      dimension scovrec2(0:numenrec,0:2*numangrec+1)
      real      scovej1,scovej2,scovrec1,scovrec2
      dimension xlim1(2),ylim1(2),xlim2(2),ylim2(2)
      integer   xlim1,ylim1,xlim2,ylim2
      real      ddxLAB
      integer   iymax,iymaxp1
      integer   ix,iy,iymod,iys
      dimension ddxejlabdis(0:numpar,0:numen2,0:numangcont)
      real      ddxejlabdis
      real      ddxrecadd,sumddxrec,fluxadd,surfbin
      integer   iarec2,inex
c
c Initialisations (angcm=0 since first compound system decay)
c
      Zcomp=0
      Ncomp=0
      compmass=nucmass(Zcomp,Ncomp)*amu
      vcm=sqrt(2.*Erecinit/compmass)
      angcm=0.
c
c loop over ejectile type
c
      do 110 type=0,6
        if (parskip(type)) goto 110
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        SS=S(Zcomp,Ncomp,type)
        ejectmass=parmass(type)*amu
        recoilmass=nucmass(Zix,Nix)*amu
        if (type.ne.0) then
          ratio=recoilmass/(ejectmass*(ejectmass+recoilmass))
        endif
c
c loop over residual nucleus discrete states
c
        if (type.eq.k0) then
            xsCMlev0=xselastot
          else
            xsCMlev0=0.
        endif
        do 120 nex=0,Nlast(Zix,Nix,0)
          sumang=0.
          if (nex.ne.0) xsCMlev0=0.
          xsCMlev=xsCMlev0+xsdisc(type,nex)
          if (nex.gt.nexmax(type)) goto 120
c
c Initialise ddxejlabdis
c
          if (flaglabddx) then
            do ix=0,numen2
              do iy=0,numangcont
                ddxejlabdis(type,ix,iy)=0.
              enddo
            enddo
          endif
c
c determine the residual nucleus bin associated with nex
c
          if (nex.eq.0) then
            Exrec1=Ex(Zix,Nix,0)
            Exrec2=0.5*(Ex(Zix,Nix,1)+Ex(Zix,Nix,0))
            goto 151
          endif
          if (nex.eq.nexmax(type)) then
            Exrec2=Exmax(Zix,Nix)
            Exrec1=0.5*(Exmax(Zix,Nix)+Ex(Zix,Nix,nex-1))
            goto 151
          endif
          Exrec1=0.5*(Ex(Zix,Nix,nex)+Ex(Zix,Nix,nex-1))
          Exrec2=0.5*(Ex(Zix,Nix,nex)+Ex(Zix,Nix,nex+1))
c
c Determination of ejectile energies in the CM frame
c
  151     Eejcm1=Etotal-SS-Exrec1
          Eejcm2=Etotal-SS-Exrec2
          if (Eejcm1.lt.0..or.Eejcm2.lt.0.) goto 110
c
c Determine recoil and ejectile momentum (or velocities)
c
          if (type.ne.0) then
              pejcm1=ejectmass*sqrt(2.*ratio*Eejcm1)
              pejcm2=ejectmass*sqrt(2.*ratio*Eejcm2)
              vejcm1=pejcm1/ejectmass
              vejcm2=pejcm2/ejectmass
            else
              pejcm1=Eejcm1
              pejcm2=Eejcm2
          endif
          vreccm1=pejcm1/recoilmass
          vreccm2=pejcm2/recoilmass
c
c loop over ejectile angles in the CM frame
c
          inorm=0
 155      do 160 iang=0,nangle
            cosejcm1=cosangmin(iang)
            cosejcm2=cosangmax(iang)
            sinejcm1=sinangmin(iang)
            sinejcm2=sinangmax(iang)
c
c Total flux that must be spread in the LAB frame
c
            ddxCM=discad(type,nex,iang)
            if (inorm.eq.0) then
              sumang=sumang+ddxCM*abs(cosejcm1-cosejcm2)*twopi
            endif
c
c Renormalisation of angular distribution if inorm=1
c
            angnorm=1.
            if (inorm.eq.0) goto 160
            if (iang.eq.0) then
              if (type.ne.k0) then
                if (sumang.gt.1.e-14) then
                  angnorm=xsdisc(type,nex)/sumang
                endif
              else
                if (nex.eq.0.and.sumang.gt.1.e-14) then
                  angnorm=(xsdisc(type,nex)+xselasinc)/sumang
                else
                  if (sumang.gt.1.e-14) then
                    angnorm=xsdisc(type,nex)/sumang
                  endif
                endif
              endif
            endif
            fluxCM=ddxCM*dcosang(iang)*angnorm
            if (fluxCM.eq.0.) goto 160
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
     +            areaejlab(type,iejlab(type),iarec2)/renorm
                  ddxejlab(type,iejlab(type),iarec2)=
     +            ddxejlab(type,iejlab(type),iarec2)*renorm
                enddo
              endif
c
c calculate lab bins occupation for the ejectile in the LAB
c We assume the image in the LAB of a CM triangle is a triangle
c This is theoretically wrong but should be a good approximation
c
              call labsurface(Zcomp,Ncomp,type,1,
     +                        Eejlab11,cosejlab11,sinejlab11,
     +                        Eejlab12,cosejlab12,sinejlab12,
     +                        Eejlab21,cosejlab21,sinejlab21,
     +                        labsurf1,xlim1,ylim1,scovej1,scovrec1)
              call labsurface(Zcomp,Ncomp,type,1,
     +                        Eejlab21,cosejlab21,sinejlab21,
     +                        Eejlab22,cosejlab22,sinejlab22,
     +                        Eejlab12,cosejlab12,sinejlab12,
     +                        labsurf2,xlim2,ylim2,scovej2,scovrec2)
              labsurf=labsurf1+labsurf2
              if (labsurf.eq.0.) goto 160
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
     +              ddxejlabdis(type,ix,iys)=ddxejlabdis(type,ix,iys)+
     +                        ddxLAB/areaejlab(type,ix,iymod)*
     +                        scovej1(ix,iymod)
                enddo
              enddo
              do  ix=xlim2(1),xlim2(2)
                do iy=ylim2(1),ylim2(2)
                  iymod=mod(iy,iymaxp1)
                  if (iymod.le.nanglecont) then
                      iys=iymod
                    else
                      iys=iymax-iymod
                  endif
                  if (areaejlab(type,ix,iymod).gt.0.)
     +              ddxejlabdis(type,ix,iys)=ddxejlabdis(type,ix,iys)+
     +                        ddxLAB/areaejlab(type,ix,iymod)*
     +                        scovej2(ix,iymod)
                enddo
              enddo
            endif
c
c We check if Recoil energies are greater than the maximum
c value calculated. If it is so we accordingly modify the
c maximum recoil energy, the recoil area array as well as
c the recoil spectrum array in the lab
c
c nanglerec   : number of recoil angles
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
     +          areareclab(Zix,Nix,maxenrec,iarec2)/renorm
                do inex=0,maxex(Zix,Nix)
                  ddxrec(Zix,Nix,inex,maxenrec,iarec2)=
     +            ddxrec(Zix,Nix,inex,maxenrec,iarec2)*renorm
                enddo
              enddo
            endif
c
c calculate lab bins occupation for the recoil in the LAB
c We assume the image in the LAB of a CM triangle is a triangle
c This is theoretically wrong but should be a good approximation
c
            call labsurface(Zcomp,Ncomp,type,0,
     +                      Ereclab11,cosreclab11,sinreclab11,
     +                      Ereclab12,cosreclab12,sinreclab12,
     +                      Ereclab21,cosreclab21,sinreclab21,
     +                      labsurf1,xlim1,ylim1,scovej1,scovrec1)
            call labsurface(Zcomp,Ncomp,type,0,
     +                      Ereclab21,cosreclab21,sinreclab21,
     +                      Ereclab22,cosreclab22,sinreclab22,
     +                      Ereclab12,cosreclab12,sinreclab12,
     +                      labsurf2,xlim2,ylim2,scovej2,scovrec2)
            labsurf=labsurf1+labsurf2
            if (labsurf.eq.0.) goto 160
            ddxLAB=fluxCM/labsurf
c
c store recoil flux in lab ddx array
c
            sumddxrec=0.
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
                  ddxrec(Zix,Nix,nex,ix,iys)=ddxrec(Zix,Nix,nex,ix,iys)+
     +              ddxrecadd
                  sumddxrec=sumddxrec+fluxadd
                endif
              enddo
            enddo
            do  ix=xlim2(1),xlim2(2)
              do  iy=ylim2(1),ylim2(2)
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
                  ddxrec(Zix,Nix,nex,ix,iys)=ddxrec(Zix,Nix,nex,ix,iys)+
     +              ddxrecadd
                  sumddxrec=sumddxrec+fluxadd
                endif
              enddo
            enddo
            ddxrectot(Zix,Nix,nex)=ddxrectot(Zix,Nix,nex)+sumddxrec
  160     continue
c
c The renormalisation factor is calculated => we redo loop 160
c
          if (inorm.eq.0) then
            inorm=1
            goto 155
          endif
c
c We add the renormalized discrete component to the total ddxejlab array
c
          if (flaglabddx) then
            renorm=1.0
            if (sumang.ne.xsCMlev.and.sumang.ge.1.e-30)
     +        renorm=xsCMlev/sumang
            do ix=1,iejlab(type)
              do iy=0,nanglecont
                ddxejlab(type,ix,iy)=ddxejlab(type,ix,iy)+renorm*
     +                               ddxejlabdis(type,ix,iy)
              enddo
            enddo
          endif
  120   continue
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
