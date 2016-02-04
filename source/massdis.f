      subroutine massdis
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 15, 2013
c | Task  : Fission fragment yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      include "gef.cmb"
      character*90 gefpath
      real         fiseps,Exfis(1000),xsfis(1000),fisepsA,partfisxs,
     +             sumpre,sumpost,xsfpapre(nummass),xsfpapost(nummass),
     +             Fmulti,beldm1(136,203),ushell1(136,203)
      integer      iz,ia,in,i,j,gefwrite,Zcomp,Ncomp,Z,A,Zix,Nix,nexend,
     +             iskip,istep,nex,nen,lengefpath
c
c ************************** Mass yields *******************************
c
c fiseps     : limit for fission cross section per nucleus
c Rfiseps    : ratio for limit for fission cross section per nucleus
c xsfistot   : total fission cross section
c nummass    : number of masses
c xsfpapre   : pre-neutron emission cross section
c xsfpapost  : post-neutron emission corrected cross section
c yieldapre  : pre-neutron emission fission yield
c yieldapost : post-neutron emission corrected fission yield
c numelem    : number of elements
c xsfpzapre  : pre-neutron emission isotopic cross section
c xsfpzapost : post-neutron emission corrected isotopic cross section
c yieldzapre : pre-neutron emission isotopic yield
c yieldzapost: post-neutron emission corrected isotopic yield
c fymodel    : fission yield model, 1: Brosa 2: GEF
c path       : directory containing structure files to be read
c gefpath    : path for GEF files
c lengefpath : character length of path for GEF files
c flagoutfy  : flag for output detailed fission yield calculation
c gefwrite   : integer for output detailed fission yield calculation
c
c Initialization
c
      fiseps=Rfiseps*xsfistot
      do 10 ia=1,nummass
        xsfpapre(ia)=0.
        xsfpapost(ia)=0.
        yieldapre(ia)=0.
        yieldapost(ia)=0.
        nupre(ia)=0.
        nupost(ia)=0.
        do 10 iz=1,numelem
          xsfpzapre(iz,ia)=0.
          xsfpzapost(iz,ia)=0.
          yieldzapre(iz,ia)=0.
          yieldzapost(iz,ia)=0.
   10 continue
      do 15 in=1,nummass-numelem
        yieldnpre(in)=0.
        yieldnpost(in)=0.
   15 continue
      do 17 i=1,numnu
        Pdisnu(i)=0.
   17 continue
      nubar=0.
      if (fymodel.eq.2) then
        gefpath=path(1:lenpath)//'fission/gef/'
        lengefpath=lenpath+12
        if (flagoutfy) then
          gefwrite=1
        else
          gefwrite=0
        endif
c
c Read nuclear structure information for GEF
c
        open (unit=4,status='unknown',
     +    file=gefpath(1:lengefpath)//'beldm.dat')
        read(4,*) beldm1
        close(4)
        do i=1,203
          do j=1,136
           beldm(i,j)=beldm1(j,i)
          end do
        end do
        open (unit=4,status='unknown',
     +    file=gefpath(1:lengefpath)//'ushell.dat')
        read(4,*) ushell1
        close(4)
        do i=1,203
          do j=1,136
           ushel(i,j)=ushell1(j,i)
          end do
        end do
        open (unit=4,status='unknown',
     +    file=gefpath(1:lengefpath)//'nucprop.dat')
        do i=1,3885
          read(4,*) (RNucTab(i,j),j=1,8)
        end do
        close(4)
      endif
c
c Loop over nuclides
c
c Zcomp      : charge number index for compound nucleus
c maxZ       : maximal number of protons away from the initial
c              compound nucleus
c Ncomp      : neutron number index for compound nucleus
c maxN       : maximal number of neutrons away from the initial
c              compound nucleus
c ZZ,Z       : charge number of residual nucleus
c AA,A       : mass number of residual nucleus
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c maxex      : maximum excitation energy bin for compound nucleus
c fiseps     : limit for fission cross section per excitation energy bin
c iskip,istep: help variables
c xsbinary   : cross section from initial compound to residual nucleus
c Ex         : excitation energy
c partfisxs  : partial fission cross section
c fisfeedex  : fission contribution from excitation energy bin
c brosafy    : subroutine for fission fragment yields based on Brosa
c              model
c disa       : normalised fission fragment mass yield per excitation
c              energy bin
c disacor    : normalised fission product mass yield per excitation
c              energy bin
c disaz      : normalised fission fragment isotope yield
c              per excitation energy bin
c disazcor   : normalised fission product isotope yield
c              per excitation energy bin
c gefran     : number of random events for GEF calculation
c
      do 20 Zcomp=0,maxZ
        do 20 Ncomp=0,maxN
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          Zix=Zindex(Zcomp,Ncomp,0)
          Nix=Nindex(Zcomp,Ncomp,0)
          if (xsfeed(Zcomp,Ncomp,-1).lt.fiseps) goto 20
          if (Zcomp.eq.0.and.Ncomp.eq.0) then
            nexend=maxex(Zcomp,Ncomp)+1
          else
            nexend=maxex(Zcomp,Ncomp)
          endif
          fisepsA=fiseps/max(3*maxex(Zcomp,Ncomp),1)
          iskip=0
          istep=4
          if (fymodel.eq.2) then
            do 30 nex=1,1000
              Exfis(nex)=0.
              xsfis(nex)=0.
   30       continue
            nen=0
          endif
          do 40 nex=nexend,0,-1
            if (nex.eq.maxex(Zcomp,Ncomp)+1) then
              excfis=Etotal
              partfisxs=xsbinary(-1)
            else
              if (mod(iskip,istep).ne.0) then
                iskip=iskip+1
                goto 40
              endif
              if (nex-istep+1.lt.0) goto 40
              if (Ex(Zcomp,Ncomp,nex-istep+1).ge.30.) then
                partfisxs=0.
                do 50 i=0,istep-1
                  partfisxs=partfisxs+fisfeedex(Zcomp,Ncomp,nex-i)
   50           continue
                if (partfisxs.ne.0) then
                  excfis=0.
                  do 60 i=0,istep-1
                    excfis=excfis+fisfeedex(Zcomp,Ncomp,nex-i)*
     +                Ex(Zcomp,Ncomp,nex-i)
   60             continue
                  excfis=excfis/partfisxs
                endif
                iskip=1
              else
                excfis=Ex(Zcomp,Ncomp,nex)
                partfisxs=fisfeedex(Zcomp,Ncomp,nex)
              endif
            endif
            if (partfisxs.gt.fisepsA) then
c
c Brosa
c
              if (fymodel.eq.1) then
                call brosafy(Zix,Nix)
                do 70 ia=1,A
                  xsfpapre(ia)=xsfpapre(ia)+disa(ia)*partfisxs
                  xsfpapost(ia)=xsfpapost(ia)+disacor(ia)*partfisxs
                  do 70 iz=1,Z
                    xsfpzapre(iz,ia)=xsfpzapre(iz,ia)+
     +                disaz(ia,iz)*partfisxs
                    xsfpzapost(iz,ia)=xsfpzapost(iz,ia)+
     +                disazcor(ia,iz)*partfisxs
   70             continue
              else
c
c GEF
c
                nen=nen+1
                Exfis(nen)=excfis
                xsfis(nen)=partfisxs
              endif
            endif
   40     continue
          if (fymodel.eq.2.and.A.le.350) then
            call geftalys(real(Z),real(A),nen,Exfis,xsfis,gefwrite,
     +        gefran)
            do 90 ia=1,A
              xsfpapre(ia)=xsfpapre(ia)+ysum(ia)
              xsfpapost(ia)=xsfpapost(ia)+ysump(ia)
              if (ia.le.200) then
                do 100 iz=1,Z
                  xsfpzapre(iz,ia)=xsfpzapre(iz,ia)+yAZ(ia,iz)
                  xsfpzapost(iz,ia)=xsfpzapost(iz,ia)+yAZp(ia,iz)
  100           continue
              endif
   90       continue
            if (xsfistot.gt.0.) then
              Fmulti=Ncomp*xsfeed(Zcomp,Ncomp,-1)
              do 110 i=1,numnu
                if (ann_sum(i).gt.0.)
     +            Pdisnu(i)=Pdisnu(i)+(Fmulti+ann_sum(i))/xsfistot
  110         continue
              do 120 ia=1,A
                if (anpre_sum(ia).gt.0.)
     +            nupre(ia)=nupre(ia)+(Fmulti+anpre_sum(ia))/xsfistot
                if (anpost_sum(ia).gt.0.)
     +            nupost(ia)=nupost(ia)+(Fmulti+anpost_sum(ia))/xsfistot
  120         continue
              nubar=nubar+(Fmulti+anMean_sum)/xsfistot
            endif
          endif
   20 continue
c
c Normalization to fission yields (sum = 2)
c
      sumpre=0.
      sumpost=0.
      do 210 ia=1,Atarget
        sumpre=sumpre+xsfpapre(ia)
        sumpost=sumpost+xsfpapost(ia)
  210 continue
      sumpre=0.5*sumpre
      sumpost=0.5*sumpost
      do 220 iz=1,Ztarget
        do 230 ia=iz+1,Atarget
          if (xsfpzapre(iz,ia).eq.0.) goto 230
          in=ia-iz
          if (in.gt.nummass-numelem) goto 230
          if (sumpre.gt.0.) yieldzapre(iz,ia)=xsfpzapre(iz,ia)/sumpre
          if (sumpost.gt.0.) yieldzapost(iz,ia)=
     +      xsfpzapost(iz,ia)/sumpost
          yieldapre(ia)=yieldapre(ia)+yieldzapre(iz,ia)
          yieldapost(ia)=yieldapost(ia)+yieldzapost(iz,ia)
          yieldnpre(in)=yieldnpre(in)+yieldzapre(iz,ia)
          yieldnpost(in)=yieldnpost(in)+yieldzapost(iz,ia)
  230   continue
  220 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
