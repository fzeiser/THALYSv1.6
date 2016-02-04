      subroutine labsurface(Zcomp,Ncomp,type,ityp,
     +                      x1,cos1,sin1,
     +                      x2,cos2,sin2,
     +                      x3,cos3,sin3,
     +                      totsurf,limx,limy,sparte,spartr)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : December 22, 2011
c | Task  : Calculation of the way the surface triangle defined by its
c |         summits (x1,y1),(x2,y2) and (x3,y3) is distributed in a
c |         given bidimensional grid depending on ityp.
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c ityp    : =0 means ejectile, =1 means recoil
c totsurf : triangle surface
c limx    : x limits covered by the triangle
c limy    : y limits covered by the triangle
c sparte  : array containing the partial surface (filled if ityp=0)
c spartr  : array containing the partial surface (filled if ityp=1)
c
      include   "talys.cmb"
      integer   Zcomp,Ncomp,type,ityp,Zix,Nix
      real      x1,cos1,sin1,x2,cos2,sin2,x3,cos3,sin3
      real      totsurf
      dimension limx(2),limy(2)
      integer   limx,limy
      dimension sparte(numen2,0:2*numangcont+1)
      real      sparte
      dimension spartr(0:numenrec,0:2*numangrec+1)
      real      spartr
c
c ************************** Local variables ***************************
c
c recoilmode     : type of recoil calculation (1=most exact,2=simplest)
c yy1,yy2,yy3    : help variables
c yc1,yc2,yc3    : cosine of yy1,yy2,yy3
c sumsurf        : sum of partial surfaces
c ixmax          : maximum x-loop index
c ix             : x-loop counter
c xl,xu          : lower and upper x-bin coordinates
c ix1,ix2,ix3    : index locating the summits in the x-axis grid
c iymax          : maximum y-loop index
c iy             : y-loop counter
c yl,yu          : lower and upper y-bin coordinates
c iy1,iy2,iy3    : index locating the summits in the y-axis grid
c numbinscovered : total number of covered bins
c oneovernumbins : help variable to improve speed
c surfloc        : help variable to improve speed
c sumsurf        : sum of the partial surfaces
c epsx           : x-accuracy (used to distinguish 2 points)
c epsy           : y-accuracy (used to distinguish 2 points)
c
      logical   belongs
      integer   recoilmode
      real      angle1,angle2,angle3
      real      yy1,yy2,yy3,yc1,yc2,yc3
      integer   is1,is2,is3
      integer   nt1,nt2,nt3
      real      add1,add2,add3
      integer   ixmax,ix
      real      xl,xu,yl,yu
      real      xl1,xu1
      real      xl2,xu2
      real      xl3,xu3
      integer   ix1,ix2,ix3
      integer   iymax,iy
      real      yl1,yu1
      real      yl2,yu2
      real      yl3,yu3
      integer   iy1,iy2,iy3
      integer   ixmin,iymin,iysup
      integer   numbinscovered
      real      oneovernumbins
      integer   iystore
      real      surfloc
      real      sumsurf
      real      epsx
      integer   nturn
      real      epsy
      real      surfbin
      real      renorm
c
c ************************** Initializations ***************************
c
c basic initialisations
c
      recoilmode=1
      Zix=Zindex(Zcomp,Ncomp,type)
      Nix=Nindex(Zcomp,Ncomp,type)
c
c We here calculate the lab angle resulting from the coupling
c of the angle of the considered CM with respect to the beam axis
c (angcm) with the emission angles in this CM.
c
      ix1=1
      ix2=1
      ix3=1
      iy1=0
      iy2=0
      iy3=0
      cos1=max(cos1,-1.)
      cos1=min(cos1,1.)
      angle1=acos(cos1)
      if (sin1.lt.0.) angle1=twopi-angle1
      cos2=max(cos2,-1.)
      cos2=min(cos2,1.)
      angle2=acos(cos2)
      if (sin2.lt.0.) angle2=twopi-angle2
      cos3=max(cos3,-1.)
      cos3=min(cos3,1.)
      angle3=acos(cos3)
      if (sin3.lt.0.) angle3=twopi-angle3
      yy1=angle1+angcm
      yy2=angle2+angcm
      yy3=angle3+angcm
      yc1=cos(yy1)
      yc2=cos(yy2)
      yc3=cos(yy3)
      is1=int(yy1/pi)
      is2=int(yy2/pi)
      is3=int(yy3/pi)
      nt1=int(yy1/twopi)
      nt2=int(yy2/twopi)
      nt3=int(yy3/twopi)
      add1=2.*is1
      add2=2.*is2
      add3=2.*is3
      yc1=-add1+sgn(is1)*yc1
      yc2=-add2+sgn(is2)*yc2
      yc3=-add3+sgn(is3)*yc3
c
c Reference area calculated
c
      totsurf=0.5*((x1-x3)*(yc2-yc3)-(x2-x3)*(yc1-yc3))
      totsurf=abs(totsurf)
c
c if the lab area is weak, we use the simplest approximation
c
      if (totsurf.lt.1.e-14) then
        totsurf=1.e-14
        recoilmode=2
      endif
c
c ****************** Ejectile surface calculation **********************
c
      if (ityp.ne.0) then
c
c determine x and y loop limits
c
c xloop
c
        ixmax=iejlab(type)
        if (x1.gt.Eejlabmax(type,ixmax)) then
          Eejlabmax(type,ixmax)=x1
          ix1=ixmax
        else
          do 10 ix=1,ixmax
            xl1=Eejlabmin(type,ix)
            xu1=Eejlabmax(type,ix)
            if (belongs(x1,xl1,xu1)) then
              ix1=ix
              go to 11
            endif
   10     continue
        endif
   11   if (x2.gt.Eejlabmax(type,ixmax)) then
          Eejlabmax(type,ixmax)=x2
          ix2=ixmax
        else
          do 20 ix=1,ixmax
            xl2=Eejlabmin(type,ix)
            xu2=Eejlabmax(type,ix)
            if (belongs(x2,xl2,xu2)) then
              ix2=ix
              go to 21
            endif
   20     continue
        endif
   21   if (x3.gt.Eejlabmax(type,ixmax)) then
          Eejlabmax(type,ixmax)=x3
          ix3=ixmax
        else
          do 30 ix=1,ixmax
            xl3=Eejlabmin(type,ix)
            xu3=Eejlabmax(type,ix)
            if (belongs(x3,xl3,xu3)) then
              ix3=ix
              go to 31
            endif
   30     continue
        endif
c
c yloop
c
   31   iymax=2*nanglecont+1
        do 40 iy=0,iymax
          yl1=-4.*nt1
          yu1=-4.*nt1
          if (iy.le.nanglecont) then
              yl1=yl1+cosangcontmin(iy)
              yu1=yu1+cosangcontmax(iy)
            else
              yl1=yl1-cosangcontmin(iy)-2.
              yu1=yu1-cosangcontmax(iy)-2.
          endif
          if (belongs(yc1,yl1,yu1)) then
            iy1=iy+nt1*(iymax+1)
            go to 41
          endif
   40   continue
   41   do 50 iy=0,iymax
          yl2=-4.*nt2
          yu2=-4.*nt2
          if (iy.le.nanglecont) then
              yl2=yl2+cosangcontmin(iy)
              yu2=yu2+cosangcontmax(iy)
            else
              yl2=yl2-cosangcontmin(iy)-2.
              yu2=yu2-cosangcontmax(iy)-2.
          endif
          if (belongs(yc2,yl2,yu2)) then
            iy2=iy+nt2*(iymax+1)
            go to 51
          endif
   50   continue
   51   do 60 iy=0,iymax
          yl3=-4.*nt3
          yu3=-4.*nt3
          if (iy.le.nanglecont) then
              yl3=yl3+cosangcontmin(iy)
              yu3=yu3+cosangcontmax(iy)
            else
              yl3=yl3-cosangcontmin(iy)-2.
              yu3=yu3-cosangcontmax(iy)-2.
          endif
          if (belongs(yc3,yl3,yu3)) then
            iy3=iy+nt3*(iymax+1)
            go to 61
          endif
   60   continue
   61   ixmin=min(ix1,ix2)
        ixmin=min(ixmin,ix3)
        ixmax=max(ix1,ix2)
        ixmax=max(ixmax,ix3)
        iymin=min(iy1,iy2)
        iymin=min(iymin,iy3)
        iysup=max(iy1,iy2)
        iysup=max(iysup,iy3)
        limx(1)=ixmin
        limx(2)=ixmax
        limy(1)=iymin
        limy(2)=iysup
c
c SPECIAL CASES
c
c if the three summits are inside the same bin
c
        if ((ixmin.eq.ixmax).and.(iymin.eq.iysup)) then
          iymin=mod(iymin,iymax+1)
          sparte(max(ixmin,1),iymin)=totsurf
          return
        endif
c
c simplest approximation
c
        if (recoilmode.eq.2) then
          numbinscovered=(ixmax-ixmin+1)*(iysup-iymin+1)
          oneovernumbins=1./float(numbinscovered)
          surfloc=totsurf*oneovernumbins
          do ix=ixmin,ixmax
            do iy=iymin,iysup
              iystore=mod(iy,iymax+1)
              sparte(max(ix,1),iystore)=surfloc
            enddo
          enddo
          return
        endif
c
c GENERAL CASE
c
        sumsurf=0.
        do 100 ix=ixmin,ixmax
          xl=Eejlabmin(type,ix)
          xu=Eejlabmax(type,ix)
          epsx=(xu-xl)/100000.
          do 110 iy=iymin,iysup
            iystore=mod(iy,iymax+1)
            nturn=iy/(2*nanglecont+2)
            yl=-4.*nturn
            yu=-4.*nturn
            if (iystore.le.nanglecont) then
              yl=yl+cosangcontmin(iystore)
              yu=yu+cosangcontmax(iystore)
            else
              yl=yl-cosangcontmin(iystore)-2.
              yu=yu-cosangcontmax(iystore)-2.
            endif
            epsy=abs(yu-yl)/100000.
            surfbin=areaejlab(type,ix,iystore)
            call binsurface(x1,yc1,x2,yc2,x3,yc3,xl,xu,yl,yu,pi,twopi,
     +                      epsx,epsy,surfloc,surfbin)
            sparte(max(ix,1),iystore)=surfloc
            sumsurf=sumsurf+surfloc
  110     continue
  100   continue
c
c Renormalisation of ejectiles
c
        if (totsurf.ne.sumsurf) then
          if (sumsurf.ne.0.) then
            renorm=totsurf/sumsurf
            sumsurf=0.
            do ix=ixmin,ixmax
              do iy=iymin,iysup
                iystore=mod(iy,iymax+1)
                sparte(max(ix,1),iystore)=sparte(max(ix,1),iystore)*
     +            renorm
                sumsurf=sumsurf+sparte(max(ix,1),iystore)
              enddo
            enddo
          else
            numbinscovered=(ixmax-ixmin+1)*(iysup-iymin+1)
            oneovernumbins=1./float(numbinscovered)
            surfloc=totsurf*oneovernumbins
            sumsurf=0.
            do ix=ixmin,ixmax
              do iy=iymin,iysup
                iystore=mod(iy,iymax+1)
                sparte(max(ix,1),iystore)=surfloc
                sumsurf=sumsurf+surfloc
              enddo
            enddo
          endif
        endif
        return
      else
c
c ********************* Recoil surface calculation *********************
c
c determine x and y loop limits
c
c xloop
c
        if (x1.gt.Erecmax(Zix,Nix,maxenrec)) then
          Erecmax(Zix,Nix,maxenrec)=x1
          ix1=maxenrec
        else
          do 210 ix=0,maxenrec
            xl1=Erecmin(Zix,Nix,ix)
            xu1=Erecmax(Zix,Nix,ix)
            if (belongs(x1,xl1,xu1)) then
              ix1=ix
              go to 211
            endif
  210     continue
        endif
  211   if (x2.gt.Erecmax(Zix,Nix,maxenrec)) then
          Erecmax(Zix,Nix,maxenrec)=x2
          ix2=maxenrec
        else
          do 220 ix=0,maxenrec
            xl2=Erecmin(Zix,Nix,ix)
            xu2=Erecmax(Zix,Nix,ix)
            if (belongs(x2,xl2,xu2)) then
              ix2=ix
              go to 221
            endif
  220     continue
        endif
  221   if (x3.gt.Erecmax(Zix,Nix,maxenrec)) then
          Erecmax(Zix,Nix,maxenrec)=x3
          ix3=maxenrec
        else
          do 230 ix=0,maxenrec
            xl3=Erecmin(Zix,Nix,ix)
            xu3=Erecmax(Zix,Nix,ix)
            if (belongs(x3,xl3,xu3)) then
              ix3=ix
              go to 231
            endif
  230     continue
        endif
c
c yloop
c
  231   iymax=2*nanglerec+1
        do 240 iy=0,iymax
          yl1=-4.*nt1
          yu1=-4.*nt1
          if (iy.le.nanglerec) then
              yl1=yl1+cosrecmin(iy)
              yu1=yu1+cosrecmax(iy)
            else
              yl1=yl1-cosrecmin(iy)-2.
              yu1=yu1-cosrecmax(iy)-2.
          endif
          if (belongs(yc1,yl1,yu1)) then
            iy1=iy+nt1*(iymax+1)
            go to 241
          endif
  240   continue
  241   do 250 iy=0,iymax
          yl2=-4.*nt2
          yu2=-4.*nt2
          if (iy.le.nanglerec) then
              yl2=yl2+cosrecmin(iy)
              yu2=yu2+cosrecmax(iy)
            else
              yl2=yl2-cosrecmin(iy)-2.
              yu2=yu2-cosrecmax(iy)-2.
          endif
          if (belongs(yc2,yl2,yu2)) then
            iy2=iy+nt2*(iymax+1)
            go to 251
          endif
  250   continue
  251   do 260 iy=0,iymax
          yl3=-4.*nt3
          yu3=-4.*nt3
          if (iy.le.nanglerec) then
              yl3=yl3+cosrecmin(iy)
              yu3=yu3+cosrecmax(iy)
            else
              yl3=yl3-cosrecmin(iy)-2.
              yu3=yu3-cosrecmax(iy)-2.
          endif
          if (belongs(yc3,yl3,yu3)) then
            iy3=iy+nt3*(iymax+1)
            go to 261
          endif
  260   continue
  261   ixmin=min(ix1,ix2)
        ixmin=min(ixmin,ix3)
        ixmax=max(ix1,ix2)
        ixmax=max(ixmax,ix3)
        iymin=min(iy1,iy2)
        iymin=min(iymin,iy3)
        iysup=max(iy1,iy2)
        iysup=max(iysup,iy3)
        limx(1)=ixmin
        limx(2)=ixmax
        limy(1)=iymin
        limy(2)=iysup
c
c SPECIAL CASES
c
c if the three summits are inside the same bin
c
        if ((ixmin.eq.ixmax).and.(iymin.eq.iysup)) then
          iymin=mod(iymin,iymax+1)
          spartr(ixmin,iymin)=totsurf
          return
        endif
c
c simplest approximation
c
        if (recoilmode.eq.2) then
          numbinscovered=(ixmax-ixmin+1)*(iysup-iymin+1)
          oneovernumbins=1./float(numbinscovered)
          surfloc=totsurf*oneovernumbins
          do ix=ixmin,ixmax
            do iy=iymin,iysup
              iystore=mod(iy,iymax+1)
              spartr(ix,iystore)=surfloc
            enddo
          enddo
          return
        endif
c
c GENERAL CASE
c
        sumsurf=0.
        do 300 ix=ixmin,ixmax
          xl=Erecmin(Zix,Nix,ix)
          xu=Erecmax(Zix,Nix,ix)
          epsx=(xu-xl)/100000.
          do 310 iy=iymin,iysup
            iystore=mod(iy,iymax+1)
            nturn=iy/(2*nanglerec+2)
            yl=-4.*nturn
            yu=-4.*nturn
            if (iystore.le.nanglerec) then
              yl=yl+cosrecmin(iystore)
              yu=yu+cosrecmax(iystore)
            else
              yl=yl-cosrecmin(iystore)-2.
              yu=yu-cosrecmax(iystore)-2.
            endif
            epsy=abs(yu-yl)/100000.
            surfbin=areareclab(Zix,Nix,ix,iystore)
            call binsurface(x1,yc1,x2,yc2,x3,yc3,xl,xu,yl,yu,pi,twopi,
     +                      epsx,epsy,surfloc,surfbin)
            spartr(ix,iystore)=surfloc
            sumsurf=sumsurf+surfloc
  310     continue
  300   continue
c
c Renormalisation of recoils
c
        if (totsurf.ne.sumsurf) then
          if (sumsurf.ne.0.) then
            renorm=totsurf/sumsurf
            sumsurf=0.
            do ix=ixmin,ixmax
              do iy=iymin,iysup
                iystore=mod(iy,iymax+1)
                spartr(ix,iystore)=spartr(ix,iystore)*renorm
                sumsurf=sumsurf+spartr(ix,iystore)
              enddo
            enddo
          else
            numbinscovered=(ixmax-ixmin+1)*(iysup-iymin+1)
            oneovernumbins=1./float(numbinscovered)
            surfloc=totsurf*oneovernumbins
            sumsurf=0.
            do ix=ixmin,ixmax
              do iy=iymin,iysup
                iystore=mod(iy,iymax+1)
                spartr(ix,iystore)=surfloc
                sumsurf=sumsurf+surfloc
              enddo
            enddo
          endif
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
