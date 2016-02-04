      subroutine wkb(Z,A,Zix,Nix,nbar)
c
c +---------------------------------------------------------------------
c | Author: Roberto Capote, Arjan Koning, and Stephane Goriely
c | Date  : December 13, 2013
c | Task  : Initialization of WKB approximation for fission
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nsmooth
      integer Z,A,nbar,nbarr,i,j,n1,Find_Extrem,Zix,Nix,ii
      real    rmiu,centr,heigth,width,ucentr,uheigth,uwidth,uexc0,uexc1,
     +        tdir,dE,uexc,phase(2*numbar),tff(2*numbar)
      character*9 filewkb
c
c ******* Finding maxima and minima of the deformation energy curve ****
c
      nsmooth=5
      rmiu = 0.054*A**(5./3.)
c initializes iiextr() and nextr
      nextr = Find_Extrem(nsmooth)
      nextr=min(nextr,2*numbar-1)
      nbarr = nextr/2 + 1
      nbar=nbarr
c     nwell = nextr/2
c     Fitting parabola
      filewkb='         '
      write(filewkb,'("wkb",2i3.3)') Z,A
      if (flagfisout) then
        open (unit=22,status='unknown',file=filewkb)
        write(22,'("Z=",i3," A=",i3)') Z,A
      endif
      do j=1,2*numbar
        Vheight(j)=0.
        Vwidth(j)=0.
      enddo
      do j=1,nextr
        CALL ParabFit(iiextr(j),3,rmiu,betafis,vfis,
     &  centr,  heigth,  width, ucentr, uheigth, uwidth)
        if(width.LT.0.05d0) CYCLE ! Skipping very narrow peaks
        if (flagfisout) then
          write(22,*) ' Def: ',betafis(iiextr(j)),
     &      ' (',centr,' +/- ',ucentr,')'
          write(22,*) ' Heigth :',vfis(iiextr(j)),
     &      ' (',heigth,' +/- ',uheigth,')'
          write(22,*) ' Width :', width, ' +/- ', uwidth
          write(22,*) '*******************************************'
        endif
c--------------------------------------------------
c       Initializes parabola's parameters
c
c       The real height of the barrier is used (not the fitted
c       parabola height)
c
        Vheight(j) = vfis(iiextr(j))
        Vwidth(j) = width
        if (mod(j,2).eq.1) then
          ii=(j+1)/2
          fbarrier(Zix,Nix,ii)=vfis(iiextr(j))
          fwidth(Zix,Nix,ii)=width
        endif
c       The real position of the barrier is used (not the fitted
c       parabola centr)
        Vpos(j) = betafis(iiextr(j))
c--------------------------------------------------
      enddo
      Vpos(nextr+1) = 100.d0
C
c     PAUSE
C==================================================================
C     Calculating the limit for the energy loop
C     (Depends on barriers' height and depth)
C
c     ustep = 0.01d0
      uexc0 = 1.d0
      uexc1 = 7.d0
      uexc0 = Vheight(2)
      if(nextr.eq.5) uexc0 = min(uexc0,Vheight(4))
      uexc1 =     max(Vheight(1),Vheight(3))
      if(nextr.eq.5) uexc1 =  max(uexc1,Vheight(5))
C==================================================================
      if (flagfisout)
     &  write(22,'(1x,A5,2x,18(A10,1x))')  ' Uexc', '   Tdir   ',
     &  'TaTb/Ta+Tb','    Ta    ','    Tb    ', '    Tc    '
C
C     ENERGY LOOP
C
      n1=3*nbinswkb/4
      uexc=0.
      Uwkb(Zix,Nix,0)=0.
      do j=1,nbar
        Twkb(Zix,Nix,0,j)=0.
      enddo
      do i=1,nbinswkb
        if (i.le.n1) then
          dE=uexc1/n1
        else
          dE=0.2
        endif
        uexc=uexc+dE
        CALL WKBFIS(Z,A-Z,uexc, tff, phase, tdir)
        Uwkb(Zix,Nix,i)=uexc
        do j=1,nbar
          Twkb(Zix,Nix,i,j)=tff(2*j-1)
          if (flagfisout)
     +      write(22,'("i=",i2," E=",f8.3," j=",i2," T=",1p,e12.3)')
     +      i,uexc,j,Twkb(Zix,Nix,i,j)
        enddo
      enddo
      if (flagfisout) close (22)
      END

C===================================================================
      subroutine wkbfis(iz,in,uexcit, tff, phase, tdir)
C
C     Roberto Capote, IAEA/NDS
C     e-mail: r.capotenoy@iaea.org, 24 August 2006, v1.00
C             rcapotenoy@yahoo.com
C
C     Mihaela Sin, University of Bucharest
C     e-mail: mihaela.sin@gmail.com
C
C     For theory behind see
C     M. Sin, R. Capote, A. Ventura et al, Phys.Rev.C74, 014608 (2006)
C
C     WKBFIS calculates WKB momentum integrals by Gauss-Legendre method
C     for any given one-dimensional barrier and a given excitation
C     energy.
C
C     * uexcit is the excitation energy
C     * tfis is the fission transmission coefficient (see eq.(4) PRC)
C     * tff(extrema), phase(extrema) are the transmission coefficients
C       and phase integrals for all extrema
C       tff(1)= Ta, tff(3)= Tb, tff(5)= Tc (for triple humped barrier)
C
c     IMPLICIT none
      include "talys.cmb"
      INTEGER iz,in
c     include 'vdeform.inc'
      real uexcit, tdir, tff(2*numbar), phase(2*numbar)
      real tdirv(2*numbar)
C-------------------------------------------------------------------
      real Uexc, Smiu
      INTEGER K
      COMMON /VARGS/ Uexc, Smiu, K

c     real pi
c     PARAMETER (pi = 3.14159259d0)
      real dmom,abserr,rmiu,epsa,epsb,dummy,eps,deps,vdef,phasecal
      INTEGER j,ieps

C
C     FUNCTIONS
      real   Fmoment, FmomentParab, GaussLegendre41, FindIntersect
      EXTERNAL FmomentParab, Fmoment

      Uexc   = uexcit

      rmiu = 0.054d0*(iz+in)**(5.d0/3.d0)
C     Below is an square root of (MIU divided by 2)
      smiu = sqrt(0.5d0*rmiu)

      DO k=1, 2*numbar ! barriers and wells
        phase(k) = 0.d0
        tff(k)   = 0.d0
        tdirv(k) = 0.d0
      enddo
C-----------------------------------------------------------------------
C     Momentum integrals are calculated
      if (flagfisout) then
        write(22,*)
        write(22,*) ' Excitation Energy : ',uexcit
      endif

      DO k=1, nextr ! barriers and wells

          if(mod(k,2).eq.1) then

C         BARRIERS

            if(uexcit.ge.Vheight(k)) then
C
C           For excitation energies above barriers the transmission
C           coefficient is calculated by Hill-Wheeler formula
C
            if (Vwidth(k).gt.0) then
              dmom = pi*(Vheight(k) - uexcit)/Vwidth(k)
            else
              dmom=-50.
            endif
            phase(k)   = min(dmom,50.)
            tff(k) = 1.d0/(1.d0 + DEXP(2.d0*dmom))
            if (flagfisout)
     &      write(22,'(1x,A6,I2,A10,d10.3,A3,d10.3,A15)')
     &      ' BARR ',k,'  Mom.Int=',dmom,' T=',tff(k),' (Hill-Wheeler)'

            else
C
C           For excitation energies below barriers the transmission
C           coefficient is calculated by WKB method
C
            epsa = FindIntersect(uexcit,iiextr(k-1),iiextr(k),.false.)
            epsb = FindIntersect(uexcit,iiextr(k),iiextr(k+1),.false.)
C
C           Calculating phase integral for real shape
C
            dmom = GaussLegendre41(Fmoment,epsa,epsb,abserr)
            if(flagfisout.and.dmom.gt.0.d0.and.abserr.gt.dmom*0.03) then
              write(22,*) ' WARNING: For extremum ',k,
     &                   ' phase integral is not accurate (',
     &        abserr/dmom*100.d0,' %)'
            endif
            phase(k) = min(dmom,50.)
            tff(k) = 1.d0/(1.d0 + DEXP(2.d0*phase(k)))
C
            if (flagfisout) then
         write(22,'(1x,A6,I2, A10,f7.4,A4,f7.4,1x,A9,d10.3,A3,d10.3)')
     &      ' BARR ',k,' at beta2 ',epsa,' to ',epsb,
     &      ' Mom.Int=',phase(k),' T=',tff(k)
         deps=(epsb-epsa)/50.
         eps=epsa-deps
         phasecal=0.
         do ieps=1,50
           eps=eps+deps
           if (mod(ieps,5).eq.0)
     &     write(22,'(10x,"eps=",f6.3," Vdef-E=",f8.3)') eps,
     &       vdef(eps)-Uexc
           if (vdef(eps).ge.Uexc.and.vdef(eps+deps).ge.Uexc)
     &     phasecal=phasecal+2.*smiu*((vdef(eps)-Uexc)**0.5+
     &       (vdef(eps+deps)-Uexc)**0.5)/2.*deps
         enddo
         if (phasecal.lt.30.)
     &     write(22,'(10x," Kcal=",f10.4," Tcal=",1p,e12.4)')
     &     phasecal,1./(1.+exp(2.*phasecal))
         endif

            endif

          else

C         WELLS

            if(uexcit.LE.Vheight(k)) then
C
C           Excitation energies below the well
C
            phase(k) = 0.d0

            else
C
C           For excitation energies above the well the transmission
C           coefficient is calculated by WKB method
C
            epsa = FindIntersect(uexcit,iiextr(k-1),iiextr(k),.true.)
            epsb = FindIntersect(uexcit,iiextr(k),iiextr(k+1),.true.)
C
C           Calculating phase integral for real shape
C
            dmom = GaussLegendre41(Fmoment,epsa,epsb,abserr)
         if (flagfisout.and.dmom.gt.0.d0.and.abserr.gT.dmom*0.03)
     &        write(22,*) ' WARNING: For extremum ',k,
     &        ' phase integral is not accurate (',
     &        abserr/dmom*100.d0,' %)'
            phase(k) = min(dmom,50.)
            if (flagfisout)
     &    write(22,'(1x,A6,I2, A10,f7.4,A4,f7.4,1x,A9,d10.3,A3,d10.3)')
     &      ' WELL ',k,' at beta2 ',epsa,' to ',epsb,
     &      ' Mom.Int=',phase(k)

            endif

          endif

      ENDDO

C
C     Fission transmission for double/triple humped barrier
C
C     Iteration over barriers
C
      if (nextr.gt.0) tdirv(nextr) = tff(nextr)
      do k=nextr-2,1,-2
         dmom = (1.d0 - tff(k))*(1.d0 - tdirv(k+2))
         if(k.gt.1) then
           tdirv(k) = tff(k)*tdirv(k+2) /
     &     (1.d0 + 2.d0*sqrt(dmom)*cos(2.d0*phase(k+1)) + dmom)
         else
           tdirv(1) = tff(1)*tdirv(3) /
     &     (1.d0 + 2.d0*sqrt(dmom)*cos(2.d0*phase(2)) + dmom)
         endif
      enddo
        tdir = tdirv(1)
        dummy = 1.d0/(1.d0/tff(1) + 1.d0/tff(3))
        if(nextr.gt.3)
     &  dummy = 1.d0/(1.d0/tff(1) + 1.d0/tff(3) + 1.d0/tff(5))

      if (flagfisout) write(22,'(1x,f5.2,2x,21(d10.3,1x))')
     &   uexc,tdir,dummy,(tff(j),j=1,nextr),(phase(j),j=1,nextr)

      RETURN
      END

      real function Fmoment(Eps)
C
C     Integrand (To be called from Gauss-Legendre integration routine)
C
C     To be defined as external function
C
      IMPLICIT NONE
      real Eps, Vdef
      real Uexc, Smiu
      INTEGER K
      COMMON /VARGS/ Uexc, Smiu, K
      Fmoment = 2.d0*SMIU* sqrt( abs (UEXC - Vdef(Eps)) )
      return
      end

      real function FmomentParab(Eps)
C
C     Integrand (To be called from Gauss-Legendre integration routine)
C
C     To be defined as external function
C
      IMPLICIT NONE
      real Eps, VdefParab
      real Uexc, Smiu
      INTEGER K
      COMMON /VARGS/ Uexc, Smiu, K
      FmomentParab = 2.d0*Smiu* sqrt( abs (Uexc - VdefParab(Eps)) )
      return
      end

      real function VdefParab(EPS)
C
C     This function calculates parabolic shape of the deformation energy
C
C     Called by gaussian integration
C
c     IMPLICIT NONE
      include "talys.cmb"
      real EPS
c     INCLUDE 'vdeform.inc'
      real Uexc, Smiu
      INTEGER K
      COMMON /VARGS/ Uexc, Smiu, K
      VdefParab = Vheight(K) + (-1)**K*(SMIU*Vwidth(K)*(EPS-Vpos(K)))**2
      return
      end

      real function Vdef(EPS)
C
C     This function calculates real shape of the deformation energy
C     by linear interpolation to obtain the value of the barrier Vdef
C     at deformation EPS (needed to integrte this function)
C
C     Called by gaussian integration
C
c     IMPLICIT NONE
      include "talys.cmb"
      real EPS
c     INCLUDE 'vdeform.inc'
C
C     Local variables
      INTEGER idef
      real vi, vip, ei, eip
C
C     The four lines below is a very simple (and inefficient) search
C     Should be replaced by efficient search routine to locate the
C     element idef of the array betafis()
      idef=1
      do while (EPS.GT.betafis(idef) .and. idef.LE.nbeta)
        idef = idef + 1
      enddo
      if (idef.ne.1) idef = idef - 1

      vi  = vfis(idef)
      vip = vfis(idef+1)
      ei  = betafis(idef)
      eip = betafis(idef+1)

      if(ei.eq.eip) then
C       Special case treated here to avoid division by zero
C       We assume that in this case vi = vip
        Vdef = vi
        return
      endif

      Vdef = vi + (EPS-ei)/(eip-ei)*(vip-vi)
      return
      end

      real FUNCTION FindIntersect(uexc,ja,jb,iswell)
C
C     Calculates solutions (beta2(i)) of the equation
C         V(beta2) = Excitation Energy
C
C     If solution not found, then assign first or last point
C     of the interval depending on whether we are solving the equation
C     in the right(left) side of the barrier(well) case
C
      include "talys.cmb"
c     IMPLICIT NONE
      INTEGER ja, jb
          real uexc
      LOGICAL iswell
c     INCLUDE 'vdeform.inc'
C
C     Local variables
      INTEGER j, is0, is1
      real slope

      is0 = -1
      IF(ABS(uexc-vfis(ja)).EQ.uexc-vfis(ja)) is0 = 1
      DO j=ja,jb
C      checking for a sign change in Uexc-Vdef
       is1 = -1
       IF(ABS(uexc-vfis(j)).EQ.uexc-vfis(j)) is1 = 1
       IF(is1.EQ.is0) CYCLE
C      Sign of (Uexc-vfis(j)) changed, calculating the precise value
C      of the deformation EPS at which Uexc = Vdef
       FindIntersect = betafis(j-1) + (betafis(j)-betafis(j-1))*
     >       (uexc-vfis(j-1))/(vfis(j)-vfis(j-1))
       RETURN
      ENDDO
C
C     Below is the analysis if intersection not found in [ja,jb]
C
      slope = vfis(jb) - vfis(ja)
      IF(iswell) then
C       WELLS
        IF(slope . ge. 0) then
          FindIntersect = betafis(jb) ! ascending
        ELSE
          FindIntersect = betafis(ja) ! descending
        ENDIF
      ELSE
C       BARRIERS
        IF(slope . ge. 0) then
          FindIntersect = betafis(ja) ! ascending
        ELSE
          FindIntersect = betafis(jb) ! descending
        ENDIF
      ENDIF
      RETURN
      END

      INTEGER FUNCTION Find_Extrem (Nsmooth)
C
C     Find all extremums of the smoothed deformation energy curve
C
C     Nsmooth - Number of smoothing points
C
      include "talys.cmb"
c     IMPLICIT NONE
      INTEGER Nsmooth
c     INCLUDE 'vdeform.inc'
      LOGICAL logmin, logmax
      INTEGER j,k
      INTEGER iext
      iext=0


      Find_Extrem = 0
      iiextr(0)=1
CC----------------------------------------------------------------------
C     determination of the minima   and maxima
      do j=nsmooth+1 , nbeta-nsmooth

        logmax=.true.
        do k=j-nsmooth,j+nsmooth
          if (k.eq.j) cycle
          if (vfis(k).gt.vfis(j)) logmax=.false.
        enddo
        if (logmax) then
          iext=iext+1
          iiextr(iext)=j
        endif
        if (iext.eq.2*numbar-1) exit
        logmin=.true.
        do k=j-nsmooth,j+nsmooth
          if (k.eq.j.or.k.lt.1) cycle
          if (vfis(k).lt.vfis(j)) logmin=.false.
        enddo
        if (logmin) then
          iext=iext+1
          iiextr(iext)=j
        endif
        if (iext.eq.2*numbar-1) exit
      enddo
      Find_Extrem = iext
      iiextr(iext+1)= nbeta
      return
      end
      real FUNCTION GaussLegendre41(F,Ea,Eb,ABSERR)
      IMPLICIT NONE
      real F
      real Eb,Ea,ABSERR
      real wg(10),xgk(21),wgk(21)
      real CENTR1,HLGTH1,RESG1,RESK1
      INTEGER J,JTW,JTWM1
      real ABSC,FSUM,ABSCM1
      EXTERNAL F
      SAVE WG, XGK, WGK
C
C     THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C     BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C     CORRESPONDING WEIGHTS ARE GIVEN.
C
C     XG - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
C     WG - WEIGHTS OF THE 20-POINT GAUSS RULE
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      DATA WG  (  1) / 0.0176140071 3915211831 1861962351 853 D0 /
      DATA WG  (  2) / 0.0406014298 0038694133 1039952274 932 D0 /
      DATA WG  (  3) / 0.0626720483 3410906356 9506535187 042 D0 /
      DATA WG  (  4) / 0.0832767415 7670474872 4758143222 046 D0 /
      DATA WG  (  5) / 0.1019301198 1724043503 6750135480 350 D0 /
      DATA WG  (  6) / 0.1181945319 6151841731 2377377711 382 D0 /
      DATA WG  (  7) / 0.1316886384 4917662689 8494499748 163 D0 /
      DATA WG  (  8) / 0.1420961093 1838205132 9298325067 165 D0 /
      DATA WG  (  9) / 0.1491729864 7260374678 7828737001 969 D0 /
      DATA WG  ( 10) / 0.1527533871 3072585069 8084331955 098 D0 /
C
      DATA XGK (  1) / 0.9988590315 8827766383 8315576545 863 D0 /
      DATA XGK (  2) / 0.9931285991 8509492478 6122388471 320 D0 /
      DATA XGK (  3) / 0.9815078774 5025025919 3342994720 217 D0 /
      DATA XGK (  4) / 0.9639719272 7791379126 7666131197 277 D0 /
      DATA XGK (  5) / 0.9408226338 3175475351 9982722212 443 D0 /
      DATA XGK (  6) / 0.9122344282 5132590586 7752441203 298 D0 /
      DATA XGK (  7) / 0.8782768112 5228197607 7442995113 078 D0 /
      DATA XGK (  8) / 0.8391169718 2221882339 4529061701 521 D0 /
      DATA XGK (  9) / 0.7950414288 3755119835 0638833272 788 D0 /
      DATA XGK ( 10) / 0.7463319064 6015079261 4305070355 642 D0 /
      DATA XGK ( 11) / 0.6932376563 3475138480 5490711845 932 D0 /
      DATA XGK ( 12) / 0.6360536807 2651502545 2836696226 286 D0 /
      DATA XGK ( 13) / 0.5751404468 1971031534 2946036586 425 D0 /
      DATA XGK ( 14) / 0.5108670019 5082709800 4364050955 251 D0 /
      DATA XGK ( 15) / 0.4435931752 3872510319 9992213492 640 D0 /
      DATA XGK ( 16) / 0.3737060887 1541956067 2548177024 927 D0 /
      DATA XGK ( 17) / 0.3016278681 1491300432 0555356858 592 D0 /
      DATA XGK ( 18) / 0.2277858511 4164507808 0496195368 575 D0 /
      DATA XGK ( 19) / 0.1526054652 4092267550 5220241022 678 D0 /
      DATA XGK ( 20) / 0.0765265211 3349733375 4640409398 838 D0 /
      DATA XGK ( 21) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0030735837 1852053150 1218293246 031 D0 /
      DATA WGK (  2) / 0.0086002698 5564294219 8661787950 102 D0 /
      DATA WGK (  3) / 0.0146261692 5697125298 3787960308 868 D0 /
      DATA WGK (  4) / 0.0203883734 6126652359 8010231432 755 D0 /
      DATA WGK (  5) / 0.0258821336 0495115883 4505067096 153 D0 /
      DATA WGK (  6) / 0.0312873067 7703279895 8543119323 801 D0 /
      DATA WGK (  7) / 0.0366001697 5820079803 0557240707 211 D0 /
      DATA WGK (  8) / 0.0416688733 2797368626 3788305936 895 D0 /
      DATA WGK (  9) / 0.0464348218 6749767472 0231880926 108 D0 /
      DATA WGK ( 10) / 0.0509445739 2372869193 2707670050 345 D0 /
      DATA WGK ( 11) / 0.0551951053 4828599474 4832372419 777 D0 /
      DATA WGK ( 12) / 0.0591114008 8063957237 4967220648 594 D0 /
      DATA WGK ( 13) / 0.0626532375 5478116802 5870122174 255 D0 /
      DATA WGK ( 14) / 0.0658345971 3361842211 1563556969 398 D0 /
      DATA WGK ( 15) / 0.0686486729 2852161934 5623411885 368 D0 /
      DATA WGK ( 16) / 0.0710544235 5344406830 5790361723 210 D0 /
      DATA WGK ( 17) / 0.0730306903 3278666749 5189417658 913 D0 /
      DATA WGK ( 18) / 0.0745828754 0049918898 6581418362 488 D0 /
      DATA WGK ( 19) / 0.0757044976 8455667465 9542775376 617 D0 /
      DATA WGK ( 20) / 0.0763778676 7208073670 5502835038 061 D0 /
      DATA WGK ( 21) / 0.0766007119 1799965644 5049901530 102 D0 /
C Integrating from Ea to Eb, converting to symmetric grid from -1 to +1
      CENTR1 = 0.5D+00*(Ea+Eb)
      HLGTH1 = 0.5D+00*(Eb-Ea)
C
C     COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO
C     THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR USING ONLY 21 points
C
      RESG1 = 0.0D+00
      RESK1 = WGK(21)*F(CENTR1)
      DO J=1,10
        JTW = J*2
        JTWM1 = JTW-1
        ABSC = HLGTH1*XGK(JTW)
        FSUM = F(CENTR1-ABSC) + F(CENTR1+ABSC)
        ABSCM1 = HLGTH1*XGK(JTWM1)
        RESG1 = RESG1+WG(J)*FSUM
        RESK1 = RESK1+WGK(JTW)*FSUM+WGK(JTWM1)*
     &              (F(CENTR1-ABSCM1)+F(CENTR1+ABSCM1))
      ENDDO

      GaussLegendre41 = RESK1*HLGTH1
      ABSERR = ABS((RESK1-RESG1)*HLGTH1)

      RETURN
      END
      SUBROUTINE ParabFit(Imax,Npfit,RMIU,EPS,VDEFORM,
     &                     CENTR,  HEIGTH,  WIDTH,
     &                    uCENTR, uHEIGTH, uWIDTH)
      IMPLICIT NONE
      INTEGER Imax  ! Position of the maximum
      INTEGER Npfit ! Number of points to fit
      REAL   RMIU   ! MIU = 0.054d0*ANUC**(5.d0/3.d0)
      REAL   EPS(*) ! Deformation Array          (X)
      REAL   VDEFORM(*)! Deformation Energy Array (Y)
      REAL    CENTR,  HEIGTH,  WIDTH ! BARRIER'S PARAMETER
      REAL   uCENTR, uHEIGTH, uWIDTH ! BARRIER'S PARAMETER UNCERTAINTIES

C     Local variables
      integer i,npts,ierr
      real   x(1024),y(1024),a(3),siga(3),chisqr

      npts = 0
      do i = Imax - Npfit, Imax + Npfit
        npts = npts + 1
        x(npts) = EPS(i) - EPS(imax)
        y(npts) = VDEFORM(i)
      enddo
C
C FITS A POLYNOMIAL OF THE FORM:
C   Y=A(1)+A(2)*X+A(3)*X**2+...+A(NTERMS)*X**(NTERMS-1)
C THOUGH THE NPTS POINTS (X,Y)
C IERR=0 OK
C IERR=2 SINGULAR MATRIX
C IERR=3 REQUIREMENTS FOR USE OF THIS ROUTINE NOT FULLFILLED
C
      CALL WPLFT(npts,3,ierr,1,x,y,a,siga,chisqr)
C     do i = 1,npts
C       write(*,'(1x,F5.3,2x,d15.8,3x,d15.8)') x(i), y(i),
C    &              a(1)+a(2)*x(i)+a(3)*x(i)**2
C     enddo
C     write(*,*) a(1),a(2),a(3)
C     write(*,*) sigmaa(1),sigmaa(2),sigmaa(3)
C
c write(22,*) ' PARABOLIC FITTING CHISQR = ',CHISQR
       HEIGTH = a(1)
      uHEIGTH = siga(1)
       WIDTH  =  SQRT(2.d0*abs(   a(3))/RMIU)
      uWIDTH  =  SQRT(2.d0*abs(siga(3))/RMIU)
       CENTR  = EPS(imax)
      IF(a(3).ne.0.d0) CENTR = 0.5d0 * a(2) / a(3)
      uCENTR  = ABS(CENTR) * SQRT((siga(2)/a(2))**2 + (siga(3)/a(3))**2)
      CENTR = EPS(imax) + CENTR
      return
      end

      SUBROUTINE WPLFT(NPTS,NTERMS,IERR,MODE,X,Y,A,SIGMAA,CHISQR)

C FITS A POLYNOMIAL OF THE FORM:
C   Y=A(1)+A(2)*X+A(3)*X**2+...+A(NTERMS)*X**(NTERMS-1)
C THOUGH THE NTPS POINTS (X,Y)
C MODE 0 = WEIGTH 0 != UNWEIGHTED
C IERR=0 OK
C IERR=2 SINGULAR MATRIX
C IERR=3 REQUIREMENTS FOR USE OF THIS ROUTINE NOT FULLFILLED
C REQUIREMENTS ARE:
C        MAX NUM TERMS = 12
C        MIN NUM TERMS = 2
C        NTPS >= NTERMS
        INTEGER MAXDEG
        PARAMETER (MAXDEG=3)

        INTEGER NPTS, NTERMS, IERR, MODE
        real X(NPTS), Y(NPTS), W(NPTS)
        real A(NTERMS), SIGMAA(NTERMS), CHISQR
        DOUBLE PRECISION ARRAY(MAXDEG,MAXDEG)

        DOUBLE PRECISION XMEAN(MAXDEG)
        DOUBLE PRECISION SIGMAX(MAXDEG), WMEAN, FNPTS
        DOUBLE PRECISION SUM,YMEAN,DET,FREE1,SIGMA,VARNCE,CHISQ
        DOUBLE PRECISION R(MAXDEG)
        INTEGER WEIGHT
        INTEGER DEG,I,J,K

        WEIGHT=0
        IERR = 0
        DEG=NTERMS-1
        IF(DEG.LT.1.OR.DEG.GT.MAXDEG.OR.NPTS.LT.NTERMS) THEN
          IERR=3
          RETURN
        ENDIF

C       INITIALIZE SUMS AND ARRAYS
        SUM = 0.D0
        YMEAN=0.D0
        SIGMA = 0.D0
        CHISQ = 0.D0
        A(1) = 0.
        SIGMAA(1) = 0.
        DO 28 J=1,DEG
        XMEAN(J)=0.D0
        SIGMAX(J)=0.D0
        R(J)=0.D0
        A(1+J) = 0.
        SIGMAA(1+J) = 0.
C       DO 28 K=1,DEG
C28      ARRAY(J,K)=0.D0
        DO 28 K=1,NTERMS
28      ARRAY(J,K)=0.D0


C       ACCUMULATE WEIGHTED SUMS

        DO 50 I = 1, NPTS
        W(I) = ABS(Y(I))
        IF( MODE.NE.WEIGHT ) W(I) = 1.
        SUM = SUM + DBLE(W(I))
C       YMEAN = Sy
        YMEAN = YMEAN + DBLE(W(I)*Y(I))
C       XMEAN(1)=Sx, XMEAN(2)=Sxx
C       S=NPTS*(NPTS+1)/2
        DO 50 J = 1, DEG
50      XMEAN(J) = XMEAN(J) + DBLE(W(I)*X(I)**J)
        YMEAN = YMEAN / SUM
        DO 53 J = 1, DEG
53      XMEAN(J) = XMEAN(J) / SUM
        FNPTS = NPTS
        WMEAN = SUM / FNPTS
        DO 57 I = 1, NPTS
57      W(I) = W(I) / WMEAN

        DO 67 I=1,NPTS
        SIGMA=SIGMA+DBLE(W(I)*(Y(I)-YMEAN)**2)
        DO 67 J=1,DEG
        SIGMAX(J)=SIGMAX(J)+DBLE(W(I)*(X(I)**J-XMEAN(J))**2)
        R(J)=R(J)+DBLE(W(I)*(X(I)**J-XMEAN(J))*(Y(I)-YMEAN))
        DO 67 K=1,J
67      ARRAY(J,K)=ARRAY(J,K)+
     1             DBLE(W(I)*(X(I)**J-XMEAN(J))*(X(I)**K-XMEAN(K)))
        FREE1 = NPTS-1
        SIGMA=sqrt(SIGMA/FREE1)
        DO 78 J=1,DEG
        SIGMAX(J)=sqrt(SIGMAX(J)/FREE1)
        R(J)=R(J)/(FREE1*SIGMAX(J)*SIGMA)
        DO 78 K=1,J
        ARRAY(J,K)=ARRAY(J,K)/(FREE1*SIGMAX(J)*SIGMAX(K))
78      ARRAY(K,J)=ARRAY(J,K)

C       INVERT SYMMETRIC MATRIX

        CALL MATINVM( DEG, DET, ARRAY )
        IF(DET.NE.0.D0) GO TO 101
        IERR = 2
        RETURN

C       CALCULATE COEFFICIENTS

101     A(1)=YMEAN
        DO J=1,DEG
          DO K=1,DEG
            A(J+1)=A(J+1)+R(K)*ARRAY(J,K)
          ENDDO
          A(J+1)=A(J+1)*SIGMA/SIGMAX(J)
          A(1)=A(1)-A(J+1)*XMEAN(J)
        ENDDO

        DO 113 I = 1, NPTS
        YFIT = A(1)
        DO 112 J = 1, DEG
112     YFIT = YFIT + A(J+1)*X(I)**J
113     CHISQ = CHISQ + W(I)*(Y(I)-YFIT)**2
        FREEN = NPTS - NTERMS
        IF( FREEN.EQ.0 ) FREEN = 1.
        CHISQR = CHISQ*WMEAN/FREEN

C       CALCULATE UNCERTAINTIES

        IF( MODE.EQ.WEIGHT ) THEN
          VARNCE = 1./WMEAN
        ELSE
          VARNCE = CHISQR
        ENDIF
        DO 133 J = 1, DEG
        SIGMAA(1+J) = ARRAY(J,J)*VARNCE/(FREE1*SIGMAX(J)**2)
133     SIGMAA(1+J) = SQRT(SIGMAA(1+J))
        SIGMAA(1) = VARNCE / FNPTS
        DO 145 J = 1, DEG
        DO 145 K = 1, DEG
145     SIGMAA(1) = SIGMAA(1) + VARNCE*XMEAN(J)*XMEAN(K)*ARRAY(J,K) /
     1              (FREE1*SIGMAX(J)*SIGMAX(K))

        SIGMAA(1) = SQRT(SIGMAA(1))
        DO 160 J = 1, DEG
        DO 160 K = 1, DEG
160     ARRAY(J,K) = ARRAY(J,K)*VARNCE/(FREE1*SIGMAX(J)*SIGMAX(K))

        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c       real FUNCTION SPOLFT(X,NTERMS,SIGMAA,ARRAY)
c       real X, SIGMAA
c       INTEGER NTERMS
c       INTEGER MAXDEG
c       PARAMETER (MAXDEG=3)
c       real ARRAY(MAXDEG,MAXDEG)

c       INTEGER DEG

c       DEG = NTERMS-1
c       S2YFIT = SIGMAA
c       DO 160 J = 1, DEG
c       S2YFIT = S2YFIT + X**(J+J) * ARRAY(J,J)
c       DO 160 K = J+1, DEG
c160     S2YFIT = S2YFIT + 2.*X**(J+K) * ARRAY(J,K)
c       SPOLFT = SQRT( S2YFIT )
c       RETURN
c       END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        SUBROUTINE MATINVM( NORDER, DET, ARRAY )

        INTEGER MAXDEG
        PARAMETER (MAXDEG=3)
        INTEGER NORDER,IK(MAXDEG),JK(MAXDEG)
        DOUBLE PRECISION DET
        DOUBLE PRECISION ARRAY(MAXDEG,MAXDEG)

        INTEGER I,J,K,L
        DOUBLE PRECISION AMAX,DSAVE

        DET=1.D0
        DO 100 K=1,NORDER

C       FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX

        AMAX=0.D0
21      DO 30 I=K,NORDER
        DO 30 J=K,NORDER
        IF(abs(AMAX).GT.abs(ARRAY(I,J))) GO TO 30
        AMAX=ARRAY(I,J)
        IK(K)=I
        JK(K)=J
30      CONTINUE

C       INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)

        IF(AMAX.NE.0.D0) GO TO 41
        DET=0.D0
        RETURN
41      I=IK(K)
        IF(I-K)21,51,43
43      DO 50 J=1,NORDER
        DSAVE=ARRAY(K,J)
        ARRAY(K,J)=ARRAY(I,J)
50      ARRAY(I,J)=-DSAVE
51      J=JK(K)
        IF(J-K)21,61,53
53      DO 60 I=1,NORDER
        DSAVE=ARRAY(I,K)
        ARRAY(I,K)=ARRAY(I,J)
60      ARRAY(I,J)=-DSAVE

C       ACCUMULATE ELEMENTS OF INVERSE MATRIX

61      DO 70 I=1,NORDER
        IF(I-K)63,70,63
63      ARRAY(I,K)=-ARRAY(I,K)/AMAX
70      CONTINUE
        DO 80 I=1,NORDER
        DO 80 J=1,NORDER
        IF(I-K)74,80,74
74      IF (J-K)75,80,75
75      ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
80      CONTINUE
        DO 90 J=1,NORDER
        IF(J-K)83,90,83
83      ARRAY(K,J)=ARRAY(K,J)/AMAX
90      CONTINUE
        ARRAY(K,K)=1./AMAX
100     DET=DET*AMAX

C       RESTORE ORDENING OF MATRIX

        DO 130 L=1,NORDER
        K= NORDER-L+1
        J=IK(K)
        IF(J-K)111,111,105
105     DO 110 I=1,NORDER
        DSAVE=ARRAY(I,K)
        ARRAY(I,K)=-ARRAY(I,J)
110     ARRAY(I,J)=DSAVE
111     I=JK(K)
        IF(I-K)130,130,113
113     DO 120 J=1,NORDER
        DSAVE=ARRAY(K,J)
        ARRAY(K,J)=-ARRAY(I,J)
120     ARRAY(I,J)=DSAVE
130     CONTINUE
        RETURN
        END
