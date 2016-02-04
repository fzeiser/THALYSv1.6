      subroutine brosafy(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 26, 2006
c | Task  : Fission fragment yields based on Brosa model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist
      character*4      gschar
      character*8      filen
      character*90     gsfile,precfile,barfile
      integer          Zix,Nix,Z,A,Zbrosa,amassmax,numoff,k,i,massdif,
     +                 index1,index2,amassar(1:7),amassdum,iloop,
     +                 numtempsl,numtempst,numtempst2
      double precision transm,sl,st,st2,stot,trcof
      real             mn,mp,noff(8),bindgs(1:9,1:2),binddum,bar,width,
     +                 rmass,temps(9),Tmpmax,Tmp,crel,elsc,hm,ET,BFT,
     +                 HMT,ELT,fraction,Edefo,somtot,bindsc,ignatyuk,ald
      real             bfinter(1:9,1:2),hwinter(1:9,1:2),bf_sl(9),
     +                 bf_st(9),bf_st2(9),bfsplin_sl(9),bfsplin_st(9),
     +                 bfsplin_st2(9),hw_sl(9),hw_st(9)
      logical dont
      real fmass(nummass), fmass_sl(nummass), fmass_st(nummass),
     +     fmass_st2(nummass),fmasscor_sl(nummass),fmasscor_st(nummass),
     +     fmasscor_st2(nummass),fmasscor(nummass)
      real fmz(nummass,numelem),fmzcor(nummass,numelem),
     +     fmz_sl(nummass,numelem),fmz_st(nummass,numelem),
     +     fmz_st2(nummass,numelem),fmzcor_sl(nummass,numelem),
     +     fmzcor_st(nummass,numelem),fmzcor_st2(nummass,numelem)
      real hmsplin(9),elsplin(9), Esplin(9)
      real hmneck_sl(9),hmneck_st(9), elneck_sl(9),elneck_st(9),
     +     Eneck_sl(9), Eneck_st(9)
      real Eneck(9), hmneck(9), elneck(9)
      real Eneckinter(9,2), hmneckinter(9,2), elneckinter(9,2)
      data temps/0.0,0.3,0.6,0.9,1.2,1.6,2.0,2.5,3.0/
c
c Zix            : charge number index for residual nucleus
c Nix            : neutron number index for residual nucleus
c ZZ,Z           : charge number of residual nucleus
c AA,A           : mass number of residual nucleus
c Zbrosa         : charge number within range of Brosa tables
c parmass        : mass of particle in a.m.u.
c mn             : neutron mass in MeV
c mp             : proton mass in MeV
c noff           : number of neutrons away from heaviest isotope given
c                  in structure data base of a certain element
c numoff         : number of isotopes for which information is present
c                  in structure data base for a certain element
c bf_sl/st/st2   : SL, ST I and ST II outer fission barrier heights
c bfsplin        : barrier height splin fit parameters
c bfsplin_sl,....: SL barrier height splin fit parameters
c bindgs         : ground-state binding energy
c gschar         : structure database ground-state binding energy file
c gsfile         : full path for structure database file
c
c     constants
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
c
c Brosa parameters available between Z=72 and Z=96, outside this range
c we adopt the values of the boundary nuclides
c
      Zbrosa=max(72,min(Z,96))
      mn=parmass(1)*amu
      mp=parmass(2)*amu
c
c     up to Np(Z=93) 7 isotopes per element calculated, above 6.
c
      noff(1)=0
      noff(2)=2
      noff(3)=4
      noff(4)=6
      noff(5)=8
      if (Z.lt.94) then
        noff(6)=11
        noff(7)=14
        noff(8)=0
        numoff=7
      else
        noff(6)=10
        noff(7)=0
        noff(8)=0
        numoff=6
      endif
c
c     initialize
c
      do 5, k=1,9
        bf_sl(k)=0.
        bf_st(k)=0.
        bf_st2(k)=0.
        bfsplin_sl(k)=0.
        bfsplin_st(k)=0.
        bfsplin_st2(k)=0.
 5    continue
      do 7, k=1,9
         do 7, i=1,2
         bindgs(k,i)=0.
 7    continue
c
c     reading ground state binding energies
c
      gschar='z000'
      write(gschar(2:4),'(i3.3)') Zbrosa
      gsfile=path(1:lenpath)//'fission/brosa/groundstate/'//gschar
      open (unit=2,status='old',file=gsfile)
c
c amassmax: heaviest isotope for which parameters have been calculated
c amassar : mass array
c massdif : mass difference between amassmax and the nucleus studied
c index1,2: denote array elements used to interpolate
c
c      determine position nucleus in array
c
      read(2,'(4x,i4)') amassmax
      rewind(2)
      do 8,k=1,numoff
        amassar(k)=amassmax-noff(k)
 8    continue
      massdif=amassmax-A
      do 9, k=1,numoff
         if(massdif.GE.noff(k).and.massdif.LE.noff(k+1))then
            if(massdif.EQ.noff(k))then
               index1=k
               index2=-1
            else
               if(massdif.EQ.noff(k+1))then
                  index1=k+1
                  index2=-1
               else
                  index1=k
                  index2=k+1
               endif
            endif
            goto 9
         endif
         if(massdif.GT.noff(numoff))then
            index1=numoff
            index2=-1
         else
            if(massdif.LT.noff(1))then
               index1=1
               index2=-1
            endif
         endif
 9    continue
c
c amassdum: dummy mass variable
c binddum : dummy binding energy variable
c hwsplin : barrier width splin fit parameters
c hw      : barrier width
c bfinter : barrier heights used for interpolation
c hwinter : barrier widths used for interpolation
c filen   : file name
c barfile : path of barrier parameter file in structure data base
c lexist  : logical to determine existence of barfile
c numtemp : running variable denoting the numtempth temperature
c tmpmax  : maximal temperatures for which structure data is available
c
c     read in ground state binding energy as function of T
c
      i=1
      k=1
      if(index2.EQ.-1)then
 10      read(2,'(4x,i4,15x,f15.5)', end=11) amassdum,binddum
         if(amassdum.EQ.amassar(index1))then
            bindgs(i,1)=abs(binddum)
            i=i+1
         endif
         goto 10
 11      continue
         close (unit=2)
      else
 12      read(2,'(4x,i4,15x,f15.5)', end=13) amassdum, binddum
         if(amassdum.EQ.amassar(index1))then
            bindgs(i,1)=abs(binddum)
            i=i+1
         else
            if(amassdum.EQ.amassar(index2))then
               bindgs(k,2)=abs(binddum)
               k=k+1
            endif
         endif
         goto 12
 13      continue
         close (unit=2)
      endif
c
c      superlong (sl), standard 1 (st), standard 2 (st2) loop
c
      do 1000, iloop=1,3
         do 18, k=1,9
            bfsplin(k)=0.
            hwsplin(k)=0.
            bf(k)=0.
            hw(k)=0.
            do 18, i=1,2
               bfinter(k,i)=0.
               hwinter(k,i)=0.
 18      continue
         if(iloop.EQ.1)then
            filen='z000.sl '
         else
            if(iloop.EQ.2)then
               filen='z000.st '
            else
               filen='z000.st2'
            endif
         endif
         write(filen(2:4),'(i3.3)') Zbrosa
         barfile=path(1:lenpath)//'fission/brosa/barrier/'//filen
         inquire (file=barfile,exist=lexist)
         if (lexist)then
           open (unit=10,status='old',file=barfile)
c
c     if fission modes exists:
c     read and calculate barrier parameters height, width, and position
c     (in calculating the width, the mass at T=0 is used, since
c     the dependence on T is very small)
c
            i=1
            k=1
 20         read(10,'(4x,i4,15x,2f15.5)',end=21) amassdum,bar,width
            if(amassdum.EQ.amassar(index1))then
               rmass=(1/(clight*clight))*(Z*mp+(A-Z)*mn-
     +              bindgs(1,1))
               bfinter(i,1)=bar
               hwinter(i,1)=hbar*sqrt(width*1.D+30/rmass)
               i=i+1
            endif
            if(index2.NE.-1)then
               if(amassdum.EQ.amassar(index2))then
                  rmass=(1/(clight*clight))*(Z*mp+(A-Z)*mn-
     +                 bindgs(1,2))
                  bfinter(k,2)=bar
                  hwinter(k,2)=hbar*sqrt(width*1.D+30/rmass)
                  k=k+1
               endif
            endif
            goto 20
 21         continue
           close (unit=10)
         endif
c
c      interpolate, if necessary, barrier parameters
c
         numtemp=1
         if(index2.EQ.-1)then
            do 31, k=1,9
               if(bfinter(k,1).NE.0.)then
                  bf(k)=bfinter(k,1)
                  hw(k)=hwinter(k,1)
                  numtemp=k
               endif
 31         continue
         else
            if(abs(noff(index1)-noff(index2)).EQ.3)then
               do 32, k=1,9
                  if(bfinter(k,1).NE.0.and.bfinter(k,2).NE.0.)then
                     massdif=abs(amassar(index1)-A)
                     bf(k)= bfinter(k,1)+
     +                    (massdif/3.)*(bfinter(k,2)-bfinter(k,1))
                     hw(k)= hwinter(k,1)+
     +                    (massdif/3.)*(hwinter(k,2)-hwinter(k,1))
                     numtemp=k
                  endif
 32            continue
            else
               do 33, k=1,9
                  if(bfinter(k,1).NE.0.and.bfinter(k,2).NE.0.)then
                     bf(k)=0.5*(bfinter(k,1)+bfinter(k,2))
                     hw(k)=0.5*(hwinter(k,1)+hwinter(k,2))
                     numtemp=k
                  endif
 33            continue
            endif
         endif
c
c      fit spline to barrier parameters
c
         if(iloop.EQ.2)then
            bf(numtemp+1)=bf_sl(numtemp+1)
            hw(numtemp+1)=hw_sl(numtemp+1)
            numtemp=numtemp+1
         endif
         if(iloop.EQ.3)then
            bf(numtemp+1)=bf_st(numtemp+1)
            hw(numtemp+1)=hw_st(numtemp+1)
            numtemp=numtemp+1
         endif
         tmpmax=temps(numtemp)
         if(bf(1).eq.0)numtemp=9
         call spline(temps,bf,numtemp,2.e+30,2.e+30,bfsplin)
         call spline(temps,hw,numtemp,2.e+30,2.e+30,hwsplin)
c
c ald      : level density parameter
c ignatyuk : function for energy dependent level density parameter
c Tmp      : temperature
c excfis   : excitation energy at fission
c trans    : subroutine to determine transmission coefficients
c            per fission mode
c trof     : transmission coefficients per fission mode
c sl,st,st2: SL,ST I, ST II transmission coefficients
c stot     : sum of sl, st, and st2
c
c      calculate transmission coefficient trcof
c
         ald=ignatyuk(Zix,Nix,excfis,0)
         Tmp=sqrt(excfis/ald)
c
c      check existence fission mode by looking at bf(1)
c
         if(Tmp.lt.tmpmax.and.bf(1).NE.0.)then
            call trans(Zix,Nix,transm)
            trcof=transm
         else
            trcof=0.d0
         endif
c
c      fill arrays with trcof, bf, bfsplin
c
         if(iloop.EQ.1)then
            sl=trcof
            numtempsl=numtemp
            do 105, k=1,numtemp
               bf_sl(k)=bf(k)
               hw_sl(k)=hw(k)
               bfsplin_sl(k)=bfsplin(k)
 105        continue
         else
            if(iloop.EQ.2)then
               st=trcof
               numtempst=numtemp
               do 106, k=1,numtemp
                  bf_st(k)=bf(k)
                  hw_st(k)=hw(k)
                  bfsplin_st(k)=bfsplin(k)
 106           continue
            else
               st2=trcof
               numtempst2=numtemp
               do 107, k=1,numtemp
                  bf_st2(k)=bf(k)
                  bfsplin_st2(k)=bfsplin(k)
 107           continue
            endif
         endif
c
c      end first loop over fission modes sl, st, st2
c
 1000 continue
c
c      final transmission coefficients
c
      stot=sl+st+st2
      if(stot.EQ.0.)then
         sl=1.
         st=0.
         st2=0.
      else
         sl=sl/stot
         st=st/stot
         st2=st2/stot
      endif
c
c crel       : scaling factor for neck curvature
c dont       : logical for mass yield calculation
c fmass      : fission fragment mass yield
c fmasscor   : corrected fission fragment mass yield
c fmz        : fission fragment isotope yield
c fmzcor     : corrected fission fragment isotope yield
c hmneck     : mean mass heavy fragment at scission point
c elneck     : nucleus half length at scission point
c Eneck      : prescission energy
c hmneckinter: mean mass heavy fragment at scission point used
c              for interpolation
c elneckinter: nucleus half length at scission point used for
c              interpolation
c Eneckinter : prescission energy used for interpolation
c hmsplin    : mean mass heavy fragment at scission point
c              splin fit parameters
c elsplin    : nucleus half length at scission point splin
c              fit parameters
c Esplin     : prescission energy splin fit parameters
c precfile   : path with prescission shape info in structure data base
c Edefo      : excitation energy at scission
c et         : energy gain at scission with respect to ground state
c BFT        : temperature-dependent barrier height
c fraction   : fraction of potential energy gain available as internal
c              excitation energy at scission
c neck       : subroutine to determine mass and isotope yield per
c              fission mode
c
c
c      MASS DISTRIBUTION
c      second loop over fission modes starts here
c
      do 2000, iloop=1,3
         do 201, k=1,9
            bf(k)=0.
            bfsplin(k)=0.
 201     continue
         if(iloop.EQ.1)then
            dont=.false.
            numtemp=numtempsl
            if(sl.eq.0.)dont=.true.
            do 203, k=1,numtemp
               bf(k)=bf_sl(k)
               bfsplin(k)=bfsplin_sl(k)
 203        continue
            crel=0.1
            filen='z000.sl '
         else
            if(iloop.EQ.2)then
               dont=.false.
               numtemp=numtempst
               if(st.eq.0.)dont=.true.
               do 204, k=1,numtemp
                  bf(k)=bf_st(k)
                  bfsplin(k)=bfsplin_st(k)
 204           continue
               crel=0.3
               filen='z000.st '
            else
               dont=.false.
               numtemp=numtempst2
               if(st2.eq.0.)dont=.true.
               do 206, k=1,numtemp
                  bf(k)=bf_st2(k)
                  bfsplin(k)=bfsplin_st2(k)
 206           continue
               crel=0.3
               filen='z000.st2'
            endif
         endif
         write(filen(2:4),'(i3.3)') Zbrosa
c
c      initialize arrays for mass and charge yields
c
         do 400, k=1,nummass
            fmass(k)=0.
            fmasscor(k)=0.
            do 400, i=1,numelem
               fmz(k,i)=0.
               fmzcor(k,i)=0.
 400     continue
c
c      initialize arrays for neck parameters
c
         do 270, k=1,9
            hmneck(k)=0.
            elneck(k)=0.
            Eneck(k)=0.
            do 270, i=1,2
               hmneckinter(k,i)=0.
               elneckinter(k,i)=0.
               Eneckinter(k,i)=0.
 270     continue
         do 300, k=1, numtemp
            if(iloop.EQ.2.and.k.EQ.numtemp)then
               hmneck(k)=hmneck_sl(k)
               elneck(k)=elneck_sl(k)
               Eneck(k)=Eneck_sl(k)
               goto 300
            endif
            if(iloop.EQ.3.and.k.EQ.numtemp)then
               hmneck(k)=hmneck_st(k)
               elneck(k)=elneck_st(k)
               Eneck(k)=Eneck_st(k)
               goto 300
            endif
 300     continue
         if(dont)goto 17998
c
c  open file with prescission output and read values heavy fragment mass
c  prescission energy, corresponding nucleus half length, and mass heavy
c  fragment
c
       precfile=path(1:lenpath)//'fission/brosa/prescission/'//filen
       open (unit=10,status='old',file=precfile)
         i=1
         k=1
 310     read(10,'(4x,i4,15x,3f15.5)',end=311) amassdum,hm,elsc,bindsc
         if(amassdum.EQ.amassar(index1))then
            hmneckinter(i,1)=hm
            elneckinter(i,1)=elsc
            Eneckinter(i,1)=abs(bindsc)-bindgs(i,1)
            i=i+1
         endif
         if(index2.NE.-1)then
            if(amassdum.EQ.amassar(index2))then
               hmneckinter(k,2)=hm
               elneckinter(k,2)=elsc
               Eneckinter(k,2)=abs(bindsc)-bindgs(k,2)
               k=k+1
            endif
         endif
         goto 310
 311     continue
         close (unit=10)
c
c     interpolate, if necessary, prescission shape parameters
c
         if(index2.EQ.-1)then
            do 331, k=1,9
               if(hmneckinter(k,1).NE.0.)then
                  hmneck(k)=hmneckinter(k,1)*A/amassar(index1)
                  elneck(k)=elneckinter(k,1)
                  Eneck(k)=Eneckinter(k,1)
               endif
 331        continue
         else
            if(abs(noff(index1)-noff(index2)).EQ.3)then
               do 332, k=1,9
                  if(hmneckinter(k,1).NE.0.and.
     +                 hmneckinter(k,2).NE.0.)then
                     massdif=abs(amassar(index1)-A)
                     hmneck(k)= hmneckinter(k,1)+
     +                    (massdif/3.)*(hmneckinter(k,2)
     +                    -hmneckinter(k,1))
                     elneck(k)= elneckinter(k,1)+
     +                    (massdif/3.)*(elneckinter(k,2)
     +                    -elneckinter(k,1))
                     Eneck(k)= Eneckinter(k,1)+
     +                    (massdif/3.)*(Eneckinter(k,2)
     +                    -Eneckinter(k,1))
                  endif
 332            continue
            else
               do 333, k=1,9
                  if(hmneckinter(k,1).NE.0.and.
     +                 hmneckinter(k,2).NE.0.)then
                    hmneck(k)=0.5*(hmneckinter(k,1)+hmneckinter(k,2))
                    elneck(k)=0.5*(elneckinter(k,1)+elneckinter(k,2))
                    Eneck(k)=0.5*(Eneckinter(k,1)+Eneckinter(k,2))
                  endif
 333           continue
            endif
         endif
c
c      spline fitting to neck input parameters
c
         call spline(temps,hmneck,numtemp,2.e+30,2.e+30,hmsplin)
         call spline(temps,elneck,numtemp,2.e+30,2.e+30,elsplin)
         call spline(temps,Eneck,numtemp,2.e+30,2.e+30,Esplin)
c
c      calculate excitation energy at scission (including an iteration
c      to calculate Edefo and ET consistently, one iteration suffices)
c
         ald=ignatyuk(Zix,Nix,excfis,0)
         Tmp=sqrt(excfis/ald)
         call splint(temps,Eneck, Esplin,numtemp,Tmp,ET)
         call splint(temps, bf, bfsplin,numtemp,Tmp,BFT)
         fraction=1.
         if(excfis.GT.BFT)then
            Edefo=fraction*(ET+BFT)+excfis-BFT
            if(ET.LE.0.)Edefo=excfis
            ald=ignatyuk(Zix,Nix,Edefo,0)
            Tmp=sqrt(Edefo/ald)
            call splint(temps,Eneck, Esplin,numtemp,Tmp,ET)
            Edefo=fraction*(ET+BFT)+excfis-BFT
            if(ET.LE.0.)Edefo=excfis
         else
            Edefo=(ET+excfis)*fraction
            if(ET.LE.0.)Edefo=excfis
            ald=ignatyuk(Zix,Nix,Edefo,0)
            Tmp=sqrt(Edefo/ald)
            call splint(temps,Eneck, Esplin,numtemp,Tmp,ET)
            Edefo=(ET+excfis)*fraction
            if(ET.LE.0.)Edefo=excfis
         endif
c
c    check if Tmp(Edefo) does not exceed the highest value for which the
c    barrier is defined
c
         ald=ignatyuk(Zix,Nix,Edefo,0)
         Tmp=sqrt(Edefo/ald)
         if(Tmp.GT.temps(numtemp))then
            ald=ignatyuk(Zix,Nix,excfis,0)
            Tmp=sqrt(excfis/ald)
            if(Tmp.LT.temps(numtemp))then
               call splint(temps,hmneck, hmsplin,numtemp,Tmp,HMT)
               call splint(temps,elneck, elsplin,numtemp,Tmp,ELT)
               call neck(Z,A,fmass,fmasscor,fmz,fmzcor,HMT
     +              ,Edefo,ELT,crel)
            endif
         else
            call splint(temps,hmneck, hmsplin,numtemp,Tmp,HMT)
            call splint(temps,elneck, elsplin,numtemp,Tmp,ELT)
            call neck(Z,A,fmass,fmasscor,fmz,fmzcor,HMT
     +           ,Edefo,ELT,crel)
         endif
17998    continue
         if(iloop.EQ.1)then
            do 500,k=1,nummass
               fmass_sl(k)=fmass(k)*sl
               fmasscor_sl(k)=fmasscor(k)*sl
               do 500, i=1,numelem
                  fmz_sl(k,i)=fmz(k,i)*sl
                  fmzcor_sl(k,i)=fmzcor(k,i)*sl
 500        continue
            do 501, k=1,numtemp
               hmneck_sl(k)=hmneck(k)
               elneck_sl(k)=elneck(k)
               Eneck_sl(k)=Eneck(k)
 501        continue
         else
            if(iloop.EQ.2)then
               do 505,k=1,nummass
                  fmass_st(k)=fmass(k)*st
                  fmasscor_st(k)=fmasscor(k)*st
               do 505, i=1,numelem
                  fmz_st(k,i)=fmz(k,i)*st
                  fmzcor_st(k,i)=fmzcor(k,i)*st
 505           continue
               do 506, k=1,numtemp
                  hmneck_st(k)=hmneck(k)
                  elneck_st(k)=elneck(k)
                  Eneck_st(k)=Eneck(k)
 506           continue
            else
               do 510,k=1,nummass
                  fmass_st2(k)=fmass(k)*st2
                  fmasscor_st2(k)=fmasscor(k)*st2
               do 510, i=1,numelem
                  fmz_st2(k,i)=fmz(k,i)*st2
                  fmzcor_st2(k,i)=fmzcor(k,i)*st2
 510           continue
            endif
         endif
c
c      end loop over sl, st and st2 modes
c
 2000 continue
c
c disa      : normalised fission fragment mass yield per excitation
c             energy bin
c disacor   : normalised fission product mass yield per excitation
c             energy bin
c disaz     : normalised fission fragment isotope yield
c             per excitation energy bin
c disazcor  : normalised fission product isotope yield
c             per excitation energy bin
c flagffevap: flag for calculation of particle evaporation from
c             fission fragment mass yields
c
      somtot=1.
      do 4000,k=1,nummass
         disa(k)=(fmass_sl(k)+fmass_st(k)+fmass_st2(k))/somtot
         if (flagffevap) then
            disacor(k)=(fmasscor_sl(k)+fmasscor_st(k)+
     +           fmasscor_st2(k))/somtot
         else
            disacor(k)=0.
         endif
         do 4000, i=1,numelem
         disaz(k,i)=(fmz_sl(k,i)+fmz_st(k,i)+fmz_st2(k,i))/somtot
         if(flagffevap)then
            disazcor(k,i)=(fmzcor_sl(k,i)+fmzcor_st(k,i)+
     +           fmzcor_st2(k,i))/somtot
         else
            disazcor(k,i)=0.
         endif
 4000 continue
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
