      subroutine recoilinit
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : June 26, 2007
c | Task  : Initialization of basic recoil information
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c nbbins        : number of bins to calculate maximum recoil energy
c nanglerec     : number of recoil angles
c numangrec     : maximum number of recoil angles
c maxenrec      : number of recoil energies
c numenrec      : maximum number of recoil energies
c ddxrec        : array containing the lab double  differential
c                 xs of the residual nucleus
c ddxrectot     : array containing the total recoil flux in a given
c                 excitation energy bin
c specrecoil    : recoil spectrum
c ddxejlab      : array containing the ddx spectrum of light
c                 particle type in the lab
c Erecmaxmax    : maximum estimated recoil energy of the nucleus
c recoilint     : total recoil integrated over spectrum
c Enrjlabmax    : maximum ejectile energy in lab (classical case)
c enrjlabmax    : absolute maximum ejectile energy
c dangrec,angval: help variables
c angrecmin     : array for lower bin angle values
c angrecmax     : array for upper bin angle values
c cosrecmin     : array for cosine of lower bin angle values
c cosrecmax     : array for cosine of upper bin angle values
c dcosangrec    : array for deltacos
c Emaxinc       : incident lab energy
c Emaxlab       : help variable
c iejlab        : number of ejectile lab bins
c Eejlab        : center of ejectile lab bin
c dEejlab       : width of ejectile lab bin
c Eejlabmax     : maximum energy of ejectile lab bin
c Eejlabmin     : minimum energy of ejectile lab bin
c projectmass   : projectile mass
c compmass      : compound nucleus mass
c PtotCM        : compound initial momentum (or cm momentum)
c Erecinit      : first compound nucleus recoil energy
c etotcm        : compound nucleus total energy (relativistic invariant)
c
      integer   nbbins,numres
      parameter (nbbins=50,numres=5000)
      integer   iang,ierec,iex,in,iz,iej,type,nl
      real      Erecmaxmax(0:numZ,0:numN)
      real      Enrjlabmax(0:numpar)
      real      dangrec,angval
      real      Emaxinc,Emaxlab
      integer   maxentype
      real      Elast
      integer   nend
      real      ehigh
      integer   nhigh
      integer   iejmax
      real      wcos
      real      projectmass,ekinprojlab,PtotCM
      real      compmass
      integer   numZN(numres,numpar)
      real      EexCMmax(numres)
      integer   ires
      integer   iloop
      integer   io1,io2,io3,io4,io5,io6
      integer   if1,if2,if3,if4,if5,if6
      integer   numZcomp,numNcomp
      integer   numZres,numNres
      integer   numZk,numNk
      real      EexCMloc
      integer   k
      integer   numrestot
      real      dEexCM
      real      EexCM(numres,0:nbbins)
      real      PCM(numres,0:nbbins)
      integer   iexc
      real      ECMbin,vCMbin
      integer   ifinal
      real      recoilGSmass
      integer   iexfinal
      real      ratio
      real      EejCM
      real      PrecCM,vrecCM
      real      preclab,vreclab
      real      vejeccm
      real      vejeclab,eejeclab
      real      Erecmaxloc(numres),Erecbin
      real      derec
      integer   irecmaxmax(0:numZ,0:numN)
      real      wnrj
c
c **** Initialization of arrays containing the results in the lab ******
c
      do iex=0,numenrec
        do in=0,numN
          do iz=0,numZ
            Erec(iz,in,iex)=0.
            Erecmin(iz,in,iex)=0.
            Erecmax(iz,in,iex)=0.
          enddo
        enddo
      enddo
      do iang=0,numangrec
        do iej=0,numenrec
          do iex=0,numex
            do in=0,numN
              do iz=0,numZ
                ddxrec(iz,in,iex,iej,iang)=0.
              enddo
            enddo
          enddo
        enddo
      enddo
      do iang=0,2*numangcont+1
        do iej=0,numen2
          do type=0,numpar
            areaejlab(type,iej,iang)=0.
          enddo
        enddo
      enddo
      do iang=0,2*numangrec+1
        do iex=0,numenrec
          do in=0,numN
            do iz=0,numZ
              areareclab(iz,in,iex,iang)=0.
            enddo
          enddo
        enddo
      enddo
      do iex=0,numex
        do in=0,numN
          do iz=0,numZ
            ddxrectot(iz,in,iex)=0.
            specrecoil(iz,in,iex)=0.
          enddo
        enddo
      enddo
      do iang=0,numangcont
        do type=0,6
          do iej=0,numen2
            ddxejlab(type,iej,iang)=0.
          enddo
        enddo
      enddo
      do type=0,6
        do iej=0,numen2
          Eejlab(type,iej)=0.
          Eejlabmin(type,iej)=0.
          Eejlabmax(type,iej)=0.
          dEejlab(type,iej)=0.
        enddo
      enddo
      do in=0,maxN+2
        do iz=0,maxZ+2
          Erecmaxmax(iz,in)=0.
          recoilint(iz,in)=0.
        enddo
      enddo
      do type=0,6
        Enrjlabmax(type)=0.
      enddo
c
c ******** Creation of the angular lab grid for the recoil nuclei ******
c
c kalbachsep: subroutine for separation energy for kalbach systematics
c
      call kalbachsep
      dangrec=180.0/nanglerec
      do 40 iang=0,nanglerec
        angval=dangrec*real(iang)
        angrecmin(iang)=max(0.,angval-0.5*dangrec)
        angrecmax(iang)=min(180.,angval+0.5*dangrec)
        cosrecmin(iang)=cos(angrecmin(iang)*deg2rad)
        cosrecmax(iang)=cos(angrecmax(iang)*deg2rad)
        dcosangrec(iang)=abs(cosrecmin(iang)-cosrecmax(iang))
  40  continue
      do iang=nanglerec+1,2*nanglerec+1
        cosrecmax(iang)=-cosrecmax(iang-nanglerec-1)
        cosrecmin(iang)=-cosrecmin(iang-nanglerec-1)
        dcosangrec(iang)=abs(cosrecmin(iang)-cosrecmax(iang))
      enddo
c
c ********** Creation of the lab energy grid for the ejectiles *********
c
      Emaxinc=Einc+S(0,0,k0)+targetE
      do 50 type=0,6
        if (parskip(type)) goto 50
        Emaxlab=Emaxinc-S(0,0,type)
        Enrjlabmax(type)=Emaxlab
        iejlab(type)=0
        if (Emaxlab.lt.0.) go to 50
        maxentype=0
        do 60 iej=0,maxen
          if (egrid(iej).ge.Emaxlab) then
            maxentype=iej
            go to 70
          endif
   60   continue
   70   if (maxen.gt.0.and.egrid(maxen).lt.Emaxlab) maxentype=maxen
        if (maxentype.gt.0) maxentype=max(maxentype,3)
        nl=Nlast(parZ(type),parN(type),0)
        Elast=eoutdis(type,nl)-elwidth
        Elast=max(Elast,0.)
        if (Elast.gt.egrid(maxentype)) then
          nend=maxentype
        else
          call locate(egrid,0,maxentype,Elast,nend)
        endif
        do 80 iej=1,nend
          Eejlab(type,iej)=egrid(iej)
   80   continue
        ehigh=egrid(maxentype)-egrid(nend)
        nhigh=nint(10.*ehigh)
        nhigh=min(nhigh,numendisc)
        iejlab(type)=nend+nhigh
        do 90 iej=nend+1,nend+nhigh
          Eejlab(type,iej)=egrid(nend)+0.1*(iej-nend)
          if (Eejlab(type,iej).gt.Emaxlab) then
            iejlab(type)=min(iejlab(type),iej)
          endif
   90   continue
        iejmax=iejlab(type)
        do 100 iej=2,iejmax-1
          dEejlab(type,iej)=0.5*(Eejlab(type,iej+1)-Eejlab(type,iej-1))
          Eejlabmax(type,iej)=0.5*(Eejlab(type,iej)+Eejlab(type,iej+1))
          Eejlabmin(type,iej)=0.5*(Eejlab(type,iej)+Eejlab(type,iej-1))
  100   continue
        dEejlab(type,1)=Eejlab(type,1)+
     +                  0.5*(Eejlab(type,2)-Eejlab(type,1))
        Eejlabmin(type,1)=0.
        Eejlabmax(type,1)=dEejlab(type,1)
        Eejlabmin(type,iejmax)=Eejlabmax(type,max(0,iejmax-1))
        Eejlabmax(type,iejmax)=2.*Eejlab(type,iejmax)-
     +                         Eejlabmin(type,iejmax)
        Eejlabmax(type,iejmax)=Emaxlab
        dEejlab(type,iejmax)=Eejlabmax(type,iejmax)-
     +                       Eejlabmin(type,iejmax)
        Eejlab(type,iejmax)=0.5*(Emaxlab+Eejlabmin(type,iejmax))
   50 continue
c
c ************ Creation of the area array for the ejectiles ************
c
      do 110 type=0,6
        if (parskip(type)) goto 110
        do iang=0,2*nanglecont+1
          wcos=dcosangcont(iang)
          do iej=1,iejlab(type)
            areaejlab(type,iej,iang)=wcos*dEejlab(type,iej)
          enddo
          areaejlab(type,0,iang)=0.
        enddo
  110 continue
c
c *************** Calculate maximum excitation energies ****************
c ***************   for all possible residual nuclei    ****************
c
      projectmass=parmass(k0)*amu
      ekinprojlab=Einc
      PtotCM=sqrt(ekinprojlab*(ekinprojlab+2.*projectmass))
      compmass=nucmass(0,0)*amu
      Erecinit=sqrt(PtotCM**2+compmass**2)-compmass
      numZN(1,1)=0
      numZN(1,2)=0
      numZN(1,3)=0
      numZN(1,4)=0
      numZN(1,5)=0
      numZN(1,6)=0
      EexCMmax(1)=Etotal
      numrestot=1
c
c We loop over all possible residual to see if further emission of light
c particles is possible. Starting from a given nucleus defined by
c (numZcomp,numNcomp) we look in loop 150 if the nucleus reached after
c emission of an ejectile (numZres,numNres) has already been obtained
c from a previous emission. If not the number a new residual is reached
c and we have to loop again. If yes, we compare the maximum excitation
c energy obtained from the previous emission with that of the current
c emission and keep the highest.
c
  120 do 130 ires=1,numrestot
        iloop=1
        if (ires.eq.numrestot) iloop=0
        numZcomp=numZN(ires,2)+numZN(ires,3)+numZN(ires,4)+
     +         2*(numZN(ires,5)+numZN(ires,6))
        numNcomp=numZN(ires,1)+numZN(ires,3)+numZN(ires,5)+
     +         2*(numZN(ires,4)+numZN(ires,6))
        do 140 type=1,6
          if1=numZN(ires,1)
          if2=numZN(ires,2)
          if3=numZN(ires,3)
          if4=numZN(ires,4)
          if5=numZN(ires,5)
          if6=numZN(ires,6)
          if (type.eq.1) if1=if1+1
          if (type.eq.2) if2=if2+1
          if (type.eq.3) if3=if3+1
          if (type.eq.4) if4=if4+1
          if (type.eq.5) if5=if5+1
          if (type.eq.6) if6=if6+1
          numZres=numZcomp+parZ(type)
          if (numZres.gt.maxZ+2) go to 140
          if (numZres.gt.Zinit) go to 140
          numNres=numNcomp+parN(type)
          if (numNres.gt.maxN+2) go to 140
          if (numNres.gt.Ninit) go to 140
          EexCMloc=EexCMmax(ires)-S(numZcomp,numNcomp,type)
          if (EexCMloc.le.0.) go to 140
          do 150 k=1,numrestot
            io1=numZN(k,1)
            io2=numZN(k,2)
            io3=numZN(k,3)
            io4=numZN(k,4)
            io5=numZN(k,5)
            io6=numZN(k,6)
            numZk=numZN(k,2)+numZN(k,3)+numZN(k,4)+
     +             2*(numZN(k,5)+numZN(k,6))
            numNk=numZN(k,1)+numZN(k,3)+numZN(k,5)+
     +         2*(numZN(k,4)+numZN(k,6))
            if ((numZk.eq.numZres).and.(numNk.eq.numNres)) then
              EexCMmax(k)=max(EexCMmax(k),EexCMloc)
              go to 140
            endif
  150     continue
c A new residual is reached
          numrestot=numrestot+1
          iloop=1
          do k=1,6
            numZN(numrestot,k)=numZN(ires,k)
          enddo
          numZN(numrestot,type)=numZN(ires,type)+1
          EexCMmax(numrestot)=EexCMloc
  140   continue
  130 continue
      if (iloop.eq.1) go to 120
c
c ******** Define excitation energy bins for all residual nuclei *******
c
      do ires=1,numrestot
        dEexCM=EexCMmax(ires)/nbbins
        do iex=0,nbbins
          EexCM(ires,iex)=min(iex*dEexCM,EexCMmax(ires))
          EexCM(ires,iex)=max(iex*dEexCM,0.)
          PCM(ires,iex)=0.
        enddo
      enddo
c
c ********      Loop over all reactions paths to determine      ********
c ******** The maximum recoil energy possible for each residual ********
c
      PCM(1,nbbins)=PtotCM
      do 180 ires=1,numrestot
        numZcomp=numZN(ires,2)+numZN(ires,3)+numZN(ires,4)+
     +         2*(numZN(ires,5)+numZN(ires,6))
        numNcomp=numZN(ires,1)+numZN(ires,3)+numZN(ires,5)+
     +         2*(numZN(ires,4)+numZN(ires,6))
        compmass=nucmass(numZcomp,numNcomp)*amu
        do 190 iex=nbbins,0,-1
          ECMbin=EexCM(ires,iex)
          vCMbin=PCM(ires,iex)/(compmass+ECMbin)
          do 200 type=0,6
            ejectmass=parmass(type)*amu
            if1=numZN(ires,1)
            if2=numZN(ires,2)
            if3=numZN(ires,3)
            if4=numZN(ires,4)
            if5=numZN(ires,5)
            if6=numZN(ires,6)
            if (type.eq.1) if1=if1+1
            if (type.eq.2) if2=if2+1
            if (type.eq.3) if3=if3+1
            if (type.eq.4) if4=if4+1
            if (type.eq.5) if5=if5+1
            if (type.eq.6) if6=if6+1
            do 210 k=1,numrestot
              io1=numZN(k,1)
              io2=numZN(k,2)
              io3=numZN(k,3)
              io4=numZN(k,4)
              io5=numZN(k,5)
              io6=numZN(k,6)
              if ((io1.eq.if1).and.(io2.eq.if2).and.(io3.eq.if3).and.
     +            (io4.eq.if4).and.(io5.eq.if5).and.(io6.eq.if6)) then
                ifinal=k
                numZres=numZN(k,2)+numZN(k,3)+numZN(k,4)+
     +                2*(numZN(k,5)+numZN(k,6))
                numNres=numZN(k,1)+numZN(k,3)+numZN(k,5)+
     +                2*(numZN(k,4)+numZN(k,6))
                go to 220
              endif
  210       continue
c
c If the residual nucleus is not among the considered ones it means
c that it cannot be obtained because all possible residuals have been
c defined in loop 130. We thus consider another ejectile.
c
            go to 200
c
c Loop over the various excitation energies bins of the residual
c nucleus to deduce the lab recoil energy for these bins
c
  220       recoilGSmass=nucmass(numZres,numNres)*amu
            do 230 iexfinal=nbbins,0,-1
              recoilmass=recoilGSmass+EexCM(ifinal,iexfinal)
              if (type.ne.0) then
                ratio=recoilmass/(ejectmass*(ejectmass+recoilmass))
              endif
c
c Determine ejectile energy to check if emission can occur
c
              EejCM=ECMbin-EexCM(ifinal,iexfinal)
              EejCM=EejCM-S(numZcomp,numNcomp,type)
              if (abs(EejCM).le.1.0e-10) EejCM=0.
              if (EejCM.le.0.) go to 230
              if (type.ne.0) then
                  PrecCM=ejectmass*sqrt(2*ratio*EejCM)
                else
                  PrecCM=EejCM
              endif
              vrecCM=PrecCM/recoilmass
              vreclab=vCMbin+vrecCM
              preclab=recoilmass*vreclab
              PCM(ifinal,iexfinal)=max(PCM(ifinal,iexfinal),preclab)
              if (type.ne.0) then
                  vejeccm=PrecCM/ejectmass
                  vejeclab=vCMbin+vejeccm
                  eejeclab=0.5*ejectmass*vejeclab**2
                  Enrjlabmax(type)=max(Enrjlabmax(type),eejeclab)
                else
                  eejeclab=EejCM
                  Enrjlabmax(type)=max(Enrjlabmax(type),eejeclab)
              endif
  230       continue
  200     continue
  190   continue
  180 continue
c
c Maximum recoil energy for each residual
c
      do 240 ires=1,numrestot
        Erecmaxloc(ires)=0.
        numZres=numZN(ires,2)+numZN(ires,3)+numZN(ires,4)+
     +        2*(numZN(ires,5)+numZN(ires,6))
        numNres=numZN(ires,1)+numZN(ires,3)+numZN(ires,5)+
     +        2*(numZN(ires,4)+numZN(ires,6))
        recoilGSmass=nucmass(numZres,numNres)*amu
        do 250 iexc=0,nbbins
          recoilmass=recoilGSmass+EexCM(ires,iexc)
          PrecCM=PCM(ires,iexc)
          Erecbin=PrecCM**2/(2.*recoilmass)
          Erecmaxloc(ires)=max(Erecmaxloc(ires),Erecbin)
  250   continue
  240 continue
c
c Finally the recoils energy grids are defined and maximum excitation
c energies are stored
c
      do ires=1,numrestot
        io1=numZN(ires,1)
        io2=numZN(ires,2)
        io3=numZN(ires,3)
        io4=numZN(ires,4)
        io5=numZN(ires,5)
        io6=numZN(ires,6)
        numZres=numZN(ires,2)+numZN(ires,3)+numZN(ires,4)+
     +        2*(numZN(ires,5)+numZN(ires,6))
        numNres=numZN(ires,1)+numZN(ires,3)+numZN(ires,5)+
     +        2*(numZN(ires,4)+numZN(ires,6))
        Erecmaxmax(numZres,numNres)=Erecmaxloc(ires)
        irecmaxmax(numZres,numNres)=ires
      enddo
      do 460 iz=0,maxZ+2
        do 470 in=0,maxN+2
          if (Erecmaxmax(iz,in).le.0.) go to 470
          derec=Erecmaxmax(iz,in)/(maxenrec+1)
          ires=irecmaxmax(iz,in)
          if (ires.eq.0) go to 470
          do 480 iexc=0,numenrec
            Erecmin(iz,in,iexc)=iexc*derec
            Erecmax(iz,in,iexc)=(iexc+1)*derec
            Erec(iz,in,iexc)=iexc*derec+0.5*derec
 480      continue
 470    continue
 460  continue
c
c ************* Creation of the area array for the recoils *************
c
      do 500 iz=0,maxZ+2
        do 510 in=0,maxN+2
          if (Erecmaxmax(iz,in).eq.0.) go to 510
          do 520 ierec=0,maxenrec
            wnrj=Erecmax(iz,in,ierec)-Erecmin(iz,in,ierec)
            do 530 iang=0,2*nanglerec+1
              wcos=dcosangrec(iang)
              areareclab(iz,in,ierec,iang)=wnrj*wcos
  530       continue
  520     continue
  510   continue
  500 continue
c
c ******* Feeding of the first compound nucleus recoil ddx array  ******
c
      derec=Erecmaxmax(0,0)/(maxenrec+1)
      irecinit=min(int(Erecinit/derec),numenrec)
      if (areareclab(0,0,irecinit,0).ne.0.)
     +  ddxrec(0,0,maxex(0,0),irecinit,0)=xstotinc/twopi/
     +    areareclab(0,0,irecinit,0)
      ddxrectot(0,0,maxex(0,0))=xstotinc/twopi
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
