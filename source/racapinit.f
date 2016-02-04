      subroutine racapinit
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely, Xu Yi, and Arjan Koning
c | Date  : December 13, 2013
c | Task  : Initialization of radiative capture model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer i,j,acp,ka,jlev,i0,iz,ia
      integer phnumracap
      integer ptclp,holep,ptcln,holen
      integer xprty,wprty
      integer parity
      real phdens2,Eldpd,Sldpd
      real racapdamp,ignatyuk
      real Ephracap,sigparleldenp,sigparleldenn,surfwellE
      real e,spinf,sf,eopt
      double precision density,ldmdposj,ldmdnegj,ldmdpos,ldmdneg
      double precision phmdpos,phmdposj,phmdnegj
      double precision phjpmdpos,phjpmdneg,phjpmdposj,phjpmdnegj
      logical surfwellgo,lexist
      integer nexphjp,nex
      real eb, ee, Egridphjp(0:numdens)
      double precision ldbpos, ldbneg, ldepos, ldeneg
      double precision lldbpos,lldepos,lldbneg,lldeneg
      double precision ldtabpos,ldtabneg
      character*100 filespec
      dimension ldmdposj(numdensracap,0:numJph),
     +          ldmdnegj(numdensracap,0:numJph),
     +          phmdposj(numdensracap,0:numJph),
     +          phmdnegj(numdensracap,0:numJph),
     +          phjpmdposj(numdensracap,0:numJph),
     +          phjpmdnegj(numdensracap,0:numJph)
      dimension ldmdpos(numdensracap),ldmdneg(numdensracap)
ctest dimension phmdpos(numdensracap),phmdneg(numdensracap)
      dimension phmdpos(numdensracap)
      dimension phjpmdpos(numdensracap),phjpmdneg(numdensracap)
      dimension Eldpd(numdensracap),Sldpd(0:numJph)
c
c ********************** General initializations ***********************
c
c     accept the level density and ph level density, with simple calculation and gridlization.
c
c     set the grid of spin, for odd A and even A:
c
      if (mod(Ainit,2).eq.1) then
        do j=0,numJph
          Sldpd(j)=real(j)+0.5
          Sldpd(j)=real(Sldpd(j))
        end do
      end if
      if (mod(Ainit,2).eq.0) then
        do j=0,numJph
          Sldpd(j)=real(j)
          Sldpd(j)=real(Sldpd(j))
        end do
      end if
c
c     set the grid of exitation energy points: 0.125(point1), 0.375(point2)  ...
c                             (bin: 0.25) 0___point1___0.25___point2___0.5   ...
c
      do i=1,numdensracap
        Eldpd(i)=0.125+(dble(i)-1.0)*0.25
        Eldpd(i)=real(Eldpd(i))
      enddo
c
c     choose the NLD model for racap calculation (ldmodelracap=1, 2, or 3)
c
c     ldmodelracap=1: using the spin- parity- dependent ph-NLD model (pphh)
c
c     edens and nendens need to be read from phdensitytablejp.f. Therefore, if ldmodelracap
c     need to be set as a keyword in the input file, phdensitytablejp.f should be changed,
c     and talys.cmb should be changed accordingly.

      do iz=0,numZ
        do ia=0,numN
          nlevracap(iz,ia)=0
        enddo
      enddo
      if (ldmodelracap.eq.1) then
        do iz=0,numZ
          do ia=0,numN
            do nex=0,numdens
              phdenstot(iz,ia,nex)=0.
              edensphjp(iz,ia,nex)=0.
              do j=0,numJph
                do i=-1,1
                  phdensjp(iz,ia,nex,j,i)=0.
                enddo
              enddo
            enddo
          enddo
        enddo
        call phdensitytablejp(0,0)
      endif
c
c ldmodelracap may have been changed in phdensitytablejp
c
      if (ldmodelracap.eq.1) then
        do nex=1,nendens(0,0)
          Egridphjp(nex)=real(edensphjp(0,0,nex))
        enddo
        Egridphjp(0)=0.0
        do j=0,numJph
          do i=1,numdensracap
            if (Eldpd(i).lt.Egridphjp(nendens(0,0))) then
              call locate(Egridphjp,0,nendens(0,0),Eldpd(i),nexphjp)
              eb=edensphjp(0,0,nexphjp)
              ee=edensphjp(0,0,nexphjp+1)
              ldbpos=phdensjp(0,0,nexphjp,j,1)
              ldbneg=phdensjp(0,0,nexphjp,j,-1)
              ldepos=phdensjp(0,0,nexphjp+1,j,1)
              ldeneg=phdensjp(0,0,nexphjp+1,j,-1)
            else
              eb=edensphjp(0,0,nendens(0,0)-1)
              ee=edensphjp(0,0,nendens(0,0))
              ldbpos=phdensjp(0,0,nendens(0,0)-1,j,1)
              ldbneg=phdensjp(0,0,nendens(0,0)-1,j,-1)
              ldepos=phdensjp(0,0,nendens(0,0),j,1)
              ldeneg=phdensjp(0,0,nendens(0,0),j,-1)
            endif
            if (ldbpos.gt.1..and.ldepos.gt.1.) then
              lldbpos=log(ldbpos)
              lldepos=log(ldepos)
              ldtabpos=exp(lldbpos+(Eldpd(i)-eb)/(ee-eb)*
     &                 (lldepos-lldbpos))
            else
              ldtabpos=ldbpos+(Eldpd(i)-eb)/(ee-eb)*(ldepos-ldbpos)
            endif
            if (ldbneg.gt.1..and.ldeneg.gt.1.) then
              lldbneg=log(ldbneg)
              lldeneg=log(ldeneg)
              ldtabneg=exp(lldbneg+(Eldpd(i)-eb)/(ee-eb)*
     &                 (lldeneg-lldbneg))
            else
              ldtabneg=ldbneg+(Eldpd(i)-eb)/(ee-eb)*(ldeneg-ldbneg)
            endif
            phjpmdposj(i,j)=ldtabpos
            phjpmdnegj(i,j)=ldtabneg
          enddo
        enddo
        do i=1,numdensracap
          phjpmdpos(i)=0.0
          phjpmdneg(i)=0.0
          do j=0,numJph
            phjpmdpos(i)=phjpmdpos(i)+phjpmdposj(i,j)
            phjpmdneg(i)=phjpmdneg(i)+phjpmdnegj(i,j)
          enddo
        enddo
      endif
c
c     ldmodelracap=2: using the total ph-NLD model (pphh)
c
      if (ldmodelracap.eq.2) then
c
c     calculation of total ph-NLD with phmodel from 1 to 2
c
c     begin to prepare the input parameters for the function "phdens2" (phmodel=2, multi-pphh)
c
c     define the particle-hole numbers for proton and neutron in the function "phdens2"
c
        phnumracap=k0
        if (phnumracap.eq.1) then
          ptclp=0
          holep=0
          ptcln=1
          holen=1
        endif
        if (phnumracap.eq.2) then
          ptclp=1
          holep=1
          ptcln=0
          holen=0
        endif
        if (phnumracap.eq.3) then
          ptclp=1
          holep=1
          ptcln=1
          holen=1
        endif
        if (phnumracap.eq.4) then
          ptclp=1
          holep=1
          ptcln=2
          holen=2
        endif
        if (phnumracap.eq.5) then
          ptclp=2
          holep=2
          ptcln=1
          holen=1
        endif
        if (phnumracap.eq.6) then
          ptclp=2
          holep=2
          ptcln=2
          holen=2
        endif
c
c     define the other parameters in the function "phdens2"
c
        sigparleldenp=real(gp(0,0))
        sigparleldenn=real(gn(0,0))
        Ephracap=Etotal
        if (flaggshell) then
          racapdamp=ignatyuk(0,0,Ephracap,0)/alev(0,0)
          sigparleldenp=sigparleldenp*racapdamp
          sigparleldenn=sigparleldenn*racapdamp
        endif
        surfwellE=60.0
        surfwellgo=.false.
c
c     end of prepare the input parameters for the function "phdens2" (phmodel=2, multi-pphh)
c
        do i = 1,numdensracap
          phmdpos(i)=phdens2(0,0,
     &      int(ptclp),int(holep),int(ptcln),int(holen),
     &      sigparleldenp,sigparleldenn,
     &      Eldpd(i),surfwellE,surfwellgo)/2.0
ctest     phmdneg(i)=phmdpos(i)
        end do
c
c     calculation of spin-dependent ph-NLD with phmodel from 1 to 2
c
        do j = 0,numJph
          do i = 1,numdensracap
            phmdposj(i,j)=phmdpos(i)/(4.0*dble(j+1)-2.0)
            phmdnegj(i,j)=phmdposj(i,j)
          end do
        end do
      end if
c
c     ldmodelracap=3: using the spin- parity- dependent NLD model (ldmodel=5!!)
c
      if (ldmodelracap.eq.3) then
c
c     calculation of spin-dependent NLD with ldmodel=5
c
        xprty=-1
        wprty=1
        do j = 0,numJph
          do i = 1,numdensracap
            ldmdposj(i,j)=density(0,0,Eldpd(i),Sldpd(j),int(wprty),
     &        int(0),5)
            ldmdnegj(i,j)=density(0,0,Eldpd(i),Sldpd(j),int(xprty),
     &        int(0),5)
          enddo
        enddo
c
c     calculation of total NLD with ldmodel=5
c
        do i=1,numdensracap
          ldmdpos(i)=0.0
          ldmdneg(i)=0.0
          do j=0,numJph
            ldmdpos(i)=ldmdpos(i)+ldmdposj(i,j)
            ldmdneg(i)=ldmdneg(i)+ldmdnegj(i,j)
          enddo
        enddo
      endif
c
c     store values to chglpos[j] and chglneg[j]
c
cAK   do i=1,numdensracap
cAK     if (ldmodelracap.eq.2) then
cAK       chglpos(i)=dble(phmdpos(i))
cAK       chglneg(i)=dble(phmdneg(i))
cAK    elseif (ldmodelracap.eq.3) then
cAK       chglpos(i)=dble(ldmdpos(i))
cAK       chglneg(i)=dble(ldmdneg(i))
cAK     else
cAK       chglpos(i)=dble(phjpmdpos(i))
cAK       chglneg(i)=dble(phjpmdneg(i))
cAK     endif
cAK   enddo

      do j=0,numJph
        do i=1,numdensracap
          if (ldmodelracap.eq.2) then
            chglposj(i,j)=dble(phmdposj(i,j))
            chglnegj(i,j)=dble(phmdnegj(i,j))
          elseif (ldmodelracap.eq.3) then
            chglposj(i,j)=dble(ldmdposj(i,j))
            chglnegj(i,j)=dble(ldmdnegj(i,j))
          else
            chglposj(i,j)=dble(phjpmdposj(i,j))
            chglnegj(i,j)=dble(phjpmdnegj(i,j))
          endif
        enddo
      enddo
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     set the potential type (JLMB or Wood-Saxon) for racap.f:
c     JLMB and Wood-Saxon are available for neutron and proton;
c     Wood-Saxon is available for d,t,h, and a.
c
      if (k0.eq.1.or.k0.eq.2) then
        if (flagjlm) then
           racopt=3
         else
           racopt=1
         endif
       else
         racopt=1
       endif
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     calculation the final wood saxon potential parameter
c
      vncap2=0.
      rvncap2=0.
      avncap2=0.
      i0=0
      eopt=real(0.000001)
      if (racopt.eq.1) then
        call optical(i0,i0,k0,eopt)
        vncap2=v
        rvncap2=rv
        avncap2=av
      endif
c
c     calculation the final state JLMB potential
c
      do i=1,numjlm
        jlmracap2(i)=0.
      enddo
      if (racopt.eq.3) then
        call mom(0,0,dble(0.),dble(0.000001))
        do i=1,numjlm
          jlmracap2(i)=normjlm(0,0,1)*(potjlm(0,0,i,1)+
     +     normjlm(0,0,3)*(rhojlmn(0,0,i,1)+rhojlmp(0,0,i,1))/
     +     (rhojlmn(0,0,i,1)-rhojlmp(0,0,i,1))*potjlm(0,0,i,3))
        enddo
      endif
c
c     set the levels scheme for direct capture
c     nlevexpracap: number of experimental levels in the final nucleus
c
c     nlevexpracap=nlevmax2(0,0)
      nlevexpracap=Nlast(0,0,0)+1
      ispect=3
      if (ispect.eq.2) nlevexpracap=1
c
c     ispect=1, only available experimental levels are used
c     ispect=2, only theoretical levels deduced from NLDs are used,
c               except the G.S.
c     ispect=3, firstly available experimental levels are used, and then
c               theoretical levels deduced from NLDs are used for the rest energy
c               range. (default)
c
c     set the initial value of spectroscopic factor firstly, then read the experimental spectroscopic factor, if any
c
c     The initial spectroscopic factor for experimental levels is 0.06, which is
c     the mean value of all available experimental data in the database(Nucl. Data Sheet).
c     Therefore, for all nuclei to be calculated, the minimum of global deviation between
c     the experimental compiled spectroscopic factor (read from database) and the estimated
c     one (the mean value 0.06 here) is guaranteed.
c
c     The initial spectroscopic factor for theoretical levels (deduced from NLDs) is 0.5.
c     This estimation minimizes the deviation between the cross section computed by
c     ispect=2 (only theoretical levels deduced from NLDs) and the one computed by
c     ispect=3 (available experimental levels and theoretical levels deduced from NLDs) for all nuclei.
c     Note that compared to ispect=3, only theoretical transition schemes deduced by the spin and
c     parity dependent NLDs are take into account for ispect=2, i.e., the allowed transitions to all
c     discrete experimental known final states are replaced in the energy range from zero to $E_{x}$.
c
      do nex=0,numex
        if (nex.le.nlevexpracap-1) then
          spectfac(0,0,nex)=spectfacexp(0,0,nex)
        else
          spectfac(0,0,nex)=spectfacth(0,0)
        endif
      enddo
c
      if (ispect.ne.2) then
      write(filespec,'("levels/spect",a1,"/z",i3.3)') parsym(k0),
     +  ZZ(0,0,0)
      acp=AA(0,0,0)
      filespec=path(1:lenpath)//filespec
      inquire (file=filespec,exist=lexist)
      if (.not.lexist) return
      open(unit=10,file=filespec,status='old')
  100 read(10,*,end=110) i,ka,jlev
      if (ka.ne.acp) then
        do 120 i=1,jlev+1
          read(10,*,end=110)
  120   continue
        goto 100
      else
        do 130 i=1,jlev
          read(10,'(2x,i3,f8.4,f6.2,2x,i2,f10.5)') j,e,spinf,parity,sf
          do nex=0,Nlast(0,0,0)
            if (abs(edis(0,0,nex)-e).lt.1.d-1.and.jdis(0,0,nex).
     +        eq.spinf.and.parlev(0,0,nex).eq.parity) then
              spectfac(0,0,nex)=sf
              goto 130
            endif
          enddo
  130   continue
      endif
  110 close(10)
      endif
      return
      end

Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
