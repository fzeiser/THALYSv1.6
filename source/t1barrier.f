      subroutine t1barrier(Zcomp,Ncomp,J2,parity,ibar,trfis,rhof,Eex,
     +  iloop)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Marieke Duijvestijn
c | Date  : December 15, 2013
c | Task  : Fission transmission coefficient for one barrier
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80     key
      integer          Zcomp,Ncomp,J2,parity,ibar,iloop,J,itr,j2trans,
     +                 pitrans,ihill,i
      real             Eex,fbar,factor1,factor2,bfis,wfis,etrans,Eeff,
     +                 twkbint,thill,elow,emid,eup,dE1,dE2
      double precision trfis,rhof,trfisone,rho1,rho2,rho3,r1log,r2log,
     +                 r3log,rho,rhotr
c
c ********** Fission transmission coefficient for one barrier **********
c
c Zcomp : charge number index for compound nucleus
c Ncomp : neutron number index for compound nucleus
c J2    : 2 * J
c parity: parity
c ibar  : fission barrier
c trfis : transmission coefficient
c rhof  : integrated level density
c Eex   : excitation energy of entrance bin
c
      J=J2/2
      trfis=0.
      rhof=0.
c
c Correct LDM barrier height with ground state shell correction
c
c fisadjust : logical for energy-dependent fission adjustment
c adjust    : subroutine for energy-dependent parameter adjustment
c factor1,. : multiplication factor
c fbaradjust: adjustable factors for fission parameters
c             (default 1.)
c fbarrier  : height of fission barrier
c fismodelx : fission model
c nfisbar   : number of fission barrier parameters
c bfis,wfis : help variables
c deltaW    : shell correction in nuclear mass
c fwidth    : width of fission barrier
c
      if (fisadjust(Zcomp,Ncomp)) then
        key='fisbar'
        call adjust(Eex,key,Zcomp,Ncomp,0,ibar,factor1)
        key='fisbaradjust'
        call adjust(Eex,key,Zcomp,Ncomp,0,ibar,factor2)
        fbar=factor1*factor2*fbarrier(Zcomp,Ncomp,ibar)*
     +    fbaradjust(Zcomp,Ncomp,ibar)
      else
        fbar=fbarrier(Zcomp,Ncomp,ibar)
      endif
      if ((fismodelx(Zcomp,Ncomp).ge.3).or.
     +  ((fismodelx(Zcomp,Ncomp).lt.3).and.(nfisbar(Zcomp,Ncomp).eq.1)))
     +  then
        bfis=fbar-(deltaW(Zcomp,Ncomp,0)-deltaW(Zcomp,Ncomp,1))
      else
        bfis=fbar
      endif
      if (fisadjust(Zcomp,Ncomp)) then
        key='fishw'
        call adjust(Eex,key,Zcomp,Ncomp,0,ibar,factor1)
        key='fishwadjust'
        call adjust(Eex,key,Zcomp,Ncomp,0,ibar,factor2)
        wfis=factor1*factor2*fwidth(Zcomp,Ncomp,ibar)*
     +    fwidthadjust(Zcomp,Ncomp,ibar)
      else
        wfis=fwidth(Zcomp,Ncomp,ibar)
      endif
c
c 1. Discrete states
c
c nfistrrot     : number of rotational transition states for barrier
c etrans,j2trans: help variables
c efistrrot     : energy of rotational transition states
c jfistrrot     : spin of rotational transition states
c pitrans,Eeff  : help variables
c pfistrrot     : parity of rotational transition states
c thill         : Hill-Wheeler penetrability
c primary       : flag to designate primary (binary) reaction
c trfisone      : help variable
c twkbint       : WKB penetrability
c ihill         : counter for Hill-Wheeler magnitude
c numhill       : maximum number of Hill-Wheeler points
c tfisA         : transmission coefficient for Hill-Wheeler magnitude
c rhofisA       : integrated level density corresponding to tfisA
c
      do 10 itr=1,nfistrrot(Zcomp,Ncomp,ibar)
        etrans=efistrrot(Zcomp,Ncomp,ibar,itr)
        if (Eex.lt.etrans) goto 10
        j2trans=int(2.*jfistrrot(Zcomp,Ncomp,ibar,itr))
        pitrans=pfistrrot(Zcomp,Ncomp,ibar,itr)
        Eeff=Eex-etrans
        if ((J2.eq.j2trans).and.(parity.eq.pitrans)) then
          if (fismodelx(Zcomp,Ncomp).eq.5) then
            trfisone=twkbint(Eeff,ibar,Zcomp,Ncomp)
          else
            trfisone=thill(Eeff,bfis,wfis)
          endif
          if ((ibar.eq.1).and.primary.and.iloop.eq.2) then
            ihill=min(int(numhill*trfisone)+1,numhill)
            tfisA(J,parity,ihill)=tfisA(J,parity,ihill)+trfisone
            tfisA(J,parity,0)=tfisA(J,parity,0)+trfisone
            rhofisA(J,parity,ihill)=rhofisA(J,parity,ihill)+1.
          endif
          trfis=trfis+trfisone
          rhof=rhof+1.
        endif
   10 continue
c
c 2. Continuum
c
c fecont       : start of continuum energy
c nbintfis     : number of integration bins
c elow,eup,emid: help variables
c eintfis      : excitation energy for fission
c dE1,dE2      : help variables
c rho1-3       : help variables
c r1log,...    : help variables
c rhofis       : integrated level density
c rho          : integrated level density
c
      if (Eex.ge.fecont(Zcomp,Ncomp,ibar)) then
        do 20 i=1,nbintfis(ibar)-2,2
          elow=eintfis(i,ibar)
          if (elow.gt.Eex) goto 20
          emid=min(eintfis(i+1,ibar),Eex)
          eup=min(eintfis(i+2,ibar),Eex)
          dE1=emid-elow
          dE2=eup-emid
          rho1=rhofis(i,J,parity,ibar)*(1.+1.d-10)+1.d-10
          rho2=rhofis(i+1,J,parity,ibar)+1.d-10
          rho3=rhofis(i+2,J,parity,ibar)*(1.+1.d-10)+1.d-10
          r1log=log(rho1)
          r2log=log(rho2)
          r3log=log(rho3)
          if (r2log.ne.r1log.and.r2log.ne.r3log) then
            rho=(rho1-rho2)/(r1log-r2log)*dE1
     +        +(rho2-rho3)/(r2log-r3log)*dE2
          else
            rho=rho2*(dE1+dE2)
          endif
          Eeff=Eex-emid
          if (fismodelx(Zcomp,Ncomp).eq.5) then
            trfisone=twkbint(Eeff,ibar,Zcomp,Ncomp)
          else
            trfisone=thill(Eeff,bfis,wfis)
          endif
          rhotr=rho*trfisone
          trfis=trfis+rhotr
          rhof=rhof+rho
          if ((ibar.eq.1).and.primary.and.iloop.eq.2) then
            ihill=min(int(numhill*trfisone)+1,numhill)
            tfisA(J,parity,ihill)=tfisA(J,parity,ihill)+rhotr
            tfisA(J,parity,0)=tfisA(J,parity,0)+rhotr
            rhofisA(J,parity,ihill)=rhofisA(J,parity,ihill)+rho
          endif
   20   continue
        if (iloop.eq.2)
     +    tfisA(J,parity,0)=max(tfisA(J,parity,0),1.d-30)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
