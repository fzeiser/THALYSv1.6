      subroutine compound(Zcomp,Ncomp,nex,J2,parity)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 17, 2011
c | Task  : Hauser-Feshbach model for multiple emission
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,nex,J2,parity,J,iloop,type,parspin2,
     +                 pspin2,Zix,Nix,NL,nexout,Pprime,Ir,Pprimebeg,
     +                 Pprimeend,Irspin2beg,Irspin2end,J2res,l2maxhf,
     +                 pardif2,Irspin2,J2minI2,J2plusI2,lprimebeg,
     +                 lprimeend,lb,lprime,irad,l2prime,
     +                 jj2primebeg,jj2primeend,jj2prime,updown2
      real             Explus,Exmin,dE1,dE2,fisfeed,s2plus1
      double precision tfd,tf,tfu,fiscontr1,fiscontr2,fiscontr,sumIPE,
     +                 sumIP,rho,total,totalrho,factor,leftover
c
c **************************** Initialization **************************
c
c Zcomp  : charge number index for compound nucleus
c Ncomp  : neutron number index for compound nucleus
c J2     : 2 * J
c parity : parity
c denomhf: denominator for compound nucleus formula
c
      J=J2/2
      denomhf=0.
c
c ********************* Loop over outgoing channels ********************
c
c Note that this whole subroutine is performed inside a loop over
c excitation energy bins nex, compound spin J and parity P.
c The loop over all outgoing particles, energies, residual spins and
c parities is done twice, and is designated by iloop.
c
c iloop=1: Loop to determine the partial and total widths of the
c          Hauser-Feshbach formula.
c
c iloop=2: Loop to perform the actual compound nucleus calculation,
c          using the ratio partial width/total width.
c
      do 10 iloop=1,2
c
c flagfission : flag for fission
c nfisbar     : number of fission barrier parameters
c tfd,tf,tfu  : help variables
c tfisdown    : fission transmission coefficients
c transeps    : absolute limit for transmission coefficient
c tfis,tfisup : fission transmission coefficients
c Explus,Exmin: help variables
c Exmax       : maximum excitation energy for residual nucleus
c Exinc       : excitation energy of entrance bin
c dExinc      : excitation energy bin for mother nucleus
c dE1,dE2     : help variables
c fiscontr    : fission contribution
c fisfeed     : cross section from compound nucleus to fission
c feed        : feeding term for compound nucleus, created during
c               iloop=1
c xsfeed      : cross section from compound to residual nucleus
c fisfeedex   : fission contribution from excitation energy bin
c
c 1. Fission
c
        if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0) then
c
c iloop=1: The fission width for the mother excitation energy bin with
c          spin J and parity P is calculated and added to the total
c          decay width of the mother bin. We use a logarithmic
c          integration of the fission transmission coefficients.
c
          if (iloop.eq.1) then
            tfd=max(tfisdown(J,parity),transeps)
            tf=max(tfis(J,parity),transeps)
            tfu=max(tfisup(J,parity),transeps)
            Explus=min(Exmax(Zcomp,Ncomp),Exinc+0.5*dExinc)
            Exmin=max(Exinc-0.5*dExinc,0.)
            dE1=Exinc-Exmin
            dE2=Explus-Exinc
            if (tfd.eq.tf) then
              fiscontr1=tf*dE1
            else
              fiscontr1=(tf-tfd)/(log(tf)-log(tfd))*dE1
            endif
            if (tfu.eq.tf) then
              fiscontr2=tf*dE2
            else
              fiscontr2=(tfu-tf)/(log(tfu)-log(tf))*dE2
            endif
            fiscontr=(fiscontr1+fiscontr2)/(Explus-Exmin)
            if (fiscontr.le.10.*transeps) fiscontr=0.
            denomhf=fiscontr
          else
c
c iloop=2: The fission contribution for the mother excitation energy bin
c          with spin J and parity P is calculated and added to the
c          various feeding arrays.
c
            fisfeed=0.
            if (denomhf.ne.0.) then
              fisfeed=real(feed*fiscontr)
              xsfeed(Zcomp,Ncomp,-1)=xsfeed(Zcomp,Ncomp,-1)+fisfeed
              fisfeedex(Zcomp,Ncomp,nex)=fisfeedex(Zcomp,Ncomp,nex)+
     +          fisfeed
            endif
          endif
        endif
c
c 20: Photon and particle channels
c
c parskip     : logical to skip outgoing particle
c parspin2    : 2 * particle spin
c parspin     : spin of particle
c pspin2,spin2: 2 * spin of particle
c Zindex,Zix  : charge number index for residual nucleus
c Nindex,Nix  : neutron number index for residual nucleus
c Nlast,NL    : last discrete level
c
        if (iloop.eq.1) then
          do 15 Pprime=-1,1,2
            do 15 Ir=0,numJ
              do 15 type=0,6
                do 15 nexout=0,nexmax(type)
                  enumhf(type,nexout,Ir,Pprime)=0.
   15     continue
        endif
        do 20 type=0,6
          if (parskip(type)) goto 20
          if (iloop.eq.1.and.type.eq.6.and.denomhf.eq.0.) goto 20
          parspin2=int(2.*parspin(type))
          s2plus1=parspin2+1.
          pspin2=spin2(type)
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          NL=Nlast(Zix,Nix,0)
c
c 30: Sum over outgoing excitation energies
c
c sumIPE: compound contribution summed over residual spin and parity
c         and energy
c nexmax: maximum excitation energy bin for residual nucleus
c numJ  : maximal J-value
c enumhf: enumerator for compound nucleus formula
c
c This loop is over all discrete levels and continuum bins of the
c final nucleus.
c
          sumIPE=0.
          do 30 nexout=0,nexmax(type)
            if (nexout.eq.0.and.NL.eq.0) goto 30
c
c Initialization of summations
c
c Pprimebeg : start of residual parity summation
c parlev    : parity of level
c Pprimeend : end of residual parity summation
c Irspin2beg: 2 * start of residual spin summation
c jdis      : spin of level
c Irspin2end: 2 * end of residual spin summation
c J2res     : help variable
c maxJ      : maximal J-value
c l2maxhf   : 2 * lmaxhf
c lmaxhf    : maximal l-value for transmission coefficients
c sumIP     : compound contribution summed over residual spin and parity
c
c For discrete states, the begin and end points of the residual
c spin/parity summation are both set equal to the residual discrete
c level spin/parity.
c
            if (nexout.le.NL) then
              Pprimebeg=parlev(Zix,Nix,nexout)
              Pprimeend=Pprimebeg
              Irspin2beg=int(2.*jdis(Zix,Nix,nexout))
              Irspin2end=Irspin2beg
            else
c
c For the continuum, the begin and end points of the residual
c spin/parity summation are set to the maximally accessible values.
c
              Pprimebeg=-1
              Pprimeend=1
              J2res=J2+parspin2
              Irspin2beg=mod(J2res,2)
              Irspin2end=2*maxJ(Zix,Nix,nexout)
            endif
            l2maxhf=2*lmaxhf(type,nexout)
            if (iloop.eq.1) then
c
c 50: Sum over residual parity
c
c pardif2: difference between residual and compound nucleus parity
c
c The variable pardif2 is used as an indicator of parity conservation
c for the outgoing channel.
c
              do 50 Pprime=Pprimebeg,Pprimeend,2
                pardif2=abs(parity-Pprime)/2
c
c 60: Sum over residual spin
c
c Irspin2 : 2 * residual spin
c Ir      : residual spin
c rho,rho0: integrated level density
c
                do 60 Irspin2=Irspin2beg,Irspin2end,2
                  Ir=Irspin2/2
                  rho=rho0(type,nexout,Ir,Pprime)
c
c The Hauser-Feshbach formula contains the following triangular
c relations:
c |J-I| < j < J+I
c |j-s| < l < j+s
c Computationally it is more economic to do the l-summation first.
c This can be done with the following equivalent triangular relations:
c ||J-I|-s| < l < |J+I+s|
c max(|J-I|,|l-s|) < j < min((J+I),l+s)
c
c 70,80: Sum over l of outgoing channel
c
c total      : help variable
c J2minI2,...: help variables
c lprimebeg  : start of l summation
c lprimeend  : end of l summation
c
                  total=0.
                  J2minI2=abs(J2-Irspin2)
                  J2plusI2=J2+Irspin2
                  lprimebeg=abs(J2minI2-parspin2)/2
                  lprimeend=abs(J2plusI2+parspin2)/2
                  lprimeend=min(lprimeend,l2maxhf/2)
c
c We include photons as a special case, with the multipole radiation
c selection rules (irad=0: M-transition, irad=1: E-transition)
c
c lprime: l'
c irad  : variable to indicate M(=0) or E(=1) radiation
c Tgam  : gamma transmission coefficients
c
c 1. Photons
c
                  if (type.eq.0) then
                    do 70 lprime=lprimebeg,lprimeend
                      irad=1-abs(mod(lprime,2)-pardif2)
                      total=total+Tgam(nexout,lprime,irad)
   70               continue
                  else
c
c 2. Particles
c
c 70: Sum over j (jj2) of outgoing channel
c
c flagfullhf : flag for full spin dependence of transmission
c              coefficients
c jj2primebeg: 2 * start of j' summation
c jj2primeend: 2 * end of j' summation
c jj2prime   : 2 * j'
c l2prime    : 2 * l'
c updown2    : spin index for transmission coefficient
c Tjlnex     : transmission coefficients as a function of
c              particle type, energy, spin and l-value
c Tlnex      : transmission coefficients as a function of
c              particle type, energy and l-value (averaged over spin)
c
c If the parity of the residual nucleus is equal (unequal) to the parity
c of compound nucleus,i.e. pardif2=0(1), the l-value must be even (odd).
c
c A choice is amde between averaged and full spin dependence of
c the transmission coefficients in the Hauser-Feshbach model (only
c significant in case of very large spin-orbit terms in the OMP).
c
                    lb=lprimebeg+abs(mod(lprimebeg,2)-pardif2)
                    if (flagfullhf) then
                      do 80 lprime=lb,lprimeend,2
                        l2prime=2*lprime
                        jj2primebeg=max(J2minI2,abs(l2prime-parspin2))
                        jj2primeend=min(J2plusI2,l2prime+parspin2)
                        do 90 jj2prime=jj2primebeg,jj2primeend,2
                          updown2=(jj2prime-l2prime)/pspin2
                          total=total+Tjlnex(type,nexout,updown2,lprime)
   90                   continue
   80                 continue
                    else
                      do 100 lprime=lb,lprimeend,2
                        total=total+Tlnex(type,nexout,lprime)
  100                 continue
                      total=s2plus1*total
                    endif
                  endif
c
c iloop=1: The partial decay widths are stored in enumhf and also
c          added to the total decay width denomhf.
c
c totalrho: help variable
c
                  totalrho=rho*total
                  enumhf(type,nexout,Ir,Pprime)=totalrho
                  denomhf=denomhf+totalrho
   60           continue
   50         continue
            else
c
c ** Populate the outgoing channels using the compound nucleus formula *
c
c factor    : help variable
c xspop     : population cross section
c xspopex   : population cross section summed over spin and parity
c mcontrib  : contribution to emission spectrum
c xspopnuc  : population cross section per nucleus
c xspartial : emitted cross section flux per energy bin
c xsfeed    : cross section from compound to residual nucleus
c
c iloop=2: Increment of the the population arrays of the residual
c          nuclei and the mcontrib array, which will be used for the
c          interpolation of the compound emission spectra.
c
              sumIP=0.
              do 110 Pprime=Pprimebeg,Pprimeend,2
                do 120 Irspin2=Irspin2beg,Irspin2end,2
                  Ir=Irspin2/2
                  factor=real(feed*enumhf(type,nexout,Ir,Pprime))
                  xspop(Zix,Nix,nexout,Ir,Pprime)=
     +              xspop(Zix,Nix,nexout,Ir,Pprime)+factor
                  sumIP=sumIP+factor
  120           continue
  110         continue
              xspopex(Zix,Nix,nexout)=xspopex(Zix,Nix,nexout)+sumIP
              xspopex(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)-sumIP
              mcontrib(type,nex,nexout)=mcontrib(type,nex,nexout)+sumIP
              sumIPE=sumIPE+sumIP
            endif
   30     continue
          if (iloop.eq.2) then
            xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)+sumIPE
            xspartial(type,nex)=xspartial(type,nex)+sumIPE
            xsfeed(Zcomp,Ncomp,type)=xsfeed(Zcomp,Ncomp,type)+sumIPE
          endif
   20   continue
        if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0.and.iloop.eq.2)
     +    then
          xspopex(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)-fisfeed
        endif
c
c ** Create feeding term for compound nucleus decay in the second loop *
c
c Dmulti  : depletion factor for multiple preequilibrium
c leftover: remaining cross section flux after decay (should be close
c           to zero)
c
        if (iloop.eq.1) then
          if (denomhf.ne.0.) then
            feed=(1.-Dmulti(nex))*xspop(Zcomp,Ncomp,nex,J,parity)/
     +        denomhf
          else
c
c Prevent trapping of cross section in the continuum. Any flux that
c is left in the compound nucleus bin is equally distributed over the
c discrete states.
c
            NL=Nlast(Zcomp,Ncomp,0)
            feed=0.
            leftover=xspop(Zcomp,Ncomp,nex,J,parity)/(NL+1.)
            xspartial(0,nex)=xspartial(0,nex)+
     +        xspop(Zcomp,Ncomp,nex,J,parity)
            do 210 nexout=0,NL
              xspopex(Zcomp,Ncomp,nexout)=xspopex(Zcomp,Ncomp,nexout)+
     +          leftover
              Ir=int(jdis(Zcomp,Ncomp,nexout))
              Pprime=parlev(Zcomp,Ncomp,nexout)
              xspop(Zcomp,Ncomp,nexout,Ir,Pprime)=
     +          xspop(Zcomp,Ncomp,nexout,Ir,Pprime)+leftover
              mcontrib(0,nex,nexout)=mcontrib(0,nex,nexout)+leftover
  210       continue
          endif
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
