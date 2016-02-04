      subroutine partfunc
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire, Stephane Goriely, Arjan Koning
c | Date  : December 30, 2011
c | Task  : Calculate partition function
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          numdiv
      parameter (numdiv=100)
      integer          Zix,Nix,NL,ldmod,i,nex,nexout,Pprime,Ir,idiv
      real             Elev,spindeg,dex,MeVkT,Rspin
      double precision sumjp(0:numex*numdiv),sum,fex,rho,density
c
c ************************ Partition Function **************************
c
c MeVkT   : E/kT expressed in MeV/T9
c Zix     : charge number index for target nucleus
c parZ    : charge number of particle
c k0      : index of incident particle
c Nix     : neutron number index for target nucleus
c parN    : neutron number of particle
c Nlast,NL: last discrete level
c deltaEx : excitation energy bin for population arrays
c numdiv  : number of subdivisions considered for each continuum bin
c           in target nucleus
c nTmax   : effective number of temperatures for Maxwellian
c numex   : maximal number of excitation energies
c sumjp   : help variable
c maxex   : maximum excitation energy bin for residual nucleus
c nlev    : number of excited levels for nucleus
c
      MeVkT=11.605
      Zix=parZ(k0)
      Nix=parN(k0)
      NL=Nlast(Zix,Nix,0)
      ldmod=ldmodel(Zix,Nix)
      do 10 i=1,nTmax
        do 20 nex=0,numdiv*numex
          sumjp(nex)=0.
   20   continue
        sum=0.
        nex=0
        do 30 nexout=0,max(maxex(Zix,Nix)-1,nlev(Zix,Nix))
c
c For discrete states, the begin and end points of the target
c spin/parity summation are both set equal to the target discrete
c level spin/parity.
c
c Elev,Edis: excitation energy in the target nucleus
c spindeg  : (2J+1) degeneracy of target level
c jdis     : spin of level
c fex,sum  : help variables
c
          dex=deltaEx(Zix,Nix,nexout)/real(numdiv)
          Elev=edis(Zix,Nix,nexout)
          if (nexout.le.NL) then
            spindeg=2.*jdis(Zix,Nix,nexout)+1
            fex=MeVkT*Elev/T9(i)
            if (fex.gt.80.) goto 30
            sum=sum+spindeg*exp(-fex)
            sumjp(0)=spindeg*exp(-fex)
          else
c
c For decay to the continuum we use a spin and parity dependent level
c density.
c See subroutine exgrid for change in case of non-equiprobable parities.
c
c Pprime     : parity variable in target continuum bin
c Ir         : level spin variable in target continuum bin
c Rspin      : true level spin in target continuum bin
c numJ       : maximal J-value
c rho        : level density
c Atarget    : mass number of target nucleus
c partf      : integrated partition function
c targetspin2: 2 * spin of target
c
            do 40 idiv=1,numdiv
              Elev=Ex(Zix,Nix,nexout)+(real(idiv)-0.5*numdiv)*dex
              fex=MeVkT*Elev/T9(i)
              if (fex.gt.80.) goto 40
              nex=nex+1
              do 50 Pprime=-1,1,2
                do 60 Ir=0,numJ
                  Rspin=Ir+mod(Atarget,2)/2.
                  spindeg=2.*(Ir+mod(Atarget,2)/2.)+1.
                  rho=density(Zix,Nix,Elev,Rspin,Pprime,0,ldmod)
                  sumjp(nex)=sumjp(nex)+spindeg*rho*exp(-fex)
   60           continue
   50         continue
   40       continue
          endif
   30   continue
        do 70 idiv=1,nex
          sum=sum+(sumjp(idiv-1)+sumjp(idiv))/2.*dex
   70   continue
        partf(i)=sum/(targetspin2+1.)
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
