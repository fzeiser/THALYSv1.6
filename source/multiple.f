      subroutine multiple
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : December 9, 2013
c | Task  : Multiple emission
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*13 fisfile
      character*60 form1,form2
      integer      Zcomp,Ncomp,type,nen,nex,nexout,Zix,Nix,Z,N,A,odd,NL,
     +             J,idensfis,parity,J2,Jfis,iang,p,h
      real         popepsA,xspopsave(0:numex),Exm,dEx,Exmin,popepsB,
     +             xsmax,Smax,emissum(0:numpar),Eaveragesum,Eout,sumxs,
     +             xsdif,ang,kalbach,Eaverage(0:numpar)
c
c ******************** Loop over nuclei ********************************
c
c Loop over all residual nuclei, starting with the initial compound
c nucleus (Zcomp=0, Ncomp=0), and then according to decreasing Z and N.
c
c primary    : flag to designate primary (binary) reaction
c flagomponly: flag to execute ONLY an optical model calculation
c flaginitpop: flag for initial population distribution
c excitation : subroutine for excitation energy population
c flagpop    : flag for output of population
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c maxZ       : maximal number of protons away from the initial compound
c              nucleus
c maxN       : maximal number of neutrons away from the initial compound
c              nucleus
c
      primary=.false.
      if (flagomponly) return
      if (flaginitpop) call excitation
      if (flagpop)
     +  write(*,'(/" ########## MULTIPLE EMISSION ##########")')
      do 10 Zcomp=0,maxZ
        do 10 Ncomp=0,maxN
c
c Continue for this compound nucleus only when it contains sufficient
c reaction flux. If so, determine the nuclear structure for the residual
c nuclei. Initialize emission spectrum for this compound nucleus.
c Determine excitation energy grid for residual nuclei.
c
c xspopnuc  : population cross section per nucleus
c popeps    : limit for population cross section per nucleus
c parskip   : logical to skip outgoing particle
c flagspec  : flag for output of spectra
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c xsemis    : emission spectrum from compound nucleus
c xsmpeemis : multiple-preequilibrium emission spectrum from compound
c             nucleus
c numex     : maximal number of excitation energies (set in talys.cmb)
c xspartial : emitted cross section flux per energy bin
c xsmpe     : multiple-preequilibrium cross section per energy bin
c mcontrib  : contribution to emission spectrum
c mpecontrib: contribution to multiple pre-equilibrium emission spectrum
c xsmpetot  : total multiple-preequilibrium cross section
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c strucexist: flag to determine whether structure info for this nucleus
c             already exists
c structure : subroutine for nuclear structure parameters
c xspopsave : help variable for diagnosis
c Dmulti    : depletion factor for multiple preequilibrium
c exgrid    : subroutine to set excitation energy grid
c
          if (xspopnuc(Zcomp,Ncomp).lt.popeps) then
            xspopnuc(Zcomp,Ncomp)=0.
            goto 500
          endif
          do 20 type=0,6
            if (flagspec) then
              do 30 nen=0,numen
                xsemis(type,nen)=0.
                xsmpeemis(type,nen)=0.
   30         continue
            endif
            do 40 nex=0,numex+1
              xspartial(type,nex)=0.
              xsmpe(type,nex)=0.
              do 40 nexout=0,numex+1
                mcontrib(type,nex,nexout)=0.
                mpecontrib(type,nex,nexout)=0.
   40       continue
            xsmpetot(type)=0.
            if (parskip(type)) goto 20
            Zix=Zindex(Zcomp,Ncomp,type)
            Nix=Nindex(Zcomp,Ncomp,type)
            if (.not.strucexist(Zix,Nix)) then
              call structure(Zix,Nix)
              strucexist(Zix,Nix)=.true.
            endif
   20     continue
          do 50 nex=0,numex
            xspopsave(nex)=0.
            Dmulti(nex)=0.
   50     continue
          call exgrid(Zcomp,Ncomp)
c
c If a new set of transmission coefficients for each residual
c nuclide is requested, we calculate inverse cross sections and
c transmission coefficients for this compound nucleus here.
c
c flagompall: flag for new optical model calculation for all residual
c             nuclei
c basicxs   : subroutine for basic cross sections and transmission
c             coefficients
c
          if (flagompall) call basicxs(Zcomp,Ncomp)
c
c Write population of compound nucleus before its decay
c
c dExinc       : excitation energy bin for mother nucleus
c deltaEx      : excitation energy bin for population arrays
c ZZ,Z         : charge number of residual nucleus
c NN,N         : neutron number of residual nucleus
c AA,A         : mass number of residual nucleus
c odd          : odd (1) or even (0) nucleus
c strucwrite   : flag for output of nuclear structure info
c flaglevels   : flag for output of discrete level information
c levelsout    : subroutine for output of discrete levels
c flagdensity  : flag for output of level densities
c densityout   : subroutine for output of level density parameters
c flagfisout   : flag for output of fission information
c fissionparout: subroutine for output for fission parameters
c nuc          : symbol of nucleus
c Exmax        : maximum excitation energy for residual nucleus
c Nlast,NL     : last discrete level
c maxex        : maximum excitation energy bin for compound nucleus
c Ex           : excitation energy
c xspopex      : population cross section summed over spin and parity
c xspop        : population cross section
c
          dExinc=deltaEx(Zcomp,Ncomp,maxex(Zcomp,Ncomp))
          Z=ZZ(Zcomp,Ncomp,0)
          N=NN(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          odd=mod(A,2)
          if (flagpop) then
            if (.not.strucwrite(Zcomp,Ncomp)) then
              if (flaglevels) call levelsout(Zcomp,Ncomp)
              if (flagdensity) call densityout(Zcomp,Ncomp)
              if (flagfisout) call fissionparout(Zcomp,Ncomp)
              strucwrite(Zcomp,Ncomp)=.true.
            endif
            write(*,'(/" Population of Z=",i3," N=",i3," (",i3,a2,
     +        ") before decay:",1p,e12.5)') Z,N,A,nuc(Z),
     +        xspopnuc(Zcomp,Ncomp)
            NL=Nlast(Zcomp,Ncomp,0)
            if (maxex(Zcomp,Ncomp).gt.NL) then
              write(*,'(" Maximum excitation energy:",f8.3,
     +          " Discrete levels:",i3," Continuum bins:",i3,
     +          " Continuum bin size:",f8.3/)') Exmax(Zcomp,Ncomp),NL,
     +          maxex(Zcomp,Ncomp)-NL,dExinc
            else
              write(*,'(" Maximum excitation energy:",f8.3,
     +          " Discrete levels:",i3)') Exmax(Zcomp,Ncomp),NL
            endif
            write(*,'(" bin    Ex    Popul.",10("    J=",f4.1)/)')
     +        (J+0.5*odd,J=0,9)
            do 60 nex=0,maxex(Zcomp,Ncomp)
              if (nex.le.NL) then
                write(*,'(1x,i3,f8.3,1p,e10.3,1p,10e10.3)') nex,
     +            Ex(Zcomp,Ncomp,nex),xspopex(Zcomp,Ncomp,nex),
     +            (max(xspop(Zcomp,Ncomp,nex,J,1),
     +            xspop(Zcomp,Ncomp,nex,J,-1)),J=0,9)
              else
                write(*,'(1x,i3,f8.3,1p,e10.3,1p,10e10.3)') nex,
     +            Ex(Zcomp,Ncomp,nex),xspopex(Zcomp,Ncomp,nex),
     +            (xspop(Zcomp,Ncomp,nex,J,1),J=0,9)
              endif
   60       continue
          endif
c
c Loop 110: De-excitation of nucleus, starting at the highest excitation
c energy bin.
c
c popepsA : limit for population cross sections per energy
c idensfis: identifier to activate possible fission level densities
c
c Continue for this (Zcomp,Ncomp,nex) only when there is sufficient
c reaction flux in the excitation energy bin.
c The fission transmission coefficients and level densities only need
c to be calculated once, at the highest excitation energy.
c
          popepsA=popeps/max(5*maxex(Zcomp,Ncomp),1)
          idensfis=1
          if (flagpop) then
            write(*,'(/" Population of each bin of Z=",i3," N=",i3,
     +      " (",i3,a2,") before its decay"/)') Z,N,A,nuc(Z)
            write(*,'(" bin    Ex    Popul.",10("    J=",f4.1)/)')
     +        (J+0.5*odd,J=0,9)
          endif
          do 110 nex=maxex(Zcomp,Ncomp),1,-1
            dExinc=deltaEx(Zcomp,Ncomp,nex)
            if (flagpop) then
              xspopsave(nex)=xspopex(Zcomp,Ncomp,nex)
              if (nex.le.NL) then
                write(*,'(1x,i3,f8.3,1p,e10.3,1p,10e10.3)') nex,
     +            Ex(Zcomp,Ncomp,nex),xspopex(Zcomp,Ncomp,nex),
     +            (max(xspop(Zcomp,Ncomp,nex,J,1),
     +            xspop(Zcomp,Ncomp,nex,J,-1)),J=0,9)
              else
                write(*,'(1x,i3,f8.3,1p,e10.3,1p,10e10.3)') nex,
     +            Ex(Zcomp,Ncomp,nex),xspopex(Zcomp,Ncomp,nex),
     +            (xspop(Zcomp,Ncomp,nex,J,1),J=0,9)
              endif
            endif
c
c For exclusive channel cross section calculations, some variables
c need to be stored in extra arrays.
c
c flagchannels: flag for exclusive channels calculation
c popexcl     : population cross section of bin just before decay
c
            if (flagchannels) popexcl(Zcomp,Ncomp,nex)=
     +        xspopex(Zcomp,Ncomp,nex)
c
c Discrete levels decay by gamma cascade. Isomers are excluded
c from gamma cascade.
c Note that we assume that discrete levels cannot particle decay.
c This needs to be added (for light nuclei) in a future version.
c
c tau    : lifetime of state in seconds
c cascade: subroutine for gamma-ray cascade
c Exinc  : excitation energy of entrance bin
c
            if (nex.le.Nlast(Zcomp,Ncomp,0)) then
              if (tau(Zcomp,Ncomp,nex).eq.0.)
     +          call cascade(Zcomp,Ncomp,nex)
              goto 110
            endif
            if (xspopex(Zcomp,Ncomp,nex).lt.popepsA) goto 110
            Exinc=Ex(Zcomp,Ncomp,nex)
c
c For each mother excitation energy bin, determine the highest
c possible excitation energy bin nexmax for the residual nuclei.
c As reference, we take the top of the mother bin. The maximal
c residual excitation energy Exm is obtained by subtracting the
c separation energy from this. The bin in which this Exm falls
c is denoted by nexmax.
c
c Exm   : maximal attainable energy
c S     : separation energy per particle
c egrid : outgoing energy grid
c dEx   : excitation energy bin for population arrays
c Exmin : lower bound of energy bin
c nexmax: maximum excitation energy bin for residual nucleus
c
            do 120 type=1,6
              if (parskip(type)) goto 120
              Exm=Exinc+0.5*dExinc-S(Zcomp,Ncomp,type)
              if (type.gt.1) Exm=Exm-egrid(ebegin(type))
              Zix=Zindex(Zcomp,Ncomp,type)
              Nix=Nindex(Zcomp,Ncomp,type)
              do 130 nexout=maxex(Zix,Nix),0,-1
                dEx=deltaEx(Zix,Nix,nexout)
                Exmin=Ex(Zix,Nix,nexout)-0.5*dEx
                if (Exmin.lt.Exm) then
                  nexmax(type)=nexout
                  goto 120
                endif
  130         continue
              nexmax(type)=-1
  120       continue
            nexmax(0)=nex-1
c
c Prepare transmission coefficients and level densities. This can
c be done outside the loops over spin and parity.
c
c densprepare: subroutine to prepare energy grid and level density
c              information for compound nucleus
c
            call densprepare(Zcomp,Ncomp,idensfis)
            idensfis=0
c
c Deplete excitation energy bin by multiple pre-equilibrium.
c
c mulpreZN   : logical for multiple pre-equilibrium per nucleus
c flag2comp  : flag for two-component pre-equilibrium model
c multipreeq2: subroutine for two-component multiple preequilibrium
c multipreeq : subroutine for multiple preequilibrium
c
            if (mulpreZN(Zcomp,Ncomp)) then
              if (flag2comp) then
                call multipreeq2(Zcomp,Ncomp,nex)
              else
                call multipreeq(Zcomp,Ncomp,nex)
              endif
            endif
c
c Compound nucleus decay of mother excitation energy/spin/parity bin.
c
c flagcomp: flag for compound nucleus calculation
c popepsB : limit for population cross sections per spin and parity
c maxJ    : maximal J-value
c parity  : parity
c
            if (flagcomp) then
              popepsB=popepsA/(5*maxJ(Zcomp,Ncomp,nex))*0.5
              do 140 parity=-1,1,2
                do 140 J=0,maxJ(Zcomp,Ncomp,nex)
c
c Continue for this (Zcomp,Ncomp,nex,J,P) only when there is sufficient
c reaction flux in the excitation energy bin. The correct value for J
c is determined.
c
c J2         : 2 * J
c flagfission: flag for fission
c nfisbar    : number of fission barrier parameters
c tfission   : subroutine for fission transmission coefficients
c compound   : subroutine for Hauser-Feshbach model for multiple
c              emission
c tfissionout: subroutine for output of fission transmission
c              coefficients
c
                  if (xspop(Zcomp,Ncomp,nex,J,parity).lt.popepsB)
     +              goto 140
                  J2=2*J+odd
                  if (flagfission.and.nfisbar(Zcomp,Ncomp).ne.0)
     +              call tfission(Zcomp,Ncomp,nex,J2,parity)
                  call compound(Zcomp,Ncomp,nex,J2,parity)
  140         continue
              if (flagfisout) call tfissionout(Zcomp,Ncomp,nex)
            endif
  110     continue
c
c Make new population cross section per nucleus
c
c flagcompo : flag for output of cross section components
c xspoppreeq: preequilibrium population cross section per nucleus
c preeqpopex: pre-equilibrium population cross section summed over
c             spin and parity
c xspopcomp : compound population cross section per nucleus
c xspopdir  : direct population cross section per nucleus
c Fdir      : direct population fraction per nucleus
c Fpreeq    : preequilibrium population fraction per nucleus
c Fcomp     : compound population fraction per nucleus
c
          xspopnuc(Zcomp,Ncomp)=xspopex(Zcomp,Ncomp,0)
          do 210 nex=1,Nlast(Zcomp,Ncomp,0)
            if (tau(Zcomp,Ncomp,nex).ne.0.) xspopnuc(Zcomp,Ncomp)=
     +        xspopnuc(Zcomp,Ncomp)+xspopex(Zcomp,Ncomp,nex)
  210     continue
          if (flagcompo) then
            xsmax=0.
            Smax=0.
            do 220 type=1,6
              if (xsfeed(Zcomp,Ncomp,type).gt.xsmax) then
                xsmax=xsfeed(Zcomp,Ncomp,type)
                Smax=S(Zcomp,Ncomp,type)
              endif
  220       continue
            do 230 nex=0,maxex(Zcomp,Ncomp)
              if (Ex(Zcomp,Ncomp,nex).gt.Smax) goto 240
              xspoppreeq(Zcomp,Ncomp)=xspoppreeq(Zcomp,Ncomp)+
     +            preeqpopex(Zcomp,Ncomp,nex)
  230       continue
  240       xspoppreeq(Zcomp,Ncomp)=min(dble(xspoppreeq(Zcomp,Ncomp)),
     +        xspopnuc(Zcomp,Ncomp)-dble(xspopdir(Zcomp,Ncomp)))
            xspopcomp(Zcomp,Ncomp)=max(xspopnuc(Zcomp,Ncomp)-
     +        dble(xspoppreeq(Zcomp,Ncomp)-xspopdir(Zcomp,Ncomp)),0.d0)
            if (xspopnuc(Zcomp,Ncomp).gt.0.) then
              Fdir(Zcomp,Ncomp)=xspopdir(Zcomp,Ncomp)/
     +          xspopnuc(Zcomp,Ncomp)
              Fpreeq(Zcomp,Ncomp)=xspoppreeq(Zcomp,Ncomp)/
     +          xspopnuc(Zcomp,Ncomp)
              Fcomp(Zcomp,Ncomp)=xspopcomp(Zcomp,Ncomp)/
     +          xspopnuc(Zcomp,Ncomp)
            endif
          endif
c
c For exclusive channel calculation: Determine feeding terms.
c
c numZchan: maximal number of outgoing proton units in individual
c           channel description
c numNchan: maximal number of outgoing neutron units in individual
c           channel description
c feedexcl: feeding terms from compound excitation energy bin to
c           residual excitation energy bin
c
          if (flagchannels.and.Zcomp.le.numZchan.and.
     +      Ncomp.le.numNchan) then
            do 250 nex=maxex(Zcomp,Ncomp),1,-1
              do 260 type=0,6
                if (parskip(type)) goto 260
                if (Zcomp.eq.0.and.Ncomp.eq.0.and.type.gt.0.and.
     +            .not.flaginitpop) goto 260
                Zix=Zindex(Zcomp,Ncomp,type)
                Nix=Nindex(Zcomp,Ncomp,type)
                do 270 nexout=0,maxex(Zix,Nix)
                  feedexcl(Zcomp,Ncomp,type,nex,nexout)=
     +              mcontrib(type,nex,nexout)
  270           continue
  260         continue
  250       continue
          endif
          if (Zcomp.eq.0.and.Ncomp.eq.0.and.flagfission)
     +      fisfeedex(0,0,maxex(0,0)+1)=xsbinary(-1)
c
c Increase emission spectra after decay of mother excitation energy bin.
c
c compemission: subroutine compound emission spectra for continuum
c emissum     : integrated emission spectrum
c Eaverage    : average outgoing energy
c Eaveragesum : help variable
c xscomp      : compound emission spectrum
c xsmpreeq    : multiple pre-equilibrium emission spectrum
c flagddx     : flag for output of double-differential cross sections
c flagrecoil  : flag for calculation of recoils
c egrid       : outgoing energy grid
c Eout        : outgoing energy
c nanglecont  : number of angles for continuum
c ang         : angle
c anglecont   : angle in degrees for continuum
c deg2rad     : conversion factor for degrees to radians
c xscompad    : compound emission angular distribution
c xsmpreeqad  : multiple preequilibrium angular distribution
c kalbach     : Kalbach function
c Einc        : incident energy in MeV
c deltaE      : energy bin around outgoing energies
c
          if (flagspec) then
            call compemission(Zcomp,Ncomp,nex)
            do 310 type=0,6
              if (parskip(type)) goto 310
              emissum(type)=0.
              Eaverage(type)=0.
              Eaveragesum=0.
              do 320 nen=ebegin(type),eend(type)
                xscomp(type,nen)=xscomp(type,nen)+xsemis(type,nen)
                xsmpreeq(type,nen)=xsmpreeq(type,nen)+
     +            xsmpeemis(type,nen)
                if (flagddx.or.flagrecoil) then
                  Eout=egrid(nen)
                  do 330 iang=0,nanglecont
                    ang=anglecont(iang)*deg2rad
                    xscompad(type,nen,iang)=xscompad(type,nen,iang)+
     +                xsemis(type,nen)/fourpi
                    xsmpreeqad(type,nen,iang)=xsmpreeqad(type,nen,iang)+
     +                xsmpeemis(type,nen)*kalbach(type,Einc,Eout,ang)
  330             continue
                endif
                emissum(type)=emissum(type)+
     +            (xsemis(type,nen)+xsmpeemis(type,nen))*deltaE(nen)
                Eaveragesum=Eaveragesum+egrid(nen)*
     +            (xsemis(type,nen)+xsmpeemis(type,nen))*deltaE(nen)
  320         continue
              if (emissum(type).gt.0.)
     +          Eaverage(type)=Eaveragesum/emissum(type)
c
c For ENDF-6 files, exclusive (n,gn), (n,gp) ... (n,ga) spectra are
c required. This is stored here.
c
c flagendf  : flag for information for ENDF-6 file
c parinclude: logical to include outgoing particle
c xsngnspec : total (projectile,gamma-ejectile) spectrum
c
              if (flagendf.and.parinclude(0).and.
     +          Zcomp.eq.0.and.Ncomp.eq.0) then
                  do 340 nen=ebegin(0),eendhigh
                    xsngnspec(type,nen)=xsemis(type,nen)+
     +                xsmpeemis(type,nen)
  340             continue
              endif
  310       continue
          endif
c
c **** Write partial emission channels and production cross sections ***
c
c Particles
c
c sumxs: sum over emission channels
c xsdif: difference between in-flux and out-flux per bin
c
          Z=ZZ(Zcomp,Ncomp,0)
          N=NN(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          if (flagpop) then
            write(*,'(/" Emitted flux per excitation energy bin of",
     +        " Z=",i3," N=",i3," (",i3,a2,"):"/)') Z,N,A,nuc(Z)
            write(*,'(" bin    Ex    ",7(2x,a8,2x),
     +        "  Total     In - out"/)') (parname(type),type=0,6)
            do 410 nex=0,maxex(Zcomp,Ncomp)
              sumxs=0.
              do 415 type=0,6
                sumxs=sumxs+xspartial(type,nex)
  415         continue
              xsdif=xspopsave(nex)-sumxs
              write(*,'(1x,i3,f8.3,1p,9e12.5)') nex,Ex(Zcomp,Ncomp,nex),
     +          (xspartial(type,nex),type=0,6),sumxs,xsdif
  410       continue
c
c Fission
c
            if (flagfission) then
              write(*,'(/" Fission contribution from Z=",i3,
     +          " N=",i3," (",i3,a2,"):"/)') Z,N,A,nuc(Z)
              write(*,'("   Ex    Popul. "/)')
              do 420 nex=0,maxex(Zcomp,Ncomp)
                write(*,'(1x,f8.3,1p,e12.5)') Ex(Zcomp,Ncomp,nex),
     +            fisfeedex(Zcomp,Ncomp,nex)
  420         continue
            endif
c
c Multiple pre-equilibrium emission
c
c p       : particle number
c h       : hole number
c xspopph : population cross section per particle-hole configuration
c xspopph2: population cross section per two-component particle-hole
c           configuration
c
            if (mulpreZN(Zcomp,Ncomp)) then
              write(*,'(/" Multiple preequilibrium emission from ",
     +          "Z=",i3," N=",i3," (",i3,a2,"):"/)') Z,N,A,nuc(Z)
              write(*,'(61x,"Feeding terms from previous ",
     +          "particle-hole configuration"/)')
              if (.not.flag2comp) then
                write(*,'(" bin    Ex  Mpe ratio  neutron   proton",
     +            "  ",10("  ",i1,"p",i1,"h    "))')
     +            ((p,h,p=1,h),h=1,4)
                write(*,'("                      emission  emission"/)')
                do 430 nex=Nlast(Zcomp,Ncomp,0)+1,maxex(Zcomp,Ncomp)
                  write(*,'(1x,i3,f8.3,f8.5,1p,12e10.3)') nex,
     +              Ex(Zcomp,Ncomp,nex),Dmulti(nex),xsmpe(1,nex),
     +              xsmpe(2,nex),
     +              ((xspopph(Zcomp,Ncomp,nex,p,h),p=1,h),h=1,4)
                  Dmulti(nex)=0.
  430           continue
              else
                write(*,'(" bin    Ex  Mpe ratio  neutron   proton",
     +            11(2x,4i2))') 1,1,0,0, 0,0,1,1, 1,0,0,1,
     +            0,1,1,0, 1,1,1,0, 1,0,1,1, 2,1,0,0, 0,0,2,1, 2,2,0,0,
     +            0,0,2,2, 1,1,1,1
                write(*,'("                      emission  emission"/)')
                do 440 nex=Nlast(Zcomp,Ncomp,0)+1,maxex(Zcomp,Ncomp)
                  write(*,'(1x,i3,f8.3,f8.5,1p,13e10.3)') nex,
     +              Ex(Zcomp,Ncomp,nex),Dmulti(nex),xsmpe(1,nex),
     +              xsmpe(2,nex),
     +              xspopph2(Zcomp,Ncomp,nex,1,1,0,0),
     +              xspopph2(Zcomp,Ncomp,nex,0,0,1,1),
     +              xspopph2(Zcomp,Ncomp,nex,1,0,0,1),
     +              xspopph2(Zcomp,Ncomp,nex,0,1,1,0),
     +              xspopph2(Zcomp,Ncomp,nex,1,1,1,0),
     +              xspopph2(Zcomp,Ncomp,nex,1,0,1,1),
     +              xspopph2(Zcomp,Ncomp,nex,2,1,0,0),
     +              xspopph2(Zcomp,Ncomp,nex,0,0,2,1),
     +              xspopph2(Zcomp,Ncomp,nex,2,2,0,0),
     +              xspopph2(Zcomp,Ncomp,nex,0,0,2,2),
     +              xspopph2(Zcomp,Ncomp,nex,1,1,1,1)
                  Dmulti(nex)=0.
  440           continue
              endif
              write(*,'(/"  Total             ",1p,2e10.3)')
     +          xsmpetot(1),xsmpetot(2)
            endif
c
c Total decay from mother nucleus
c
c xsfeed: cross section from compound to residual nucleus
c
            write(*,'(/" Emission cross sections to residual ",
     +        "nuclei from Z=",i3," N=",i3," (",i3,a2,"):"/)')
     +        Z,N,A,nuc(Z)
            if (flagfission)
     +        write(*,'(" fission  channel",23x,":",1p,e12.5)')
     +        xsfeed(Zcomp,Ncomp,-1)
            do 450 type=0,6
              if (parskip(type)) goto 450
              Z=ZZ(Zcomp,Ncomp,type)
              N=NN(Zcomp,Ncomp,type)
              A=AA(Zcomp,Ncomp,type)
              write(*,'(1x,a8," channel to Z=",i3," N=",i3," (",i3,a2,
     +          "):",1p,e12.5)') parname(type),Z,N,A,nuc(Z),
     +          xsfeed(Zcomp,Ncomp,type)
  450       continue
c
c Total emission spectra from mother nucleus.
c
c eendhigh : last energy point for energy grid for any particle
c flagcheck: flag for output of numerical checks
c
            Z=ZZ(Zcomp,Ncomp,0)
            N=NN(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            if (flagspec) then
              write(*,'(/" Emission spectra from Z=",i3," N=",i3,
     +          " (",i3,a2,"):"/)') Z,N,A,nuc(Z)
              write(*,'("  Energy ",7(2x,a8,2x)/)') (parname(type),
     +          type=0,6)
              do 460 nen=ebegin(0),eendhigh
                write(*,'(1x,f8.3,1p,7e12.5)') egrid(nen),
     +            (xsemis(type,nen)+xsmpeemis(type,nen),type=0,6)
  460         continue
              if (flagcheck) then
                write(*,'(/" ++++++++++ CHECK OF INTEGRATED ",
     +            "EMISSION SPECTRA ++++++++++"/)')
                write(*,'(13x,"Cross section   Integrated spectrum",
     +           "  Average emission energy"/)')
                do 470 type=0,6
                  if (parskip(type)) goto 470
                  write(*,'(1x,a8,1p,2(4x,e12.5),0p,10x,f8.3)')
     +              parname(type),xsfeed(Zcomp,Ncomp,type),
     +              emissum(type),Eaverage(type)
  470           continue
              endif
            endif
c
c Final production cross section for ground state and isomer of nucleus
c after decay.
c
            write(*,'(/" Final production cross section of Z=",i3,
     +        " N=",i3," (",i3,a2,"):"/)') Z,N,A,nuc(Z)
            write(*,'(" Total       :",1p,e12.5)') xspopnuc(Zcomp,Ncomp)
            write(*,'(" Ground state:",1p,e12.5)')
     +        xspopex(Zcomp,Ncomp,0)
            do 480 nex=1,Nlast(Zcomp,Ncomp,0)
              if (tau(Zcomp,Ncomp,nex).ne.0.)
     +          write(*,'(" Level",i3,"    :",1p,e12.5)') nex,
     +            xspopex(Zcomp,Ncomp,nex)
  480       continue
          endif
c
c Fission
c
c fisfeedex  : fission contribution from excitation energy bin
c flagfisfeed: flag for output of fission per excitation bin
c fisfeedJP  : fission contribution from excitation energy bin per J,P
c
          if (flagfission.and.flagfisfeed) then
            fisfile='fis000000.nex'
            write(fisfile(4:6),'(i3.3)') Z
            write(fisfile(7:9),'(i3.3)') A
            open (unit=1,status='unknown',file=fisfile)
            write(1,'("# Reaction: ",g12.4," MeV ",a8," on Z=",i3,
     +        " N=",i3," (",i3,a2,")")') einc,parname(k0),Ztarget,
     +        Ntarget,Atarget,nuc(Ztarget)
            write(1,'("# Fission contribution from Z=",i3,
     +        " N=",i3," (",i3,a2,")")') Z,N,A,nuc(Z)
            if (Zcomp.eq.0.and.Ncomp.eq.0) then
              nen=maxex(Zcomp,Ncomp)+1
            else
              nen=maxex(Zcomp,Ncomp)
            endif
            Jfis=0
            do 422 nex=0,nen
              do 422 parity=-1,1,2
                do 422 J=0,numJ
                  if (fisfeedJP(Zcomp,Ncomp,nex,J,parity).gt.0.)
     +              Jfis=max(Jfis,J)
  422       continue
            write(1,'("# # energies =",i3)') nen+1
            write(1,'("# # spins    =",i3)') Jfis+1
            form1='("#    Ex   Population",xx(5x,i2,"+",9x,i2,"-",4x))'
            write(form1(25:26),'(i2.2)') Jfis+1
            form2='(1x,f8.3,1p,e12.5,xx(2e12.5))'
            write(form2(19:20),'(i2.2)') Jfis+1
            write(1,fmt=form1) (J,J,J=0,Jfis)
            do 425 nex=0,nen
              write(1,fmt=form2) Ex(Zcomp,Ncomp,nex),
     +          fisfeedex(Zcomp,Ncomp,nex),
     +          ((fisfeedJP(Zcomp,Ncomp,nex,J,parity),parity=-1,1,2),
     +          J=0,Jfis)
  425       continue
            close (unit=1)
          endif
c
c ******* Add binary cross sections to initial compound nucleus ********
c
c xsngnsum  : sum over total (projectile,gamma-ejectile) cross sections
c xsngn     : total (projectile,gamma-ejectile) cross section
c xsbinary  : cross section from initial compound to residual nucleus
c xsreacinc : reaction cross section for incident channel
c feedbinary: feeding from first compound nucleus
c
c This is done to keep proper track of exclusive cross sections and
c channels such as (n,gn).
c
  500     if (Zcomp.eq.0.and.Ncomp.eq.0.and..not.flaginitpop) then
            xsngnsum=0.
            do 510 type=-1,6
              if (parinclude(type)) then
                xsngn(type)=xsfeed(0,0,type)
                if (type.ne.0) xsngnsum=xsngnsum+xsngn(type)
                xsfeed(0,0,type)=xsfeed(0,0,type)+xsbinary(type)
              endif
  510       continue
            if (flagchannels) then
              popexcl(0,0,maxex(0,0)+1)=xsreacinc
              if (flagfission) fisfeedex(0,0,maxex(0,0)+1)=xsbinary(-1)
              do 520 type=0,6
                if (parskip(type)) goto 520
                Zix=Zindex(Zcomp,Ncomp,type)
                Nix=Nindex(Zcomp,Ncomp,type)
                do 530 nexout=0,maxex(Zix,Nix)
                  feedexcl(0,0,type,maxex(0,0)+1,nexout)=
     +              feedbinary(type,nexout)
  530           continue
  520         continue
            endif
          endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
