      subroutine prodyield
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 13, 2013
c | Task  : Calculate production yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          it,Zix,Nix,Z,N,A,is,isob,Zp,Ncool
      real             rfac,yfac,N0
      double precision Th,dT,Tc,prate0,pratei,lamD,C1,TT,exp1,denom,
     +                 prateP,lamPD,CP1,exp2,term,t1,t2,expo1,expo2
c
c ************ Initial condition for irradiation ***********************
c
c Ntar0    : number of original target atoms
c avogadro : Avogadro's number
c Atarget  : mass number of target nucleus
c rhotarget: target density
c Vtar     : active target volume
c
      Ntar0=avogadro/Atarget*rhotarget*Vtar
c
c ******************** Set time grid ***********************************
c
c Ntime: number of time points
c dT   : time step
c Tir  : irradiation time
c Tgrid: time
c
      Ntime=numtime/2
      Th=Tir/hoursec
      dT=Th/Ntime
      if (dT.le.0.1) dT=0.1
      if (dT.gt.0.1.and.dT.le.0.2) dT=0.2
      if (dT.gt.0.2.and.dT.le.0.5) dT=0.5
      if (dT.gt.0.5.and.dT.le.1.0) dT=1.0
      do 10 it=0,Ntime
        Tgrid(it)=it*dT
        if (Tgrid(it).gt.Th) then
          Ntime=it-1
          goto 20
        endif
   10 continue
   20 if (Tco.gt.0.) then
        Tc=Tco/hoursec
        Ncool=numtime-Ntime
        dT=log(Tc)/Ncool
        do 30 it=Ntime+1,numtime
          Tgrid(it)=Th+exp((it-Ntime)*dT)-1.
   30   continue
      endif
c
c ******************** Activity (yield) in MBq *************************
c
c prate,pratei: production rate per isotope
c prate0      : production rate for all isotopes
c maxZ        : maximal number of protons away from the initial
c               compound nucleus
c Z           : charge number of nucleus
c Zinit       : charge number of initial compound nucleus
c maxN        : maximal number of neutrons away from the initial
c               compound nucleus
c N           : neutron number of nucleus
c Ninit       : neutron number of initial compound nucleus
c A           : mass number of nucleus
c Nisotot     : number of elemental isotopes produced after
c               irradiation
c Nisomer     : number of isomers for this nuclide
c Tmax        : irradiation time with maximal yield
c Tp          : irradiation time with maximal yield per time unit
c Niso        : number of isotopes produced after irradiation
c activity    : activity of produced isotope in MBq
c yield       : yield of produced isotope in MBq/(mA.h)
c Nisorel     : fraction of number of produced isotopes per element
c Yexist      : flag for existence of yield
c lambda,lamD : decay rate per isotope
c C1,...      : constants
c hoursec     : number of seconds in an hour
c
      do 40 Zix=0,numZ
        do 50 it=0,numtime
          Nisotot(Zix,it)=0.
   50   continue
        do 60 Nix=0,numN
          do 70 is=-1,numisom
            Yexist(Zix,Nix,is)=.false.
            Tmax(Zix,Nix,is)=0.
            do 80 it=1,5
              Tp(Zix,Nix,is,it)=0.
   80       continue
            do 90 it=0,numtime
              Niso(Zix,Nix,is,it)=0.
              activity(Zix,Nix,is,it)=0.
              yield(Zix,Nix,is,it)=0.
              Nisorel(Zix,Nix,is,it)=0.
   90       continue
   70     continue
   60   continue
   40 continue
      prate0=dble(prate(-1,-1,-1))
      do 110 Zix=0,maxZ
        Z=Zinit-Zix
        do 130 Nix=0,maxN
          N=Ninit-Nix
          A=Z+N
          do 140 is=-1,Nisomer(Zix,Nix)
            if (Z.eq.Ztarget.and.A.eq.Atarget.and.is.eq.-1) then
              Niso(Zix,Nix,is,0)=dble(Ntar0)
              Nisotot(Zix,0)=Nisotot(Zix,0)+Niso(Zix,Nix,is,0)
            endif
            pratei=dble(prate(Zix,Nix,is))
            lamD=dble(lambda(Zix,Nix,is))
            C1=dble(Ntar0)*pratei
            do 170 it=1,numtime
              TT=dble(Tgrid(it)*hoursec)
c
c Depletion of target
c
              if (Z.eq.Ztarget.and.A.eq.Atarget.and.is.eq.-1) then
                if (it.le.Ntime) then
                  Niso(Zix,Nix,is,it)=dble(Ntar0)*exp(-prate0*TT)
                else
                  Niso(Zix,Nix,is,it)=Niso(Zix,Nix,is,Ntime)
                endif
              else
c
c Production and decay of other isotopes
c
c 1. Production directly from target
c
                t1=prate0*TT
                expo1=exp(-t1)
                t2=lamD*TT
                expo2=exp(-t2)
                if (it.le.Ntime) then
                  denom=lamD-prate0
                  exp1=expo1/denom-expo2/denom
                  Niso(Zix,Nix,is,it)=C1*exp1
                else
                  t2=lamD*(TT-Tir)
                  expo2=exp(-t2)
                  term=Niso(Zix,Nix,is,Ntime)*expo2
                  Niso(Zix,Nix,is,it)=term
                endif
c
c 2. Production from decay of other isotope
c
c rtyp      : type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
c lamP,lamPD: decay rate for parent isotope
c Ibeam     : beam current in mA
c
                do 180 isob=-1,1
                  Zp=Zix-isob
                  if (Zp.lt.0) goto 180
                  if ((isob.eq.-1.and.rtyp(Zp,Nix,-1).eq.1).or.
     +              (isob.eq.1.and.rtyp(Zp,Nix,-1).eq.2)) then
                    prateP=dble(prate(Zp,Nix,is))
                    lamPD=dble(lambda(Zp,Nix,is))
                    if (it.le.Ntime) then
                      if (prateP.gt.0..and.lamPD.gt.0.) then
                        CP1=dble(Ntar0)*prateP
                        t1=prateP*TT
                        expo1=exp(-t1)
                        denom=lamD-prateP
                        exp2=expo1/denom-expo2/denom
                        term=lamPD*CP1/(prateP-prate0)*(exp1-exp2)
                        Niso(Zix,Nix,is,it)=Niso(Zix,Nix,is,it)+term
                      endif
                    else
c
c Cooling only
c
                      N0=Niso(Zp,Nix,is,Ntime)
                      t1=lamPD*(TT-Tir)
                      exp1=exp(-t1)
                      t2=lamD*(TT-Tir)
                      exp2=exp(-t2)
                      denom=lamD-lamPD
                      if (denom.ne.0.) then
                        term=N0*lamPD/denom*(exp1-exp2)
                        Niso(Zix,Nix,is,it)=Niso(Zix,Nix,is,it)+term
                      endif
                    endif
                  endif
  180           continue
                activity(Zix,Nix,is,it)=lamD*Niso(Zix,Nix,is,it)*1.e-6
                yield(Zix,Nix,is,it)=abs(activity(Zix,Nix,is,it)-
     +            activity(Zix,Nix,is,it-1))/
     +            (Ibeam*dble(Tgrid(it)-Tgrid(it-1)))
              endif
              if (Niso(Zix,Nix,is,it).gt.0.) Yexist(Zix,Nix,is)=.true.
              if (is.eq.-1) Nisotot(Zix,it)=Nisotot(Zix,it)+
     +          Niso(Zix,Nix,is,it)
  170       continue
            if (.not.Yexist(Zix,Nix,is)) goto 140
            if (lamD.gt.0..and.prate0.gt.0.) then
              Tmax(Zix,Nix,is)=log(lamD/prate0)/(lamD-prate0)
            else
              Tmax(Zix,Nix,is)=1.e30
            endif
c
c Write irradiation time with maximum yield in years, days, etc.
c
c minutesec: number of seconds in a minute
c daysec   : number of seconds in a day
c yearsec  : number of seconds in a year
c
            TT=Tmax(Zix,Nix,is)
            Tp(Zix,Nix,is,1)=int(TT/yearsec)
            TT=TT-Tp(Zix,Nix,is,1)*yearsec
            Tp(Zix,Nix,is,2)=int(TT/daysec)
            TT=TT-Tp(Zix,Nix,is,2)*daysec
            Tp(Zix,Nix,is,3)=int(TT/hoursec)
            TT=TT-Tp(Zix,Nix,is,3)*hoursec
            Tp(Zix,Nix,is,4)=int(TT/minutesec)
            TT=TT-Tp(Zix,Nix,is,4)*minutesec
            Tp(Zix,Nix,is,5)=int(TT)
  140     continue
  130   continue
        do 190 Nix=0,maxN
          do 190 is=-1,Nisomer(Zix,Nix)
            if (.not.Yexist(Zix,Nix,is)) goto 190
            do 200 it=0,numtime
              if (Nisotot(Zix,it).ne.0.) Nisorel(Zix,Nix,is,it)=
     +          Niso(Zix,Nix,is,it)/Nisotot(Zix,it)
  200       continue
  190   continue
  110 continue
c
c Transform quantities to user-dependent units
c
c radiounit: unit for radioactivity: Bq, kBq, MBq, Gbq,
c            mCi, Ci or kCi
c rfac     : conversion factor for radioactivity
c yieldunit: unit for isotope yield: num (number),
c            mug (micro-gram), mg, g, or kg
c yfac     : conversion factor for isotope yield
c
      rfac=1.
      yfac=1.
      if (radiounit.eq.'bq') rfac=1.e6
      if (radiounit.eq.'kbq') rfac=1.e3
      if (radiounit.eq.'gbq') rfac=1.e-3
      if (radiounit.eq.'ci') rfac=1./3.7e4
      if (radiounit.eq.'kci') rfac=1./3.7e7
      if (radiounit.eq.'mci') rfac=1./3.7e1
      do 210 Zix=0,maxZ
        Z=Zinit-Zix
        do 220 Nix=0,maxN
          N=Ninit-Nix
          A=Z+N
          if (yieldunit.eq.'g') yfac=real(A)/avogadro
          if (yieldunit.eq.'mug') yfac=real(A)/avogadro*1.e6
          if (yieldunit.eq.'mg') yfac=real(A)/avogadro*1.e3
          if (yieldunit.eq.'kg') yfac=real(A)/avogadro*1.e-3
          do 230 is=-1,Nisomer(Zix,Nix)
            do 240 it=1,numtime
              activity(Zix,Nix,is,it)=rfac*activity(Zix,Nix,is,it)
              yield(Zix,Nix,is,it)=rfac*yield(Zix,Nix,is,it)
              Niso(Zix,Nix,is,it)=yfac*Niso(Zix,Nix,is,it)
  240       continue
  230     continue
  220   continue
  210 continue
      return
      end
Copyright (C) 2010  A.J. Koning
