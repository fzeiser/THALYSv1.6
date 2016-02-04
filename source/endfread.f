      subroutine endfread
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 1, 2013
c | Task  : Read ECIS results for incident particle on ENDF-6 energy
c |         grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72     line
      integer          infileendf,nen,Z,A,nend,mt,nen2
      real             tripathi,e,Ea,Eb,xsa,xsb,xsc,xsd,Efac,xsdift,
     +                 xsdife,enuc
      double precision xs
c
c ************ Read total, reaction and elastic cross section **********
c
c flagendfecis: flag for new ECIS calculation for ENDF-6 files
c nen6        : total number of energies
c e6,e        : energies of ENDF-6 energy grid in MeV
c k0          : index of incident particle
c coullimit   : energy limit for charged particle OMP calculation
c xs          : help variable
c xstot6      : total cross section (neutrons only) for ENDF-6 file
c xsreac6     : reaction cross section for ENDF-6 file
c xsopt6      : optical model reaction cross section for ENDF-6 file
c xselassh6   : shape elastic cross section (neutrons only) for ENDF-6
c               file
c
      if (flagendfecis) then
        open (unit=3,status='unknown',file='ecis.endfcs')
        infileendf=3
        open (unit=23,status='unknown',file='endf.cs')
   10   read(3,'(a72)',end=20) line
        write(23,'(a72)') line
        goto 10
   20   rewind 3
      else
        infileendf=23
      endif
      do 30 nen=1,nen6
        e=real(e6(nen))
        if (k0.gt.1.and.e.lt.coullimit(k0)) goto 30
        read(infileendf,'()')
        if (k0.eq.1) then
          read(infileendf,*) xs
          xstot6(nen)=max(real(xs),0.)
        endif
        read(infileendf,*) xs
        xsreac6(nen)=max(real(xs),0.)
        xsopt6(nen)=xsreac6(nen)
        if (k0.eq.1) then
          read(infileendf,*) xs
          xselassh6(nen)=max(real(xs),0.)
        endif
   30 continue
      close (unit=3,status=ecisstatus)
      close (unit=23,status=ecisstatus)
      open (unit=10,status='unknown',file='ecis.endfin')
      close (unit=10,status=ecisstatus)
c
c ********** Compound elastic contribution and normalization ***********
c
c locate     : subroutine to find value in ordered table
c xscompel    : compound elastic cross section
c pol1       : subroutine for polynomial interpolation of first order
c xselas6    : total elastic cross section (neutrons only) for ENDF-6
c              file
c xsnon6     : non-elastic cross section for ENDF-6 file
c flagrescue : flag for final rescue: normalization to data
c
        if (k0.ge.1) then
          do 110 nen=1,nen6
            e=real(e6(nen))
            if (e.le.eninc(1)) then
              xsc=xscompel6(1)
              xsd=xsnonel6(1)
            else
              call locate(eninc,1,numinc,e,nend)
              Ea=eninc(nend)
              Eb=eninc(nend+1)
              xsa=xscompel6(nend)
              xsb=xscompel6(nend+1)
              call pol1(Ea,Eb,xsa,xsb,e,xsc)
              xsa=xsnonel6(nend)
              xsb=xsnonel6(nend+1)
              call pol1(Ea,Eb,xsa,xsb,e,xsd)
            endif
            xselas6(nen)=xselassh6(nen)+xsc
            xsnon6(nen)=xsd
            xstot6(nen)=xselas6(nen)+xsnon6(nen)
c
c ************************ Adjustment factors **************************
c
c Set incident energy dependent adjustment factors (purely for
c fitting purposes).
c
c Nrescue: number of energies for adjustment factors
c Crescue: adjustment factor for this incident energy
c Erescue: energy grid for adjustment factors
c frescue: adjustment factor
c
            if (flagrescue) then
              do 120 mt=1,3
                if (Nrescue(mt,-1).eq.0) goto 120
                Crescue(mt,-1)=1.
                if (e.le.Erescue(mt,-1,1)) then
                  Crescue(mt,-1)=frescue(mt,-1,1)
                  goto 120
                endif
                if (e.ge.Erescue(mt,-1,Nrescue(mt,-1))) then
                  Crescue(mt,-1)=frescue(mt,-1,Nrescue(mt,-1))
                  goto 120
                endif
                do 130 nen2=1,Nrescue(mt,-1)-1
                  if (e.gt.Erescue(mt,-1,nen2).and.
     +              e.le.Erescue(mt,-1,nen2+1)) then
                    Efac=(e-Erescue(mt,-1,nen2))/
     +                (Erescue(mt,-1,nen2+1)-Erescue(mt,-1,nen2))
                    Crescue(mt,-1)=frescue(mt,-1,nen2)+
     +                Efac*(frescue(mt,-1,nen2+1)-frescue(mt,-1,nen2))
                    if (frescue(mt,-1,nen2+1).gt.1.e10)
     +                Crescue(mt,-1)=frescue(mt,-1,nen2)
                    if (frescue(mt,-1,nen2).gt.1.e10)
     +                Crescue(mt,-1)=frescue(mt,-1,nen2+1)
                    goto 120
                  endif
  130           continue
  120         continue
c
c Put difference in the elastic (or total) cross section
c
              if (Crescue(1,-1).ne.1..and.Crescue(1,-1).ne.0.)
     +          xsdift=xstot6(nen)*(1./Crescue(1,-1)-1.)
              if (Crescue(2,-1).ne.1..and.Crescue(2,-1).ne.0.)
     +          xsdife=xselas6(nen)*(1./Crescue(2,-1)-1.)
              if (Crescue(2,-1).ne.1..and.Crescue(2,-1).ne.0.) then
                xselas6(nen)=xselas6(nen)+xsdife
                xstot6(nen)=xstot6(nen)+xsdife
              else
                xselas6(nen)=xselas6(nen)+xsdift
                xstot6(nen)=xstot6(nen)+xsdift
              endif
            endif
  110     continue
        else
          do 150 nen=1,nen6
            xsnon6(nen)=xsreac6(nen)
  150     continue
        endif
c
c ************ Normalization with semi-empirical results ***************
c
c flagsys   : flag for reaction cross section from systematics
c ZZ,Z      : charge number of residual nucleus
c AA,A      : mass number of residual nucleus
c enuc      : incident energy in MeV per nucleon
c parA      : mass number of particle
c tripathi  : function for semi-empirical reaction cross section of
c             Tripathi et al.
c parZ      : charge number of particle
c threshnorm: normalization factor at trheshold
c
      if (flagsys(k0)) then
        Z=ZZ(0,0,k0)
        A=AA(0,0,k0)
        do 210 nen=1,nen6
          if (xsopt6(nen).eq.0.) goto 210
          e=real(e6(nen))
          enuc=e/parA(k0)
          xs=tripathi(parZ(k0),parA(k0),Z,A,enuc)
          if (xs.eq.0.) xs=xsopt6(nen)*threshnorm(k0)
          xsnon6(nen)=xs
          if (k0.eq.1) xselas6(nen)=xselas6(nen)+xsopt6(nen)-xs
  210   continue
      endif
c
c **************** Write total cross sections to file ******************
c
c parsym    : symbol of particle
c Atarget   : mass number of target nucleus
c Ztarget   : charge number of target nucleus
c Starget   : symbol of target nucleus
c numinclow : number of incident energies below Elow
c fxsnonel  : non-elastic cross section for incident channel
c fxselastot: total elastic cross section (neutrons only) for
c             incident channel
c fxstotinc : total cross section (neutrons only) for incident channel
c
      open (unit=1,status='unknown',file='endf.tot')
      write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +  parsym(k0),Atarget,Starget
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# # energies =",i4)') nen6+numinclow
      write(1,'("#    E        Non-elastic Elastic     Total")')
      do 310 nen=1,numinclow
        write(1,'(1p,e12.5,2x,3e11.4)') eninc(nen),fxsnonel(nen),
     +    fxselastot(nen),fxstotinc(nen)
  310 continue
      do 320 nen=1,nen6
        write(1,'(1p,e12.5,2x,3e11.4)') e6(nen),xsnon6(nen),
     +    xselas6(nen),xstot6(nen)
  320 continue
      close (unit=1)
c
c Total cross sections only
c
      open (unit=1,status='unknown',file='endftot.tot')
      write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +  parsym(k0),Atarget,Starget
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# # energies =",i4)') nen6+numinclow
      write(1,'("#    E      Cross section")')
      do 410 nen=1,numinclow
        write(1,'(1p,e12.5,2x,e11.4)') eninc(nen),fxstotinc(nen)
  410 continue
      do 420 nen=1,nen6
        write(1,'(1p,e12.5,2x,e11.4)') e6(nen),xstot6(nen)
  420 continue
      close (unit=1)
c
c Elastic cross sections only
c
      open (unit=1,status='unknown',file='endfel.tot')
      write(1,'("# ",a1," + ",i3,a2," Elastic cross sections")')
     +  parsym(k0),Atarget,Starget
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# # energies =",i4)') nen6+numinclow
      write(1,'("#    E      Cross section")')
      do 510 nen=1,numinclow
        write(1,'(1p,e12.5,2x,e11.4)') eninc(nen),fxselastot(nen)
  510 continue
      do 520 nen=1,nen6
        write(1,'(1p,e12.5,2x,e11.4)') e6(nen),xselas6(nen)
  520 continue
      close (unit=1)
c
c Nonelastic cross sections only
c
      open (unit=1,status='unknown',file='endfnon.tot')
      write(1,'("# ",a1," + ",i3,a2," Nonelastic cross sections")')
     +  parsym(k0),Atarget,Starget
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# # energies =",i4)') nen6+numinclow
      write(1,'("#    E      Cross section")')
      do 610 nen=1,numinclow
        write(1,'(1p,e12.5,2x,e11.4)') eninc(nen),fxsnonel(nen)
  610 continue
      do 620 nen=1,nen6
        write(1,'(1p,e12.5,2x,e11.4)') e6(nen),xsnon6(nen)
  620 continue
      close (unit=1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
