        subroutine integral
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Calculate effective cross section for integral spectrum
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80  fluxfile
      integer       i,istat,nen,Nspec,Nxs,nen0,is,k
      real          Efluxup(0:numenin),Eflux(0:numenin),fspec(numenin),
     +              fluxsum,Exs(0:numenin),xs(0:numenin),xseff,Efl,Ea1,
     +              Eb1,xsa,xsb,xsf,xsexp,ratio
c
c ********* Read integral spectrum from experimental database **********
c
c parsym  : symbol of particle
c k0      : index of incident particle
c Atarget : mass number of target nucleus
c Starget : symbol of target nucleus
c Ztarget : charge number of target nucleus
c Nflux   : number of reactions with integral data
c path    : directory containing structure files to be read
c lenpath : length of pathname
c fluxfile: file with experimental integral spectrum
c fluxname: name of experimental flux
c Nspec   : number of spectral energies
c fspec   : spectrum values
c
      open (unit=1,status='unknown',file='integral.dat')
      write(1,'("# ",a1," + ",i3,a2,
     +  ": Effective cross sections from integral data")')
     +  parsym(k0),Atarget,Starget
      write(1,'("# Channel      Flux          Eff. c.s. (b)",
     +  " Exp. c.s. (b)       Ratio")')
      do 10 i=1,Nflux
        fluxfile=path(1:lenpath)//'integral/flux/spectrum.'//fluxname(i)
        open (unit=2,status='old',file=fluxfile,iostat=istat)
        if (istat.eq.0) then
          read(2,'(1x,i4)') Nspec
          do 110 nen=Nspec,1,-1
            read(2,*) nen0,Efluxup(nen),fspec(nen)
  110     continue
          close (unit=2)
c
c Determine middle of histograms and energy bins. Calculate the
c integral of the flux for normalization.
c
          fluxsum=0.
          Efluxup(0)=Efluxup(1)
          do 120 nen=1,Nspec
            Eflux(nen)=0.5*(Efluxup(nen-1)+Efluxup(nen))*1.e-6
            fluxsum=fluxsum+fspec(nen)
  120     continue
        else
          write(*,'(" TALYS-warning: integral spectrum file ",a80,
     +      "does not exist")') fluxfile
          goto 10
        endif
c
c *************** Read cross sections from TALYS output files **********
c
c Nxs: number of incident energies
c
        is=1
  200   open (unit=3,status='old',file=xsfluxfile(i),iostat=istat)
        if (istat.eq.0) then
          Exs(0)=0.
          xs(0)=0.
          read(3,'(///,14x,i3,/)') Nxs
          do 210 nen=1,Nxs
            read(3,*) Exs(nen),xs(nen)
  210     continue
          close (unit=3)
        else
          if (is.lt.numlev) then
            is=is+1
            do 220 k=1,20
              if (xsfluxfile(i)(k:k).eq.'L') then
                write(xsfluxfile(i)(k+1:k+2),'(i2.2)') is
                goto 200
              endif
  220       continue
          endif
          write(*,'(" TALYS-warning: cross section file ",a80,
     +      "does not exist")') xsfluxfile(i)
          goto 10
        endif
c
c *************** Calculate effective cross section by folding *********
c
c xseff: effective cross section
c xsexp: experimental effective cross section
c
c Interpolate cross section grid on the grid of the integral spectrum
c
        xseff=0.
        do 310 nen=1,Nspec
          Efl=Eflux(nen)
          if (Efl.ge.Exs(0).and.Efl.le.Exs(Nxs)) then
            call locate(Exs,0,Nxs,Efl,nen0)
            Ea1=Exs(nen0)
            Eb1=Exs(nen0+1)
            xsa=xs(nen0)
            xsb=xs(nen0+1)
            call pol1(Ea1,Eb1,xsa,xsb,Efl,xsf)
            xseff=xseff+xsf*fspec(nen)
          endif
  310   continue
        xseff=0.001*xseff/fluxsum
        xsexp=integralexp(i)
        if (xsexp.gt.0.) then
          ratio=xseff/xsexp
        else
          ratio=1.
        endif
        write(1,'(2a15,1p,2e12.5,0p,f15.5)') xsfluxfile(i),fluxname(i),
     +    xseff,xsexp,ratio
   10 continue
      close (unit=1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
