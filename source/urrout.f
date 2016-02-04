      subroutine urrout
c
c +---------------------------------------------------------------------
c | Author: Gilles Noguere and Arjan Koning
c | Date  : October 23, 2013
c | Task  : Output of unresolved resonance parameters in separate files
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*3  lstring
      character*20 urrfile
      character*35 ufor2
      character*51 ufor1
      integer      odd,l,J,type,k,nen
      real         x(0:numl,0:numJ),y(0:numl,0:numJ),xx(4)
c
c General output file
c
c numinclow  : number of incident energies below Elow
c flagurrendf: flag for URR info for ENDF
c nin        : counter for incident energy
c Ztarget    : charge number of target nucleus
c Atarget    : mass number of target nucleus
c Starget    : symbol of target nucleus
c Ltarget    : excited level of target
c jdis       : spin of level
c xscaptherm : thermal capture cross section
c S          : neutron separation energy
c Rprime     : potential scattering radius
c eninc,Einc : incident energy in MeV
c odd        : odd (1) or even (0) nucleus
c lurr       : maximal orbital angular momentum for URR calculation
c lminU,lmaxU: minimal and maximal orbital angular momentum
c JminU,JmaxU: minimal and maximal total angular momentum
c Dlj        : mean resonance spacing per J,l value
c Dl         : mean resonance spacing per l value
c Purrlj     : (l,j) parity for URR calculation
c strengthlj : (l,j) neutron strength function
c strengthl  : l neutron strength function
c urrwidth   : channel width in URR
c numinc     : number of incident energies
c
      if (nin.eq.numinclow+1) then
        flagurrendf=.true.
        open (unit=21,status='unknown',file='urr.dat')
        write(21,'("#")')
        write(21,'("# Resonance parameters for Z=",i4," A=",i4," (",
     +    i3,a2,") Target spin=",f4.1)') Ztarget,Atarget,Atarget,
     +    Starget,jdis(0,1,Ltarget)
        write(21,'("# Thermal capture cross section=",1pe12.5," mb",
     +    "   Sn=",e12.5," MeV")') xscaptherm,S(0,0,1)
        write(21,'("#")')
      endif
      write(21,'("#  Einc[MeV]=",1pe12.5)') Einc
      write(21,'("# Rprime[fm]=",1pe12.5)') Rprime
      write(21,'("# l   J  P   D(l)[eV]   D(l,J)[eV]    S(l)        ",
     +  "S(l,J)   Gx(l,J)[eV] Gn(l,J)[eV] Gg(l,J)[eV] Gf(l,J)[eV]",
     +  " Gp(l,J)[eV] Ga(l,J)[eV]")')
      odd=mod(Atarget+1,2)
      do 10 l=0,lurr
        do 20 J=JminU(l),JmaxU(l)
          write(21,'(i3,f5.1,1x,a1,1p,10e12.5)') l,J+0.5*odd,
     +      cparity(Purrlj(l,J)),Dl(l),Dlj(l,J),strengthl(l),
     +      strengthlj(l,J),urrwidth(1,l,J),urrwidth(3,l,J),
     +      urrwidth(0,l,J),urrwidth(-1,l,J),urrwidth(2,l,J),
     +      urrwidth(6,l,J)
   20   continue
   10 continue
      write(21,'("#")')
      if (nin.eq.numinc) close (unit=21)
c
c Output of (l,J) dependent widths, spacings and strength functions
c in separate files
c
c urrfile: file with URR parameters
c numl   : maximal number of l-values
c nulj   : (l,j) degree of freedom
c x,y    : help variables
c Q      : Q-value
c
      do 110 type=-1,6
        if (type.eq.5.and.Q(2).le.0.) goto 110
        if (type.eq.6.and.Q(6).le.0.) goto 110
        do 120 l=0,lurr
          do 130 J=0,numJ
            x(l,J)=0.
            y(l,J)=urrwidth(type,l,J)
            if (type.eq.-1.or.type.eq.1) x(l,J)=nulj(type,l,J)
            if (type.eq.2) y(l,J)=Dlj(l,J)
            if (type.eq.3) x(l,J)=nulj(0,l,J)
            if (type.eq.4) y(l,J)=strengthlj(l,J)
            if (type.eq.5) then
              x(l,J)=nulj(2,l,J)
              y(l,J)=urrwidth(2,l,J)
            endif
            if (type.eq.6) x(l,J)=nulj(6,l,J)
  130     continue
  120   continue
        do 140 l=lminU,min(lmaxU,lurr)
          urrfile='                    '
          lstring='l00'
          write(lstring(2:3),'(i2.2)') l
          if (type.eq.-1) urrfile='urrfiswidth.'//lstring
          if (type.eq.0) urrfile='urrgamwidth.'//lstring
          if (type.eq.1) urrfile='urrcomwidth.'//lstring
          if (type.eq.2) urrfile='urrspacinglj.'//lstring
          if (type.eq.3) urrfile='urrneuwidth.'//lstring
          if (type.eq.4) urrfile='urrneustrengthlj.'//lstring
          if (type.eq.5) urrfile='urrprowidth.'//lstring
          if (type.eq.6) urrfile='urralpwidth.'//lstring
          k=1+JmaxU(l)-JminU(l)
          ufor1='("# E [MeV]   ",nn(" J      nu        width     "))'
          write(ufor1(17:18),'(i2)') k
          ufor2='(1p,e11.3,0p,nn(f4.1,1p,2e12.5,0p))'
          write(ufor2(14:15),'(i2)') k
          if (.not.urrexist(type,l)) then
            urrexist(type,l)=.true.
            open (unit=1,status='unknown',file=urrfile)
            write(1,'("# ",a1," + ",i3,a2," l-value:",i2)')
     +        parsym(k0),Atarget,Starget,l
            write(1,'("# URR parameters for ENDF-6 format")')
            if (type.eq.-1) write(1,'("# Average fission width")')
            if (type.eq.0) write(1,'("# Average radiation width")')
            if (type.eq.1) write(1,'("# Average competitive width")')
            if (type.eq.2) write(1,'("# Mean level spacing per l,J")')
            if (type.eq.3) write(1,'("# Reduced neutron width")')
            if (type.eq.4) write(1,'("# Neutron strength function")')
            if (type.eq.5) write(1,'("# Average proton width")')
            if (type.eq.6) write(1,'("# Average alpha width")')
            write(1,'("# # energies =",i3)') numinc
            if (type.eq.2) write(ufor1(38:44),'("spacing")')
            if (type.eq.4) write(ufor1(38:45),'("strength")')
            write(1,fmt=ufor1)
            do 150 nen=1,numinclow
              write(1,fmt=ufor2) eninc(nen),
     +          (J+0.5*odd,x(l,J),y(l,J),J=JminU(l),JmaxU(l))
  150       continue
            do 160 nen=numinclow+1,nin-1
              write(1,'(1p,e11.3," !!! not calculated")') eninc(nen)
  160       continue
          else
            open (unit=1,status='old',file=urrfile)
            do 170 nen=1,nin+4
              read(1,*,end=200,err=200)
  170       continue
          endif
          write(1,fmt=ufor2) Einc,
     +      (J+0.5*odd,x(l,J),y(l,J),J=JminU(l),JmaxU(l))
  200     close (unit=1)
  140   continue
  110 continue
c
c Output of l-dependent mean level spacing and neutron strength
c function in separate file
c
      do 210 type=7,8
        do 220 l=lminU,min(lmaxU,lurr)
          urrfile='                    '
          lstring='l00'
          write(lstring(2:3),'(i2.2)') l
          if (type.eq.7) then
            urrfile='urrspacingl.'//lstring
            xx(1)=Dl(l)
          else
            urrfile='urrneustrengthl.'//lstring
            xx(1)=strengthl(l)
          endif
          if (.not.urrexist(type,l)) then
            urrexist(type,l)=.true.
            open (unit=1,status='unknown',file=urrfile)
            write(1,'("# ",a1," + ",i3,a2," l-value:",i2)')
     +        parsym(k0),Atarget,Starget,l
            write(1,'("# URR parameters for ENDF-6 format")')
            if (type.eq.7) then
              write(1,'("# Mean level spacing per l")')
            else
              write(1,'("# Neutron strength function per l")')
            endif
            write(1,'("# # energies =",i3)') numinc
            if (type.eq.7) then
              write(1,'("# E [MeV]     spacing")')
            else
              write(1,'("# E [MeV]     strength")')
            endif
            do 230 nen=1,numinclow
              write(1,'(1p,e11.3,e12.5)') eninc(nen),xx(1)
  230       continue
            do 240 nen=numinclow+1,nin-1
              write(1,'(1p,e11.3," !!! not calculated")') eninc(nen)
  240       continue
          else
            open (unit=1,status='old',file=urrfile)
            do 250 nen=1,nin+4
              read(1,*,end=300,err=300)
  250       continue
          endif
          write(1,'(1p,e11.3,e12.5)') Einc,xx(1)
  300     close (unit=1)
  220   continue
  210 continue
c
c Output of URR cross section from NJOY method and TALYS
c
c flagurrnjoy : normalization of URR parameters with NJOY method
c RprimeU     : potential scattering radius
c
      do 310 type=9,11
        if (type.ne.10.and..not.flagurrnjoy) goto 310
        if (type.eq.9) then
          urrfile='urrnjoy.tot         '
          do 320 j=1,4
            xx(j)=xsurrN(j)
  320     continue
        endif
        if (type.eq.10) then
          urrfile='urrtalys.tot        '
          do 330 j=1,4
            xx(j)=xsurrT(j)
  330     continue
        endif
        if (type.eq.11) then
          urrfile='urrratio.tot        '
          do 340 j=1,4
            if (xsurrN(j).gt.0.) then
              xx(j)=xsurrT(j)/xsurrN(j)
            else
              xx(j)=1.
            endif
  340     continue
        endif
        if (.not.urrexist(type,0)) then
          urrexist(type,0)=.true.
          open (unit=1,status='unknown',file=urrfile)
          write(1,'("# ",a1," + ",i3,a2)') parsym(k0),Atarget,Starget
          if (type.eq.9) then
            write(1,'("# URR cross section with formalism from NJOY")')
            write(1,'("# Rprime[fm] = ",1pe12.5)') RprimeU
          endif
          if (type.eq.10) then
            write(1,'("# URR cross section with formalism from TALYS")')
            write(1,'("# Rprime[fm] = ",1pe12.5)') Rprime
          endif
          if (type.eq.11) then
            write(1,'("# URR cross section ratio TALYS:NJOY")')
            write(1,'("# R ratio    = ",1pe12.5)') Rprime/RprimeU
          endif
          write(1,'("# energies   = ",i3)') numinc
          write(1,'("# E [MeV]   total       elastic    ",
     +      " capture     fission")')
          do 350 nen=1,numinclow
            write(1,'(1p,e11.3,4e12.5)') eninc(nen),xx(1),xx(2),xx(4),
     +        xx(3)
  350     continue
          do 360 nen=numinclow+1,nin-1
            write(1,'(1p,e11.3," !!! not calculated")') eninc(nen)
  360     continue
        else
          open (unit=1,status='old',file=urrfile)
          do 370 nen=1,nin+4
            read(1,*,end=400,err=400)
  370     continue
        endif
        write(1,'(1p,e11.3,4e12.5)') Einc,xx(1),xx(2),xx(4),xx(3)
  400   close (unit=1)
  310 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
