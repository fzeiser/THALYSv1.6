      subroutine normalization
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 23, 2012
c | Task  : Normalize cross sections to experimental or evaluated data
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer mt,MTchan(nummt),is,nen,mtf,idc,Zix,Nix,type,iyield,
     +        iiso,i1,i2,mt0,mtc,mtd,nex
      real    Efac,xsdiftot,xsdifelas,xsadd,xsdif,ratio,R,xsfrac,
     +        xsdifiso,ratioiso,xsdifgs,ratiogs
c
c ************************ Indices for MT numbers **********************
c
c Indices for MT numbers
c
c nummt : maximum number of MT numbers
c MTchan: channel index for MT number
c
      do 10 mt=1,nummt
        MTchan(mt)=-1000
   10 continue
      MTchan(4)=100000
      MTchan(11)=201000
      MTchan(16)=200000
      MTchan(17)=300000
      MTchan(22)=100001
      MTchan(23)=100003
      MTchan(24)=200001
      MTchan(25)=300001
      MTchan(28)=110000
      MTchan(29)=100002
      MTchan(30)=200002
      MTchan(32)=101000
      MTchan(33)=100100
      MTchan(34)=100010
      MTchan(35)=101002
      MTchan(36)=100102
      MTchan(37)=400000
      MTchan(41)=210000
      MTchan(42)=310000
      MTchan(44)=120000
      MTchan(45)=110001
      MTchan(102)=000000
      MTchan(103)=010000
      MTchan(104)=001000
      MTchan(105)=000100
      MTchan(106)=000010
      MTchan(107)=000001
      MTchan(108)=000002
      MTchan(109)=000003
      MTchan(111)=020000
      MTchan(112)=010001
      MTchan(113)=000102
      MTchan(114)=001002
      MTchan(115)=011000
      MTchan(116)=010100
      MTchan(117)=001001
      MTchan(152)=500000
      MTchan(153)=600000
      MTchan(154)=200100
      MTchan(155)=000101
      MTchan(156)=410000
      MTchan(157)=301000
      MTchan(158)=101001
      MTchan(159)=210001
      MTchan(160)=700000
      MTchan(161)=800000
      MTchan(162)=510000
      MTchan(163)=610000
      MTchan(164)=710000
      MTchan(165)=400001
      MTchan(166)=500001
      MTchan(167)=600001
      MTchan(168)=700001
      MTchan(169)=401000
      MTchan(170)=501000
      MTchan(171)=601000
      MTchan(172)=300100
      MTchan(173)=400100
      MTchan(174)=500100
      MTchan(175)=600100
      MTchan(176)=200010
      MTchan(177)=300010
      MTchan(178)=400010
      MTchan(179)=320000
      MTchan(180)=300002
      MTchan(181)=310001
      MTchan(182)=001100
      MTchan(183)=111000
      MTchan(184)=110100
      MTchan(185)=101100
      MTchan(186)=110010
      MTchan(187)=101010
      MTchan(188)=100110
      MTchan(189)=100101
      MTchan(190)=220000
      MTchan(191)=010010
      MTchan(192)=001010
      MTchan(193)=000011
      MTchan(194)=420000
      MTchan(195)=400002
      MTchan(196)=410001
      MTchan(197)=030000
      MTchan(198)=130000
      MTchan(199)=320001
      MTchan(200)=520000
c
c ************************ Adjustment factors **************************
c
c Set incident energy dependent adjustment factors (purely for
c fitting purposes).
c
c Crescue  : adjustment factor for this incident energy
c Nrescue  : number of energies for adjustment factors
c numinclow: number of incident energies below Elow
c Einc     : incident energy in MeV
c Erescue  : energy grid for adjustment factors
c frescue  : adjustment factor
c
      do 110 mt=1,nummt
        do 120 is=-1,1
          Crescue(mt,is)=1.
          if (Nrescue(mt,is).eq.0) goto 120
          if (nin.eq.numinclow+1.or.Einc.le.Erescue(mt,is,1)) then
            Crescue(mt,is)=frescue(mt,is,1)
            goto 140
          endif
          if (Einc.ge.Erescue(mt,is,Nrescue(mt,is))) then
            Crescue(mt,is)=frescue(mt,is,Nrescue(mt,is))
            goto 140
          endif
          do 130 nen=1,Nrescue(mt,is)-1
            if (Einc.gt.Erescue(mt,is,nen).and.
     +        Einc.le.Erescue(mt,is,nen+1)) then
              Efac=(Einc-Erescue(mt,is,nen))/
     +          (Erescue(mt,is,nen+1)-Erescue(mt,is,nen))
              Crescue(mt,is)=frescue(mt,is,nen)+
     +          Efac*(frescue(mt,is,nen+1)-frescue(mt,is,nen))
              if (frescue(mt,is,nen+1).gt.1.e10)
     +          Crescue(mt,is)=frescue(mt,is,nen)
              if (frescue(mt,is,nen).gt.1.e10)
     +          Crescue(mt,is)=frescue(mt,is,nen+1)
              goto 140
            endif
  130     continue
  140     if (Crescue(mt,is).eq.0.) Crescue(mt,is)=1.
  120   continue
  110 continue
c
c *************************** Normalization ****************************
c
c xsdif,......: difference in cross section after normalization
c xstotinc    : total cross section (neutrons only) for incident channel
c xselastot   : total elastic cross section (shape + compound)
c xsadd       : total difference in cross section after normalization
c ratio       : adjustment factor
c
c Application of normalization from rescue files. For each partial
c cross section the applied difference is stored, so that later it can
c be redistributed.
c
      if (Crescue(1,-1).ne.1.) xsdiftot=xstotinc*(1./Crescue(1,-1)-1.)
      if (Crescue(2,-1).ne.1.) xsdifelas=xselastot*(1./Crescue(2,-1)-1.)
      xsadd=0.
      do 210 mt=4,nummt
        if (Crescue(mt,-1).eq.1..and.Crescue(mt,0).eq.1..and.
     +    Crescue(mt,1).eq.1.) goto 210
        iiso=0
        xsdifiso=0.
        ratioiso=1.
        do 220 is=1,-1,-1
          if (is.ge.0.and.Crescue(mt,0).eq.1..and.Crescue(mt,1).eq.1.)
     +      goto 220
          xsdif=0.
          ratio=1./Crescue(mt,is)
c
c Fission
c
c flagfission : flag for fission
c xsfischannel: fission channel cross section
c xsfistot    : total fission cross section
c idc         : help variables
c idnum       : counter for exclusive channel
c idchannel   : identifier for exclusive channel
c
          if (flagfission.and.mt.eq.18.and.xsfistot.gt.0) then
            xsdif=xsfistot*(ratio-1.)
            xsadd=xsadd+xsdif
            xsfistot=xsfistot*ratio
            do 230 mtf=19,38
              if (mtf.gt.21.and.mtf.lt.38) goto 230
              do 240 idc=0,idnum
                if ((mtf.eq.19.and.idchannel(idc).eq.000000).or.
     +            (mtf.eq.20.and.idchannel(idc).eq.100000).or.
     +            (mtf.eq.21.and.idchannel(idc).eq.200000).or.
     +            (mtf.eq.22.and.idchannel(idc).eq.300000)) then
                  R=ratio
                  if (Crescue(mtf,is).ne.1.) R=1./Crescue(mtf,is)
                  xsfischannel(idc)=xsfischannel(idc)*R
                endif
  240         continue
  230       continue
            goto 220
          endif
c
c Partial channels
c
c xschannel: channel cross section
c
          do 250 idc=0,idnum
            if (idchannel(idc).ne.MTchan(mt)) goto 250
            if (xschannel(idc).le.0.) goto 250
c
c Find associated residual nucleus
c
            Zix=0
            Nix=0
            do 260 type=1,6
              iyield=mod(idchannel(idc),10**(7-type))/(10**(6-type))
              Zix=Zix+iyield*parZ(type)
              Nix=Nix+iyield*parN(type)
  260       continue
            if (xspopnuc(Zix,Nix).gt.0.) then
              xsfrac=xschannel(idc)/xspopnuc(Zix,Nix)
            else
              xsfrac=1.
            endif
c
c Normalize ground state
c
c ratiogs  : adjustment factor for ground state
c xschaniso: channel cross section per isomer
c xspopex  : population cross section summed over spin and parity
c
            if (is.eq.0.and.Crescue(mt,is).ne.1.) then
              ratiogs=1./Crescue(mt,is)
              xsdifgs=xschaniso(idc,0)*(ratiogs-1.)
              xschaniso(idc,0)=xschaniso(idc,0)*ratiogs
              R=1.+(ratioiso-1.)*xsfrac
              xspopex(Zix,Nix,0)=xspopex(Zix,Nix,0)*R
            endif
c
c Normalize isomer
c
c iiso      : isomeric state
c ratioiso  : adjustment factor for isomer
c xspopnuc  : population cross section per nucleus
c flagcompo : flag for output of cross section components
c xspopdir  : direct population cross section per nucleus
c xspoppreeq: preequilibrium population cross section per nucleus
c xspopcomp : compound population cross section per nucleus

c
            if (is.eq.1) then
              do 270 i1=1,numlev
                if (xschaniso(idc,i1).ne.0.) then
                  if (Crescue(mt,is).ne.1.) then
                    ratioiso=1./Crescue(mt,is)
                    xsdifiso=xschaniso(idc,i1)*(ratioiso-1.)
                    xschaniso(idc,i1)=xschaniso(idc,i1)*ratioiso
                    R=1.+(ratioiso-1.)*xsfrac
                    xspopex(Zix,Nix,i1)=xspopex(Zix,Nix,i1)*R
                  endif
                  iiso=i1
                  goto 280
                endif
  270         continue
            endif
c
c Correct ground state, isomer or total cross section
c
  280       if (is.eq.-1) then
              R=1.+(ratioiso-1.)*xsfrac
              if (Crescue(mt,1).eq.1.) then
                if (Crescue(mt,0).ne.1.) then
                  xschaniso(idc,iiso)=
     +              max(xschannel(idc)-xschaniso(idc,0),0.)
                  xspopex(Zix,Nix,iiso)=
     +              max(xspopnuc(Zix,Nix)-xspopex(Zix,Nix,0),0.d0)
                else
                  if (Crescue(mt,-1).ne.1.) then
                    xsdif=xschannel(idc)*(ratio-1.)
                    xschannel(idc)=xschannel(idc)*ratio
                    xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)*R
                    do 290 i1=0,numlev
                      xschaniso(idc,i1)=xschaniso(idc,i1)*ratio
                      xspopex(Zix,Nix,i1)=xspopex(Zix,Nix,i1)*R
  290               continue
                    if (flagcompo) then
                      xspopdir(Zix,Nix)=xspopdir(Zix,Nix)*R
                      xspoppreeq(Zix,Nix)=xspoppreeq(Zix,Nix)*R
                      xspopcomp(Zix,Nix)=xspopcomp(Zix,Nix)*R
                    endif
                  endif
                endif
              else
                if (Crescue(mt,0).ne.1.) then
                  xsdif=xsdifiso+xsdifgs
                  ratio=1.+xsdif/xschannel(idc)
                  xschannel(idc)=xschannel(idc)*ratio
                  xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)*R
                  if (flagcompo) then
                    xspopdir(Zix,Nix)=xspopdir(Zix,Nix)*R
                    xspoppreeq(Zix,Nix)=xspoppreeq(Zix,Nix)*R
                    xspopcomp(Zix,Nix)=xspopcomp(Zix,Nix)*R
                  endif
                else
                  if (Crescue(mt,-1).ne.1.) then
                    xsdif=xschannel(idc)*(ratio-1.)
                    xschannel(idc)=xschannel(idc)*ratio
                    xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)*R
                    if (flagcompo) then
                      xspopdir(Zix,Nix)=xspopdir(Zix,Nix)*R
                      xspoppreeq(Zix,Nix)=xspoppreeq(Zix,Nix)*R
                      xspopcomp(Zix,Nix)=xspopcomp(Zix,Nix)*R
                    endif
                  endif
                  xschaniso(idc,0)=
     +              max(xschannel(idc)-xschaniso(idc,iiso),0.)
                  xspopex(Zix,Nix,0)=xspopex(Zix,Nix,0)*R
                endif
              endif
            endif
            if (is.ge.0) goto 220
            xsadd=xsadd+xsdif
c
c Normalization of related cross sections and/or spectra, by the
c rescue factors.
c
c xsgamchannel: gamma channel cross section
c xsgamdischan: discrete gamma channel cross section
c xschannelsp : channel cross section spectra
c iyield      : particle yield
c xsparticle  : total particle production cross section
c multiplicity: particle multiplicity
c
            xsgamchannel(idc)=xsgamchannel(idc)*ratio
            do 310 i1=0,numlev
              do 310 i2=0,numlev
                xsgamdischan(idc,i1,i2)=xsgamdischan(idc,i1,i2)*ratio
  310       continue
            if (flagspec) then
              do 320 type=0,6
                do 330 nen=0,numen
                  xschannelsp(idc,type,nen)=xschannelsp(idc,type,nen)*
     +              ratio
  330           continue
  320         continue
            endif
            do 340 type=1,6
              iyield=mod(idchannel(idc),10**(7-type))/(10**(6-type))
              xsparticle(type)=xsparticle(type)+iyield*xsdif
              multiplicity(type)=xsparticle(type)/xsreacinc
  340       continue
c
c Normalization of discrete level and continuum cross sections
c
c xsdisc       : total cross section for discrete state
c xscompdisc   : compound cross section for discrete state
c xsdirdisc    : direct cross section for discrete state
c xsexclcont   : exclusive single channel cross section for continuum
c xsngn        : total (projectile,gamma-ejectile) cross section
c xsexclusive  : exclusive single channel cross section
c xsdisctot    : total cross section summed over discrete states
c xsdirdisctot : direct cross section summed over discrete states
c xscompdisctot: compound cross section summed over discrete states
c xscompcont   : compound cross section for continuum
c xsdircont    : direct cross section for continuum
c xsconttot    : total cross section for continuum
c xsdirect     : total direct cross section
c xsbinary     : cross section from initial compound to residual nucleus
c
c Due to possible memory limitation, we allow only individual discrete
c level adjustment for inelastic scattering.
c
            if (mt.eq.4.or.(mt.ge.103.and.mt.le.107)) then
              if (mt.eq.4) then
                type=1
                mt0=50
                mtc=91
              else
                type=mt-101
              endif
              do 410 nex=0,numlev
                R=ratio
                if (type.eq.1) then
                  mtd=mt0+nex
                  if (Crescue(mtd,is).ne.1.) R=1./Crescue(mtd,is)
                endif
                xsdisc(type,nex)=xsdisc(type,nex)*R
                xscompdisc(type,nex)=xscompdisc(type,nex)*R
                xsdirdisc(type,nex)=xsdirdisc(type,nex)*R
  410         continue
              R=ratio
              if (type.eq.1.and.Crescue(mtc,is).ne.1.)
     +          R=1./Crescue(mtc,is)
              xsexclcont(type)=xsexclcont(type)*R
              xsngn(type)=xsngn(type)*ratio
              xsexclusive(type)=xsexclusive(type)*ratio
              xsdisctot(type)=xsdisctot(type)*ratio
              xsdirdisctot(type)=xsdirdisctot(type)*ratio
              xscompdisctot(type)=xscompdisctot(type)*ratio
              xscompcont(type)=xscompcont(type)*ratio
              xsdircont(type)=xsdircont(type)*ratio
              xsconttot(type)=xsconttot(type)*ratio
              xsdirect(type)=xsdirect(type)*ratio
              xsbinary(type)=xsbinary(type)*ratio
            endif
  250     continue
  220   continue
  210 continue
c
c Put difference in the elastic (or total) cross section
c
c xselastot  : total elastic cross section (shape + compound)
c xstotadjust: total cross section adjustment
c xseladjust : elastic cross section adjustment
c xsnonel    : non-elastic cross section
c xsnonadjust: nonelastic cross section adjustment
c channelsum : sum over exclusive channel cross sections
c
      if (Crescue(2,-1).ne.1.) then
        xseladjust(nin)=xsdifelas
        xselastot=xselastot+xseladjust(nin)
        xstotadjust(nin)=xsdifelas+xsadd
        xstotinc=xstotinc+xstotadjust(nin)
      else
        xstotadjust(nin)=xsdiftot
        xstotinc=xstotinc+xstotadjust(nin)
        xseladjust(nin)=xsdiftot-xsadd
        xselastot=xselastot+xseladjust(nin)
      endif
      xsnonadjust(nin)=xsadd
      xsnonel=xsnonel+xsnonadjust(nin)
      channelsum=channelsum+xsadd
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
