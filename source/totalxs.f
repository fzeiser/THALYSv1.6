      subroutine totalxs
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 11, 2012
c | Task  : Total cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,ident,idc,Zcomp,Ncomp
c
c *********************** Specific cross sections **********************
c
c flagchannels: flag for exclusive channels calculation
c parskip     : logical to skip outgoing particle
c xsexclusive : exclusive single channel cross section
c xschannel   : channel cross section
c idnum       : counter for exclusive channel
c idchannel   : identifier for exclusive channel
c xsexclcont  : exclusive single channel cross section for continuum
c xsdisctot   : total cross section summed over discrete states
c
      if (flagchannels) then
        do 10 type=0,6
          if (parskip(type)) goto 10
          if (type.eq.0) then
            xsexclusive(0)=xschannel(0)
          else
            ident=10**(6-type)
            do 20 idc=0,idnum
              if (idchannel(idc).eq.ident) then
                xsexclusive(type)=xschannel(idc)
                goto 30
              endif
   20       continue
          endif
   30     if (xsconttot(type).eq.0.) then
            xsexclcont(type)=0.
          else
            xsexclcont(type)=max(xsexclusive(type)-xsdisctot(type),0.)
          endif
   10   continue
      endif
c
c *************** Total particle production cross sections *************
c
c xsparticle  : total particle production cross section
c Zcomp       : charge number index for compound nucleus
c maxZ        : maximal number of protons away from the initial
c               compound nucleus
c Ncomp       : neutron number index for compound nucleus
c maxN        : maximal number of neutrons away from the initial
c               compound nucleus
c xsfeed      : cross section from compound to residual nucleus
c multiplicity: particle multiplicity
c xsnonel     : non-elastic cross section
c
      do 110 type=0,6
        if (parskip(type)) goto 110
        xsparticle(type)=0.
        do 120 Zcomp=0,maxZ
          do 120 Ncomp=0,maxN
            xsparticle(type)=xsparticle(type)+xsfeed(Zcomp,Ncomp,type)
  120   continue
        if (xsnonel.ne.0.)
     +    multiplicity(type)=xsparticle(type)/xsnonel
  110 continue
c
c ******************* Total fission cross sections ********************
c
c flagfission: flag for fission
c xsfistot   : total fission cross section
c xsfeed     : cross section from compound to residual nucleus
c flagastro  : flag for calculation of astrophysics reaction rate
c xsastrofis : astrophysical fission cross section
c
      if (flagfission) then
        xsfistot=0.
        do 210 Zcomp=0,maxZ
          do 210 Ncomp=0,maxN
          xsfistot=xsfistot+xsfeed(Zcomp,Ncomp,-1)
  210   continue
        if (flagastro) xsastrofis(nin)=xsfistot
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
