      subroutine fissionout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Output of fission cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*12 rpfile
      integer      Acomp,Zcomp,Ncomp,Z,A,nen
c
c *********************** Fission cross sections ***********************
c
c Acomp    : mass number index for compound nucleus
c maxA     : maximal number of nucleons away from the initial compound
c            nucleus
c Zcomp    : charge number index for compound nucleus
c maxZ     : maximal number of protons away from the initial compound
c            nucleus
c Ncomp    : neutron number index for compound nucleus
c maxN     : maximal number of neutrons away from the initial compound
c            nucleus
c nin      : counter for incident energy
c numinclow: number of incident energies below Elow
c xsfeed   : cross section from compound to residual nucleus
c ZZ,Z     : charge number of residual nucleus
c AA,A     : mass number of residual nucleus
c nuc      : symbol of nucleus
c
      write(*,'(/" 4b. Fission cross section per fissioning",
     +  " nuclide"/)')
      do 10 Acomp=0,maxA
        do 20 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 20
          if (xsfeed(Zcomp,Ncomp,-1).ne.0.) then
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            write(*,'(1x,2i4," (",i3,a2,")",1p,e12.5)') Z,A,A,
     +        nuc(Z),xsfeed(Zcomp,Ncomp,-1)
          endif
   20   continue
   10 continue
c
c Write results to separate file
c
c filefission: flag for fission cross sections on separate file
c fisexist   : flag for existence of fission cross section
c parsym     : symbol of particle
c k0         : index of incident particle
c Atarget    : mass number of target nucleus
c Ztarget    : charge number of target nucleus
c Starget    : symbol of target nucleus
c numinc     : number of incident energies
c eninc,Einc : incident energy in MeV
c
      if (filefission) then
        do 110 Acomp=0,maxA
          do 120 Zcomp=0,maxZ
            Ncomp=Acomp-Zcomp
            if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 120
            if (xsfeed(Zcomp,Ncomp,-1).eq.0..and.
     +        .not.fisexist(Zcomp,Ncomp)) goto 120
            rpfile='rp000000.fis'
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            write(rpfile(3:8),'(2i3.3)') Z,A
            if (.not.fisexist(Zcomp,Ncomp)) then
              fisexist(Zcomp,Ncomp)=.true.
              open (unit=1,status='unknown',file=rpfile)
              write(1,'("# ",a1," + ",i3,a2,": Fission of ",i3,a2)')
     +          parsym(k0),Atarget,Starget,A,nuc(Z)
              write(1,'("#                 ")')
              write(1,'("#                 ")')
              write(1,'("# # energies =",i3)') numinc
              write(1,'("#    E         xs ")')
              do 130 nen=1,numinclow
                write(1,'(1p,e10.3,e12.5)') eninc(nen),0.
  130         continue
              do 140 nen=numinclow+1,nin-1
                write(1,'(1p,e10.3,e12.5)') eninc(nen),0.
  140         continue
            else
              open (unit=1,status='old',file=rpfile)
              do 150 nen=1,nin+4
                read(1,*)
  150         continue
            endif
            write(1,'(1p,e10.3,e12.5)') Einc,xsfeed(Zcomp,Ncomp,-1)
            close (unit=1)
  120     continue
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
