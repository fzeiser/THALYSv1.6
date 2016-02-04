      subroutine gamdisout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Output of discrete gamma-ray intensities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*19 gamfile
      integer      Zcomp,Ncomp,i1,i2,Z,A,nen
      real         Egam
c
c **************** Discrete gamma-ray intensities **********************
c
c Zcomp    : charge number index for compound nucleus
c maxZ     : maximal number of protons away from the initial compound
c            nucleus
c Ncomp    : neutron number index for compound nucleus
c maxN     : maximal number of neutrons away from the initial compound
c            nucleus
c xspopnuc : population cross section per nucleus
c popeps   : limit for population cross section per nucleus
c ZZ,Z     : charge number of residual nucleus
c AA,A     : mass number of residual nucleus
c nuc      : symbol of nucleus
c numlev   : maximum number of included discrete levels
c xsgamdis : discrete gamma-ray cross section
c Egam     : outgoing energy
c edis     : energy of level
c jdis     : spin of level
c cparity  : parity (character)
c parlev   : parity of level
c
      write(*,'(/" 10. Gamma-ray intensities")')
      do 10 Zcomp=0,maxZ
        do 10 Ncomp=0,maxN
          if (xspopnuc(Zcomp,Ncomp).lt.popeps) goto 10
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          write(*,'(/" Nuclide: ",i3,a2/)') A,nuc(Z)
          write(*,'("     Initial level          Final level",
     +      "     Gamma Energy  Cross section "/)')
          write(*,'("  no.  J/Pi    Ex         no.  J/Pi    Ex"/)')
          do 30 i1=0,numlev
            do 30 i2=0,i1
              if (xsgamdis(Zcomp,Ncomp,i1,i2).eq.0.) goto 30
              Egam=edis(Zcomp,Ncomp,i1)-edis(Zcomp,Ncomp,i2)
              write(*,'(1x,i3,2x,f4.1,a1,f8.4,"  --->",i3,2x,f4.1,a1,
     +          f8.4,f11.5,1p,e15.5)') i1,jdis(Zcomp,Ncomp,i1),
     +          cparity(parlev(Zcomp,Ncomp,i1)),edis(Zcomp,Ncomp,i1),i2,
     +          jdis(Zcomp,Ncomp,i2),cparity(parlev(Zcomp,Ncomp,i2)),
     +          edis(Zcomp,Ncomp,i2),Egam,xsgamdis(Zcomp,Ncomp,i1,i2)
   30     continue
          write(*,'(/"  Total",47x,1p,e15.5)') xsgamdistot(Zcomp,Ncomp)
   10 continue
c
c Write results on separate files
c
c filegamdis: flag for gamma-ray intensities on separate file
c gamexist  : flag for existence of gamma production cross section
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c nuc       : symbol of nucleus
c Ztarget   : charge number of target nucleus
c Starget   : symbol of target nucleus
c Ethresh   : threshold incident energy for residual nucleus
c numinc    : number of incident energies
c eninc,Einc: incident energy in MeV
c numinclow : number of incident energies below Elow
c nin       : counter for incident energy
c
      if (filegamdis) then
        do 110 Zcomp=0,maxZ
          do 110 Ncomp=0,maxN
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            do 120 i1=0,numlev
              do 120 i2=0,i1
                if (xsgamdis(Zcomp,Ncomp,i1,i2).eq.0..and.
     +            .not.gamexist(Zcomp,Ncomp,i1,i2)) goto 120
                Egam=edis(Zcomp,Ncomp,i1)-edis(Zcomp,Ncomp,i2)
                gamfile='gam000000L00L00.tot'
                write(gamfile(4:9),'(2i3.3)') Z,A
                write(gamfile(11:12),'(i2.2)') i1
                write(gamfile(14:15),'(i2.2)') i2
                if (.not.gamexist(Zcomp,Ncomp,i1,i2)) then
                  gamexist(Zcomp,Ncomp,i1,i2)=.true.
                  open (unit=1,status='unknown',file=gamfile)
                  write(1,'("# ",a1," + ",i3,a2,
     +              ": Gamma-ray intensity - ",i3,a2,": Level",i3,
     +              " --> Level",i3," - gamma energy ",f11.6)')
     +              parsym(k0),Atarget,Starget,A,nuc(Z),i1,i2,Egam
                  write(1,'("# E-initial  =",f11.6," E-final=",f11.6)')
     +              edis(Zcomp,Ncomp,i1),edis(Zcomp,Ncomp,i2)
                  write(1,'("# E-threshold=",1p,e12.5)')
     +              Ethresh(Zcomp,Ncomp,i1)
                  write(1,'("# # energies =",i3)') numinc
                  write(1,'("#    E         xs")')
                  do 130 nen=1,numinclow
                    write(1,'(1p,e10.3,e12.5)') eninc(nen),0.
  130             continue
                  do 140 nen=numinclow+1,nin-1
                    write(1,'(1p,e10.3,e12.5)') eninc(nen),0.
  140             continue
                else
                  open (unit=1,status='old',file=gamfile)
                  do 150 nen=1,nin+4
                    read(1,*)
  150             continue
                endif
                write(1,'(1p,e10.3,e12.5)') Einc,
     +            xsgamdis(Zcomp,Ncomp,i1,i2)
                close (unit=1)
  120       continue
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
