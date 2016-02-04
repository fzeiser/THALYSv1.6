      subroutine recoilout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 1, 2012
c | Task  : Output of recoils
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*28     recfile
      integer          Zcomp,Ncomp,Z,A,nen
      double precision sumcm
c
c ***************************** Spectra ********************************
c
c Zcomp     : charge number index for compound nucleus
c maxZ      : maximal number of protons away from the initial compound
c             nucleus
c Ncomp     : neutron number index for compound nucleus
c maxN      : maximal number of neutrons away from the initial compound
c             nucleus
c xspopnuc  : population cross section per nucleus
c xseps     : limit for cross sections
c ZZ,Z      : charge number of residual nucleus
c AA,Z      : mass number of residual nucleus
c nuc       : symbol of nucleus
c maxenrec  : number of recoil energies
c Erec      : recoil energy
c specrecoil: recoil spectrum
c recoilint : total recoil integrated over spectrum
c parZ      : charge number of particle
c k0        : index of incident particle
c parN      : neutron number of particle
c sumcm     : total residual production in the CM frame
c xselasinc : total elastic cross section (neutrons only) for incident
c             channel
c filerecoil: flag for recoil spectra on separate file
c Starget   : symbol of target nucleus
c Einc      : incident energy in MeV
c
      write(*,'(/" 8. Recoil spectra")')
      do 10 Zcomp=0,maxZ
        do 10 Ncomp=0,maxN
          if (xspopnuc(Zcomp,Ncomp).lt.xseps) goto 10
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          write(*,'(/" Recoil Spectrum for ",i3,a2/)') A,nuc(Z)
          write(*,'("  Energy   Cross section"/)')
          do 20 nen=0,maxenrec
            write(*,'(1x,f7.3,1p,e12.5)') Erec(Zcomp,Ncomp,nen),
     +        specrecoil(Zcomp,Ncomp,nen)
   20   continue
        write(*,'(/" Integrated recoil spectrum       : ",1p,e12.5)')
     +    recoilint(Zcomp,Ncomp)
        if (Zcomp.eq.parZ(k0).and.Ncomp.eq.parN(k0)) then
          sumcm=xspopnuc(Zcomp,Ncomp)+xselasinc
        else
          sumcm=xspopnuc(Zcomp,Ncomp)
        endif
        write(*,'(" Residual production cross section: ",1p,e12.5)')
     +    sumcm
c
c Write results to separate file
c
        if (filerecoil) then
          recfile='rec000000spec000.000.tot'//natstring(iso)
          write(recfile(4:9),'(2i3.3)') Z,A
          write(recfile(14:20),'(f7.3)') Einc
          write(recfile(14:16),'(i3.3)') int(Einc)
          open (unit=1,status='unknown',file=recfile)
          write(1,'("# ",a1," + ",i3,a2,": Recoil spectrum for ",
     +      i3,a2)') parsym(k0),Atarget,Starget,A,nuc(Z)
          write(1,'("# E-incident = ",f7.3)') Einc
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') maxenrec+1
          write(1,'("# E-out   Cross section")')
          do 30 nen=0,maxenrec
            write(1,'(f7.3,1p,e12.5)') Erec(Zcomp,Ncomp,nen),
     +        specrecoil(Zcomp,Ncomp,nen)
   30     continue
          close (unit=1)
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
