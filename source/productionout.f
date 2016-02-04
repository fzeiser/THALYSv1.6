      subroutine productionout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 10, 2012
c | Task  : Output of particle production cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*13 totfile
      character*15 fisfile
      integer      type,nen
c
c ************** Total particle production cross sections **************
c
c parskip     : logical to skip outgoing particle
c parname     : name of particle
c xsparticle  : total particle production cross section
c multiplicity: particle multiplicity
c
      write(*,'(/" 3. Total particle production cross sections"/)')
      do 10 type=0,6
        if (parskip(type)) goto 10
        write(*,'(1x,a8,"=",1p,e12.5,"    Multiplicity=",e12.5)')
     +    parname(type),xsparticle(type),multiplicity(type)
c
c Write results to separate file
c
c filetotal : flag for total cross sections on separate file
c natstring : string extension for file names
c iso       : counter for isotope
c numinclow : number of incident energies below Elow
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c Ztarget   : charge number of target nucleus
c numinc    : number of incident energies
c eninc,Einc: incident energy in MeV
c
        if (filetotal) then
          totfile=' prod.tot'//natstring(iso)
          write(totfile(1:1),'(a1)') parsym(type)
          if (nin.eq.numinclow+1) then
            open (unit=1,status='unknown',file=totfile)
            write(1,'("# ",a1," + ",i3,a2," Total ",a8," production")')
     +        parsym(k0),Atarget,Starget,parname(type)
            write(1,'("# Q-value    =",1p,e12.5)') Q(type)
            write(1,'("# ")')
            write(1,'("# # energies =",i3)') numinc
            write(1,'("#    E         xs         Yield")')
            do 20 nen=1,numinclow
              write(1,'(1p,e10.3,2e12.5)') eninc(nen),0.,0.
   20       continue
          else
            open (unit=1,status='old',file=totfile)
            do 30 nen=1,nin+4
              read(1,*,end=40,err=40)
   30       continue
          endif
          write(1,'(1p,e10.3,2e12.5)') Einc,xsparticle(type),
     +      multiplicity(type)
   40     close (unit=1)
        endif
   10 continue
c
c Total fission cross section
c
c flagfission : flag for fission
c xsfistot    : total fission cross section
c
      if (flagfission) then
        write(*,'(" fission =",1p,e12.5)') xsfistot
c
c Write results to separate file
c
c filefission: flag for fission cross sections on separate file
c
        if (filefission) then
          fisfile='fission.tot'//natstring(iso)
          if (nin.eq.numinclow+1) then
            open (unit=1,status='unknown',file=fisfile)
            write(1,'("# ",a1," + ",i3,a2,"   : (",a1,",f)        ",
     +        "  Total")') parsym(k0),Atarget,Starget,parsym(k0)
            write(1,'("# ")')
            write(1,'("# ")')
            write(1,'("# # energies =",i3)') numinc
            write(1,'("#    E         xs")')
            do 110 nen=1,numinclow
              write(1,'(1p,e10.3,e12.5)') eninc(nen),0.
  110       continue
          else
            open (unit=1,status='old',file=fisfile)
            do 120 nen=1,nin+4
              read(1,*,end=130,err=130)
  120       continue
          endif
          write(1,'(1p,e10.3,e12.5)') Einc,xsfistot
  130     close (unit=1)
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
