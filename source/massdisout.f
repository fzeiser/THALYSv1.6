      subroutine massdisout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 12, 2013
c | Task  : Output of fission fragment yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  type(2)
      character*90 yieldfile,nufile,fpfile
      integer      i,iz,ia,in,nen
c
c Write results to separate files
c
c yieldfile : file with fission yields
c natstring : string extension for file names
c iso       : counter for isotope
c Einc      : incident energy in MeV
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c Ztarget   : charge number of target nucleus
c yieldapre : pre-neutron emission mass yield
c yieldapost: post-neutron emission corrected mass yield
c
      type(1)='A'
      type(2)='N'
      do i=1,2
        yieldfile='yield'//type(i)//'000.000.fis'//natstring(iso)
        write(yieldfile(7:13),'(f7.3)') Einc
        write(yieldfile(7:9),'(i3.3)') int(Einc)
        open (unit=1,status='unknown',file=yieldfile)
        write(1,'("# ",a1," + ",i3,a2,": fission yields")')
     +    parsym(k0),Atarget,Starget
        write(1,'("# E-incident = ",f7.3)') Einc
        write(1,'("# ")')
        write(1,'("# ")')
        write(1,'("#  ",a1,"   Post-n yield   Pre-n yield")') type(i)
        if (i.eq.1) then
          do 10 ia=1,Atarget
            write(1,'(i3,1p,2e15.4)') ia,yieldapost(ia),yieldapre(ia)
  10      continue
        else
          do 20 in=1,Ntarget
            write(1,'(i3,1p,2e15.4)') in,yieldnpost(in),yieldnpre(in)
  20      continue
        endif
        close (unit=1)
      enddo
      nufile='nuA000.000.fis'//natstring(iso)
      write(nufile(4:10),'(f7.3)') Einc
      write(nufile(4:6),'(i3.3)') int(Einc)
      open (unit=1,status='unknown',file=nufile)
      write(1,'("# ",a1," + ",i3,a2,": Average prompt neutron ",
     +  "multiplicity as function of mass")') parsym(k0),Atarget,Starget
      write(1,'("# E-incident = ",f7.3," MeV")') Einc
      write(1,'("# Mean value (nubar-prompt) = ",f10.5)') nubar
      write(1,'("# ")')
      write(1,'("#  A    Post-n nu      Pre-n nu")')
      do 30 ia=1,Atarget
        write(1,'(i3,1p,2e15.4)') ia,nupost(ia),nupre(ia)
  30  continue
      close (unit=1)
      nufile='Pnu000.000.fis'//natstring(iso)
      write(nufile(4:10),'(f7.3)') Einc
      write(nufile(4:6),'(i3.3)') int(Einc)
      open (unit=1,status='unknown',file=nufile)
      write(1,'("# ",a1," + ",i3,a2,": Prompt neutron multiplicity ",
     +  "distribution ")') parsym(k0),Atarget,Starget
      write(1,'("# E-incident = ",f7.3," MeV")') Einc
      write(1,'("# Mean value (nubar-prompt) = ",f10.5)') nubar
      write(1,'("# ")')
      write(1,'("# number   nu ")')
      do 40 i=1,50
        if (Pdisnu(i).gt.0.) write(1,'(i3,1p,e15.4)') i,Pdisnu(i)
  40  continue
      close (unit=1)
c
c Write nubar
c
c nubarexist : flag for existence of nubar file
c numinc     : number of incident energies
c eninc      : incident energy in MeV
c
      nufile='nubar.tot'//natstring(iso)//'       '
      if (.not.nubarexist) then
        nubarexist=.true.
        open (unit=1,status='unknown',file=nufile)
        write(1,'("# ",a1," + ",i3,a2,
     +    ": Average prompt neutron multiplicity (nubar-prompt)")')
     +    parsym(k0),Atarget,Starget
        write(1,'("# ")')
        write(1,'("# # energies =",i3)') numinc
        write(1,'("# ")')
        write(1,'("# E-incident     nubar")')
        do 110 nen=1,nin-1
          write(1,'(1p,e10.3,e15.4)') eninc(nen),0.
  110   continue
      else
        open (unit=1,status='old',file=nufile)
        do 120 nen=1,nin+4
          read(1,*,end=130,err=130)
  120   continue
      endif
        write(1,'(1p,e10.3,e15.4)') Einc,nubar
  130   close (unit=1)
c
c Write ff/fp residual production
c
c fpexist    : flag for existence of fission product
c fpfile     : file with fission product
c yieldzapre : pre-neutron emission isotopic yield
c yieldzapost: post-neutron emission corrected isotopic yield
c
      do 210 ia=1,Atarget
        do 220 iz=1,Ztarget
          if (yieldzapre(iz,ia).lt.1.e-3.and..not.fpexist(iz,ia))
     +      goto 220
          fpfile='fp000000.tot'//natstring(iso)
          write(fpfile(3:8),'(2i3.3)') iz,ia
          if (.not.fpexist(iz,ia)) then
            fpexist(iz,ia)=.true.
            open (unit=1,status='unknown',file=fpfile)
            write(1,'("# ",a1," + ",i3,a2,": Fission product yield of ",
     +        i3,a2)') parsym(k0),Atarget,Starget,ia,nuc(iz)
            write(1,'("# ")')
            write(1,'("# # energies =",i3)') numinc
            write(1,'("# ")')
            write(1,'("# E-incident   Post-n xs       Pre-n xs ",
     +        "    Post-n yield    Pre-n yield")')
            do 230 nen=1,nin-1
              write(1,'(1p,e10.3,4e15.4)') eninc(nen),0.,0.,0.,0.
  230       continue
          else
            open (unit=1,status='old',file=fpfile)
            do 240 nen=1,nin+4
              read(1,*,end=250,err=250)
  240       continue
          endif
          write(1,'(1p,e10.3,4e15.4)') Einc,xsfpzapost(iz,ia),
     +      xsfpzapre(iz,ia),yieldzapost(iz,ia),yieldzapre(iz,ia)
  250     close (unit=1)
  220   continue
  210 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
