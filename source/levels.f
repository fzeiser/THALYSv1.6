      subroutine levels(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 9, 2013
c | Task  : Discrete levels
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*1  bas(numlev)
      character*4  levelchar
      character*90 levfile
      integer      Zix,Nix,Z,A,nlev2,ia,nlevlines,nnn,i,j,klev(numlev),
     +             ii,nb,Lis
      real         br(numlev),con(numlev)
c
c ******************** Default nuclear levels **************************
c
c For any nuclide, we first assign default ground state spins and
c parities. If there is information in the
c discrete level file, this will of course be overwritten. The index 0
c of edis, etc. represents the ground state, the index 1 the first
c excited state, etc.
c
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c ZZ,Z    : charge number of residual nucleus
c AA,A    : mass number of residual nucleus
c edis    : energy of level
c jdis    : spin of level
c parlev  : parity of level
c gsspin  : ground state spin
c gsparity: ground state parity
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      edis(Zix,Nix,0)=0.
      jdis(Zix,Nix,0)=gsspin(Zix,Nix)
      parlev(Zix,Nix,0)=gsparity(Zix,Nix)
c
c We also assign a default value to the first excited state, if it does
c not exist in the discrete level file
c
c jassign    : flag for assignment of spin
c passign    : flag for assignment of parity
c branchratio: gamma-ray branching ratio to level
c nbranch    : number of branching levels
c bassign    : flag for assignment of branching ratio
c conv       : conversion coefficient
c tau        : lifetime of state in seconds
c
      edis(Zix,Nix,1)=min(26./A,10.)
      jdis(Zix,Nix,1)=jdis(Zix,Nix,0)+2.
      parlev(Zix,Nix,1)=parlev(Zix,Nix,0)
      jassign(Zix,Nix,1)='J'
      passign(Zix,Nix,1)='P'
      branchratio(Zix,Nix,1,1)=1.
      nbranch(Zix,Nix,1)=1
      bassign(Zix,Nix,1,0)='B'
      conv(Zix,Nix,1,0)=0.
      tau(Zix,Nix,1)=0.
c
c ************************ Read nuclear levels *************************
c
c Note that nlev is the number of excited levels, i.e. excluding the
c ground state.
c
c 1. Inquire whether file is present
c
c levelfile: discrete level file
c levelchar: help variable
c disctable: table with discrete levels
c levfile  : level file
c path     : directory containing structure files to be read
c lenpath  : length of pathname
c
      if (levelfile(Zix)(1:1).ne.' ') then
        levfile=levelfile(Zix)
      else
        levelchar='z   '
        write(levelchar(2:4),'(i3.3)') Z
        if (disctable.eq.1)
     +    levfile=path(1:lenpath)//'levels/final/'//levelchar
        if (disctable.eq.2)
     +    levfile=path(1:lenpath)//'levels/exp/'//levelchar
        if (disctable.eq.3)
     +    levfile=path(1:lenpath)//'levels/hfb/'//levelchar
        inquire (file=levfile,exist=lexist)
        if (.not.lexist) return
      endif
      open (unit=2,status='unknown',file=levfile)
c
c 2. Search for the isotope under consideration
c
c nlev,nlev2 : number of excited levels for nucleus
c ia         : mass number from level file
c nlevlines  : number of lines on discrete level file for nucleus
c nnn        : number of levels in discrete level file
c flagautorot: flag for automatic rotational coupled channels
c
      do j=1,numlev
        con(j)=0.
        br(j)=0.
        bas(j)='B'
      enddo
      nlev2=0
      nnn=0
   10 read(2,'(4x,i4,2i5)',end=100) ia,nlevlines,nnn
      if (A.ne.ia) then
        do 20 i=1,nlevlines
          read(2,'()')
   20   continue
        goto 10
      endif
      nlev2=min(nnn,numlev)
      if (nnn.lt.2) flagautorot=.false.
c
c 3. Read discrete level information
c
c branchlevel: level to which branching takes place
c ENSDF      : string from original ENSDF discrete level file
c
      do 30 i=0,nlev2
        read(2,'(4x,f11.6,f6.1,3x,i2,i3,18x,e10.3,1x,2a1,a18)')
     +    edis(Zix,Nix,i),jdis(Zix,Nix,i),parlev(Zix,Nix,i),
     +    nb,tau(Zix,Nix,i),jassign(Zix,Nix,i),
     +    passign(Zix,Nix,i),ENSDF(Zix,Nix,i)
        do 40 j=1,nb
          read(2,'(29x,i3,f10.6,e10.3,5x,a1)') klev(j),br(j),con(j),
     +      bas(j)
   40   continue
c
c Branching ratio from input
c
        if (nbranch(Zix,Nix,i).gt.0) then
          do 50 j=1,nbranch(Zix,Nix,i)
            conv(Zix,Nix,i,j)=con(j)
            bassign(Zix,Nix,i,j)=bas(j)
   50     continue
        else
c
c Branching ratio from discrete level file
c
          ii=0
          do 60 j=1,nb
            if (br(j).ne.0.) then
              ii=ii+1
              branchlevel(Zix,Nix,i,ii)=klev(j)
              branchratio(Zix,Nix,i,ii)=br(j)
              conv(Zix,Nix,i,ii)=con(j)
              bassign(Zix,Nix,i,ii)=bas(j)
            endif
   60     continue
          nbranch(Zix,Nix,i)=ii
        endif
c
c Spins beyond numJ are set to numJ
c
c numJ: maximal J-value
c
        jdis(Zix,Nix,i)=min(jdis(Zix,Nix,i),real(numJ))
c
c Overwrite value of isomer for shorter-lived target level.
c
c Ltarget: excited level of target
c isomer : definition of isomer in seconds
c
        if (Ltarget.ne.0.and.Zix.eq.parZ(k0).and.Nix.eq.parN(k0).
     +    and.i.eq.Ltarget.and.tau(Zix,Nix,i).lt.isomer)
     +    isomer=tau(Zix,Nix,i)
c
c Set highest discrete level equal to isomer if that exists
c
        if (i.gt.nlev(Zix,Nix).and.tau(Zix,Nix,i).ge.isomer )
     +    nlev(Zix,Nix)=i
   30 continue
c
c Lifetimes below the isomeric definition are set to zero.
c The isomeric number is determined.
c
      do 70 i=0,nlev2
        if (tau(Zix,Nix,i).lt.isomer) tau(Zix,Nix,i)=0.
   70 continue
      if (massmodel.le.1) then
        if (jassign(Zix,Nix,0).eq.'J'.and.passign(Zix,Nix,0).eq.'P')
     +    then
          jdis(Zix,Nix,0)=gsspin(Zix,Nix)
          parlev(Zix,Nix,0)=gsparity(Zix,Nix)
        endif
        if (nlev2.eq.0) then
          jdis(Zix,Nix,1)=jdis(Zix,Nix,0)+2.
          parlev(Zix,Nix,1)=parlev(Zix,Nix,0)
        endif
      endif
c
c Read extra levels which are used only for the level density matching
c problem or for direct reactions (deformation parameters). The
c branching ratios and lifetimes are not read for these higher levels.
c
c nlevmax2,numlev2: maximum number of levels
c
      nlevmax2(Zix,Nix)=min(nnn,numlev2)
      do 80 i=nlev2+1,nlevmax2(Zix,Nix)
        read(2,'(4x,f11.6,f6.1,3x,i2,i3,29x,2a1)') edis(Zix,Nix,i),
     +    jdis(Zix,Nix,i),parlev(Zix,Nix,i),nb,jassign(Zix,Nix,i),
     +    passign(Zix,Nix,i)
        do 90 j=1,nb
          read(2,*)
   90   continue
   80 continue
  100 close (unit=2)
      nlev(Zix,Nix)=min(nlev(Zix,Nix),nlev2)
      nlev(Zix,Nix)=max(nlev(Zix,Nix),1)
c
c The maximal value of Ntop is always given by the last discrete
c level of the discrete level file.
c
c Ntop: highest discrete level for temperature matching
c
      nlevmax2(Zix,Nix)=min(nnn,numlev2)
      nlevmax2(Zix,Nix)=max(nlevmax2(Zix,Nix),1)
      Ntop(Zix,Nix,0)=min(nlevmax2(Zix,Nix),Ntop(Zix,Nix,0))
c
c Check existence of excited level of target
c
c parZ   : charge number of particle
c k0     : index for incident particle
c parN   : neutron number of particle
c
      if (Ltarget.ne.0) then
        if (Zix.eq.parZ(k0).and.Nix.eq.parN(k0).and.
     +    Ltarget.gt.nlev(Zix,Nix)) then
          write(*,'(" TALYS-error: excited level of target does",
     +      " not exist")')
          stop
        endif
      endif
c
c Determine isomeric level number
c
c Lisomer: level number of isomer
c Nisomer: number of isomers for this nuclide
c
      Lis=0
      do 110 i=1,nlev(Zix,Nix)
        if (Lis.lt.numisom.and.tau(Zix,Nix,i).ne.0.) then
          Lis=Lis+1
          Lisomer(Zix,Nix,Lis)=i
        endif
  110 continue
      Nisomer(Zix,Nix)=Lis

c Determine isomeric number of target
c
c Lisoinp: user assignment of target isomer number
c Liso   : isomeric number of target
c
      if (Lisoinp.eq.-1) then
        if (Zix.eq.parZ(k0).and.Nix.eq.parN(k0)) then
          Liso=0
          if (Ltarget.ne.0) then
            do 120 i=1,Lis
              if (Ltarget.eq.Lisomer(Zix,Nix,i)) Liso=i
  120       continue
          endif
        endif
      else
        Liso=Lisoinp
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
