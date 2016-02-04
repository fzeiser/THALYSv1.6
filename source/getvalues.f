      subroutine getvalues(class,word,Zix,Nix,type,
     +  ibar,irad,lval,igr,val,ival,cval,flagassign)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 4, 2013
c | Task  : Assign values to keywords
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      numEkey
      parameter    (numEkey=65)
      logical      flagassign,lexist
      character*1  ch
      character*80 keyword(numEkey),word(40),key,adfile,cval
      integer      class,Zix,Nix,type,ibar,irad,lval,igr,ival,i,iz,ia,
     +             k,type2,j
      real         val,Ea,Eb,Em,D,Eadj(0:numenadj),Dadj(numenadj)
c
c Several keywords may be altered over local energy ranges
c
      data (keyword(i),i=1,numEkey) /
     +  'adepthcor', 'aradialcor',
     +  'avadjust', 'avdadjust', 'avsoadjust', 'awadjust', 'awdadjust',
     +  'awsoadjust', 'cbreak', 'cknock', 'cstrip', 'ctable',
     +  'd1adjust', 'd2adjust', 'd3adjust', 'egr', 'egradjust',
     +  'epr', 'epradjust', 'etable',
     +  'fisbar', 'fisbaradjust', 'fishw', 'fishwadjust', 'ftable',
     +  'ggr', 'ggradjust', 'gnorm', 'gpr', 'gpradjust', 'krotconstant',
     +  'lv1adjust', 'lvadjust', 'lvsoadjust', 'lw1adjust',
     +  'lwadjust', 'lwsoadjust', 'm2constant', 'ptable', 'rcadjust',
     +  'rspincut', 'rvadjust', 'rvdadjust', 'rvsoadjust', 'rwadjust',
     +  'rwdadjust', 'rwsoadjust', 's2adjust', 'sgr', 'sgradjust',
     +  'spr', 'spradjust', 'tljadjust', 'v1adjust', 'v2adjust',
     +  'v3adjust', 'v4adjust', 'vso1adjust', 'vso2adjust', 'w1adjust',
     +  'w2adjust', 'w3adjust', 'w4adjust', 'wso1adjust', 'wso2adjust'/
c
c ************************ Read values for keywords ********************
c
c Each keyword is characterized by a certain order of parameter
c and value input. They are distinguished by different classes.
c
c Classes:
c
c  1: keyword Z A real-value [optional: local adjustment]
c  2: keyword Z A integer-value
c  3: keyword Z A real-value barrier [optional: local adjustment]
c  4: keyword Z A integer-value barrier [optional: local adjustment]
c  5: keyword Z A real-value rad-type l-val [optional: local adjustment]
c  6: keyword particle-type real-value [optional: local adjustment]
c  7: keyword particle-type integer-value
c  8: keyword particle-type real value integer-value
c     [optional: local adjustment]
c  9: keyword value [optional: local adjustment]
c 10: keyword Z filename
c
c flagassign: flag to assign value or not
c word      : words on input line
c key       : keyword
c val       : real value
c ival      : integer value
c cval      : character value
c
      flagassign=.false.
      key=word(1)
      ch=word(2)(1:1)
      val=0.
      ival=0
      cval='                                                           '
c
c Z,A dependent keywords
c
c class: input class
c iz   : charge number
c ia   : mass number
c Zix  : charge number index for residual nucleus
c numZ : maximal number of protons away from the initial
c        compound nucleus
c Nix  : neutron number index for residual nucleus
c numN : maximal number of neutrons away from the initial
c        compound nucleus
c Zinit: charge number of initial compound nucleus
c Ninit: neutron number of initial compound nucleus
c
      if (class.le.5) then
        read(word(2),*,end=100,err=100) iz
        read(word(3),*,end=100,err=100) ia
        if (class.eq.2.or.class.eq.4) then
          read(word(4),*,end=100,err=100) ival
        else
          read(word(4),*,end=100,err=100) val
        endif
        Zix=Zinit-iz
        Nix=Ninit-ia+iz
        if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) goto 200
        i=4
c
c Z,A dependent keywords with possible fission barriers
c
c ibar  : fission barrier
c numbar: number of fission barriers
c
        if (class.eq.3.or.class.eq.4) then
          read(word(5),*,end=10,err=40) ibar
          if (ibar.lt.0.or.ibar.gt.numbar) goto 110
          i=5
        endif
c
c Z,A dependent keywords with gamma parameters
c
c irad  : variable to indicate M(=0) or E(=1) radiation
c lval  : multipolarity
c numgam: maximum number of l-values for gamma multipolarity
c igr   : giant resonance
c
        if (class.eq.5) then
          ch=word(5)(1:1)
          if (ch.eq.'m'.or.ch.eq.'e') then
            if (ch.eq.'m') irad=0
            if (ch.eq.'e') irad=1
            read(word(5)(2:2),*,end=100,err=100) lval
            if (lval.lt.1.or.lval.gt.numgam) goto 120
            read(word(6),*,end=10,err=100) igr
            if (igr.lt.1.or.igr.gt.2) goto 130
            i=6
          else
            read(word(5),*,end=100,err=100) lval
            if (lval.lt.0.or.lval.gt.numl) goto 120
            i=5
          endif
        endif
      endif
c
c Particle type dependent keywords
c
c parsym: symbol of particle
c
      if (class.ge.6.and.class.le.8) then
        do 20 type2=0,6
          if (ch.eq.parsym(type2)) then
            type=type2
            goto 30
          endif
   20   continue
        goto 180
   30   if (class.eq.7) then
          read(word(3),*,end=100,err=100) ival
        else
          read(word(3),*,end=100,err=100) val
        endif
        i=3
        if (class.eq.8) then
          read(word(4),*,end=10,err=40) lval
          if (lval.lt.0.or.lval.gt.numl) goto 120
          i=4
        endif
      endif
c
c Simple dependent keywords
c
      if (class.eq.9) then
        read(word(2),*,end=100,err=100) val
        i=2
      endif
c
c Keywords with input files
c
      if (class.eq.10) then
        read(word(2),*,end=100,err=100) iz
        Zix=Zinit-iz
        if (Zix.lt.0.or.Zix.gt.numZ) then
          goto 200
        else
          cval=word(3)
        endif
        i=3
      endif
c
c Local energy-dependent adjustment of parameters
c
c numEkey   : number of keywords
c Ea        : start energy of local adjustment
c Eb        : end energy of local adjustment
c Em        : intermediate energy of local adjustment
c D         : depth of local adjustment
c Nadjust   : number of adjustable parameters
c adjustkey : keyword for local adjustment
c adjustfile: file for local adjustment
c adjustix  : local adjustment index
c adjustpar : local adjustment parameters
c nenadjust : number of tabulated energies of local adjustment
c Eadjust   : tabulated energy of local adjustment
c Dadjust   : tabulated depth of local adjustment
c
   40 Ea=0.
      Eb=0.
      Em=0.
      D=0.
      adfile='                                                         '
      do 50 k=1,numEkey
        if (key.eq.keyword(k)) goto 60
   50 continue
      goto 10
   60 Eadj(0)=0.
      k=1
      if ((word(i+1)(1:1).ge.'a'.and.word(i+1)(1:1).le.'z').or.
     +  (word(i+1)(1:1).ge.'A'.and.word(i+1)(1:1).le.'Z')) then
        read(word(i+1),*,end=10,err=100) adfile
        inquire (file=adfile,exist=lexist)
        if (.not.lexist) goto 140
        open (unit=1,status='unknown',file=adfile)
   70   read(1,*,end=80,err=150) Eadj(k),Dadj(k)
        if (Eadj(k).eq.Eadj(k-1)) goto 70
        if (Eadj(k).lt.Eadj(k-1)) goto 160
        k=k+1
        if (k.gt.numenadj) goto 170
        goto 70
   80   close (1)
      else
        read(word(i+1),*,end=10,err=10) Ea
        read(word(i+2),*,end=10,err=10) Eb
        read(word(i+3),*,end=10,err=10) Em
        read(word(i+4),*,end=10,err=10) D
      endif
      if (Nadjust.lt.numadj) Nadjust=Nadjust+1
      nenadjust(Nadjust)=k-1
      do 90 j=1,nenadjust(Nadjust)
        Eadjust(Nadjust,j)=Eadj(j)
        Dadjust(Nadjust,j)=Dadj(j)
   90 continue
      adjustkey(Nadjust)=key
      adjustfile(Nadjust)=adfile
      adjustix(Nadjust,1)=Zix
      adjustix(Nadjust,2)=Nix
      adjustix(Nadjust,3)=type
      if (class.eq.5) then
        adjustix(Nadjust,4)=lval
      else
        adjustix(Nadjust,4)=ibar
      endif
      adjustpar(Nadjust,1)=Ea
      adjustpar(Nadjust,2)=Eb
      adjustpar(Nadjust,3)=Em
      adjustpar(Nadjust,4)=D
   10 flagassign=.true.
      return
c
c Error and warning messages
c
  100 write(*,'(" TALYS-error: Wrong input for: ",a80)') key
      stop
  110 write(*,'(" TALYS-error: 0(1) <= fission barrier <=",i3,
     +  ", ibar index out of range: ",a80)') numbar,key
      stop
  120 write(*,'(" TALYS-error: 0 <= multipole radiation <= ",i1,
     +  ", lval index out of range: ",a80)') numgam,key
      stop
  130 write(*,'(" TALYS-error: 0 <= resonance number <= 2",
     +  ", igr index out of range: ",a80)') key
      stop
  140 write(*,'(" TALYS-error: parameter file ",a80,
     +  " does not exist for keyword ",a80)') adfile,key
      stop
  150 write(*,'(" TALYS-error: parameter file ",a80,
     +  " has wrong format for keyword ",a80)') adfile,key
      stop
  160 write(*,'(" TALYS-error: parameter file ",a80,
     +  " must have energies in increasing order for keyword ",a80)')
     +  adfile,key
      stop
  170 write(*,'(" TALYS-error: parameter file ",a80,
     +  " has more than ",i6," energies for keyword ",a80)')
     +  adfile,numenadj,key
      stop
  180 write(*,'(" TALYS-error: wrong particle symbol for: ",a80)') key
      stop
  200 write(*,'(" TALYS-warning: Z,N index out of range,",
     +  " keyword ignored: ",a80)') key
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
