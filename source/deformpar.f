      subroutine deformpar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 14, 2013
c | Task  : Deformation parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist,first2,first3,first4
      character*1  colltype1,deftype1,leveltype1
      character*4  defchar
      character*90 deffile
      integer      Zix,Nix,Z,N,A,ia,ndisc,i,natpar,k,iirot,idef,irot,ii,
     +             nex,distance,k2,odd,vibband1,lband1,Kmag1,iphonon1,
     +             nrotlev,type,ibar
      real         deform1(numrotcc),dspin,R
c
c ************************ Read deformation parameters *****************
c
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c ZZ        : charge number of residual nucleus
c NN        : neutron number of residual nucleus
c AA        : mass number of residual nucleus
c colltype  : type of collectivity (S, V, R or A)
c disctable : table with discrete levels
c defchar   : help variable
c deffile   : deformation parameter file
c deformfile: deformation parameter file
c path      : directory containing structure files to be read
c lenpath   : length of pathname
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      colltype(Zix,Nix)='S'
      if (disctable.eq.3) goto 150
      defchar='z   '
      write(defchar(2:4),'(i3.3)') Z
c
c 1. Inquire whether file is present
c
      if (deformfile(Zix)(1:1).ne.' ') then
        deffile=deformfile(Zix)
      else
        deffile=path(1:lenpath)//'deformation/'//defchar
      endif
      inquire (file=deffile,exist=lexist)
      if (.not.lexist) goto 150
      open (unit=2,status='old',file=deffile)
c
c 2. Search for the isotope under consideration
c
c ia          : mass number from deformation parameter file
c ndisc       : number of lines on discrete level file for nucleus
c colltype1,..: help variables
c
   10 read(2,'(4x,2i4,2(3x,a1))',end=140) ia,ndisc,colltype1,deftype1
      if (A.ne.ia) then
        do 20 i=1,ndisc
          read(2,'()')
   20   continue
        goto 10
      endif
c
c Initialization
c
c numrotcc: number of rotational deformation parameters
c deform1 : deformation parameter
c
      do 30 i=1,numrotcc
        deform1(i)=0.
   30 continue
c
c 3. Read deformation parameters
c
c flagspher : flag to force spherical optical model
c iirot,....: help variables
c deftype   : deformation length (D) or parameter (B)
c nex       : level number
c numlev2   : maximum number of levels
c indexlevel: level index
c maxrot    : number of included excited rotational levels
c leveltype : type of level (rotational (R) or vibrational (V))
c vibband   : band number of level
c maxband   : highest vibrational level added to rotational model
c iphonon   : phonon (1 or 2)
c Kmag,Kband: magnetic quantum number
c defpar    : deformation parameter
c lband     : angular momentum
c nrot      : number of deformation parameters for rotational nucleus
c rotpar    : deformation parameters for rotational nucleus
c
c The coupling scheme is read from the nuclear structure database.
c Different actions are performed if the requested calculation is for a
c spherical (S), vibrational (V) or a rotational (R) nucleus.
c
      if (colltype1.ne.'R'.and.colltype1.ne.'A'.and.colltype1.ne.'V')
     +  colltype1='S'
      if (flagspher.and.colltype1.eq.'A') then
        iirot=1
      else
        iirot=numrotcc
      endif
      if (flagspher) colltype1='S'
      colltype(Zix,Nix)=colltype1
      deftype(Zix,Nix)=deftype1
      idef=0
      irot=0
      ii=0
      do 110 i=1,ndisc
        read(2,'(i4,3x,a1,4i4,4f9.5)') nex,leveltype1,vibband1,lband1,
     +    Kmag1,iphonon1,(deform1(k),k=1,iirot)
        if (nex.le.numlev2) then
          idef=idef+1
          indexlevel(Zix,Nix,i)=nex
          if (leveltype1.eq.'R') then
            irot=irot+1
            if (irot.gt.maxrot+1) leveltype1='D'
          endif
          leveltype(Zix,Nix,nex)=leveltype1
          vibband(Zix,Nix,i)=vibband1
          if (colltype(Zix,Nix).eq.'R'.and.vibband1.gt.maxband)
     +      leveltype1='D'
          iphonon(Zix,Nix,i)=max(iphonon1,1)
          if (leveltype1.eq.'R'.and.rotpar(Zix,Nix,1).eq.0.) then
c
c Default: read deformation parameters from mass table
c
            if (nex.eq.0) then
              if (deform1(1).eq.0.) then
                deform1(1)=beta2(Zix,Nix,0)
                deform1(2)=beta4(Zix,Nix)
              endif
            endif
            do 120 k=1,numrotcc
              rotpar(Zix,Nix,k)=deform1(k)
              if (deform1(k).eq.0.) then
                nrot(Zix,Nix)=k-1
                goto 130
              endif
  120       continue
          endif
          if (leveltype1.eq.'V'.and.defpar(Zix,Nix,vibband1).eq.0.) then
            lband(Zix,Nix,vibband1)=lband1
            Kband(Zix,Nix,vibband1)=Kmag1
            defpar(Zix,Nix,vibband1)=deform1(1)
          endif
  130     if (colltype(Zix,Nix).eq.'S'.or.leveltype1.eq.'D') then
c
c For DWBA, we only include natural parity states.
c
c natpar : natural parity
c jdis   : spin of level
c parlev : parity of level
c deform : deformation parameter
c indexcc: level index for coupled channel
c ndef   : number of collective levels
c
            natpar=sgn(int(jdis(Zix,Nix,nex)))
            if (parlev(Zix,Nix,nex).eq.natpar)
     +        deform(Zix,Nix,nex)=deform1(1)
            if (i.gt.1.and.i.le.numrotcc+1.and.leveltype1.eq.'R')
     +        deform(Zix,Nix,nex)=rotpar(Zix,Nix,i-1)
          else
            ii=ii+1
            indexcc(Zix,Nix,ii)=nex
          endif
        endif
  110 continue
      ndef(Zix,Nix)=idef
  140 close (unit=2)
c
c ******************** Default deformation parameters ******************
c
c distance: number of nucleons to closest magic number
c odd     : odd (1) or even (0) nucleus
c
c Automatic assignment of rotational deformation parameters.
c Calculate distance to closest magic number as a measure for
c deformation.
c
  150 distance=1000
      do 160 k2=1,8
        distance=min(abs(N-magic(k2)),distance)
        distance=min(abs(Z-magic(k2)),distance)
  160 continue
      odd=mod(A,2)
c
c Read rotational deformation parameters
c
c flagautorot  : flag for automatic rotational coupled channels
c                calculations for A > 150
c deformmodel  : model for theoretical deformation parameters
c beta2,beta4  : deformation parameters
c dspin        : angular momentum increase for rotational band
c
      if (colltype(Zix,Nix).eq.'S'.and.(.not.flagspher).and.A.gt.150.
     +  and.distance.ge.8.and.flagautorot) then
        indexlevel(Zix,Nix,1)=0
        indexcc(Zix,Nix,1)=0
        leveltype(Zix,Nix,0)='R'
        if (odd.eq.0) then
          dspin=2.
        else
          dspin=1.
        endif
        ndef(Zix,Nix)=maxrot+1
        do 180 i=2,ndef(Zix,Nix)
          ii=i-1
          spin=jdis(Zix,Nix,0)+dspin*ii
          do 190 nex=1,numlev2
            if (spin.eq.jdis(Zix,Nix,nex).and.
     +        parlev(Zix,Nix,nex).eq.parlev(Zix,Nix,0)) then
              indexlevel(Zix,Nix,i)=nex
              indexcc(Zix,Nix,i)=nex
              leveltype(Zix,Nix,nex)='R'
              goto 180
            endif
  190     continue
  180   continue
        if (indexcc(Zix,Nix,ndef(Zix,Nix)).ne.0) colltype(Zix,Nix)='R'
        nrot(Zix,Nix)=2
        rotpar(Zix,Nix,1)=beta2(Zix,Nix,0)
        rotpar(Zix,Nix,2)=beta4(Zix,Nix)
        deftype(Zix,Nix)='B'
  200   close (unit=2)
        if (odd.eq.0) goto 400
      endif
c
c Number of rotational levels should remain below maxrot
c
c nrotlev: number of rotational levels
c
      nrotlev=0
      do 310 nex=1,numlev2
        if (leveltype(Zix,Nix,nex).eq.'R') then
          nrotlev=nrotlev+1
          if (nrotlev.gt.maxrot) leveltype(Zix,Nix,nex)='D'
        endif
  310 continue
c
c Assign vibrational deformation parameters
c
c Systematics for first 2+, 3- and 4+ vibrational states, derived
c from individual deformation parameters. Also, we assign small
c deformation parameter to all remaining discrete levels.
c
c first2-4: flag to determine first state of specific spin
c k0      : index of incident particle
c edis    : energy of level
c
      if (odd.ne.0) goto 400
      if (colltype(Zix,Nix).ne.'S') then
        first2=.false.
        first3=.false.
        first4=.false.
      else
        first2=.true.
        first3=.true.
        first4=.true.
      endif
      type=2*Zix+Nix
      do 320 k=0,numlev
        if (k.eq.0.and.type.eq.k0) goto 320
        if (colltype(Zix,Nix).ne.'S'.and.leveltype(Zix,Nix,k).eq.'V')
     +    goto 320
        if (colltype(Zix,Nix).eq.'A'.and.leveltype(Zix,Nix,k).eq.'R')
     +    goto 320
        if (colltype(Zix,Nix).eq.'R'.and.leveltype(Zix,Nix,k).eq.'R')
     +    goto 320
        if (jdis(Zix,Nix,k).eq.0.) goto 320
        if (first2.and.jdis(Zix,Nix,k).eq.2..and.
     +      parlev(Zix,Nix,k).eq.1) then
          if (leveltype(Zix,Nix,k).ne.'R'.and.deform(Zix,Nix,k).eq.0.)
     +      then
            deform(Zix,Nix,k)=0.40*exp(-0.012*A)+
     +        0.025*min(distance,5)
            if (edis(Zix,Nix,k).le.0.1) deform(Zix,Nix,k)=0.02
            if (deftype(Zix,Nix).eq.'D')
     +        deform(Zix,Nix,k)=deform(Zix,Nix,k)*1.24*(A**onethird)
          endif
          first2=.false.
          goto 320
        endif
        if (first3.and.jdis(Zix,Nix,k).eq.3..and.
     +      parlev(Zix,Nix,k).eq.-1) then
          if (leveltype(Zix,Nix,k).ne.'R'.and.deform(Zix,Nix,k).eq.0.)
     +      then
            deform(Zix,Nix,k)=0.35*exp(-0.008*A)
            if (edis(Zix,Nix,k).le.0.1) deform(Zix,Nix,k)=0.02
            if (deftype(Zix,Nix).eq.'D')
     +        deform(Zix,Nix,k)=deform(Zix,Nix,k)*1.24*(A**onethird)
          endif
          first3=.false.
          goto 320
        endif
        if (first4.and.jdis(Zix,Nix,k).eq.4..and.
     +      parlev(Zix,Nix,k).eq.1) then
          if (leveltype(Zix,Nix,k).ne.'R'.and.deform(Zix,Nix,k).eq.0.)
     +      then
            deform(Zix,Nix,k)=0.20*exp(-0.006*A)
            if (edis(Zix,Nix,k).le.0.1) deform(Zix,Nix,k)=0.02
            if (deftype(Zix,Nix).eq.'D')
     +        deform(Zix,Nix,k)=deform(Zix,Nix,k)*1.24*(A**onethird)
          endif
          first4=.false.
          goto 320
        endif
        if (deform(Zix,Nix,k).ne.0.) goto 320
        natpar=sgn(int(jdis(Zix,Nix,k)))
        if (parlev(Zix,Nix,k).eq.natpar) deform(Zix,Nix,k)=0.02
        if (deftype(Zix,Nix).eq.'D')
     +    deform(Zix,Nix,k)=deform(Zix,Nix,k)*1.24*(A**onethird)
  320 continue
c
c ************** Rigid body value for moment of inertia ****************
c
c R       : radius
c onethird: 1/3
c Irigid0 : undeformed rigid body value of moment of inertia
c parmass : mass of particle in a.m.u.
c amu     : atomic mass unit in MeV
c hbarc   : hbar.c in MeV.fm
c Irigid  : rigid body value of moment of inertia
c
  400 R=1.2*A**onethird
      Irigid0(Zix,Nix)=0.4*R*R*A*parmass(1)*amu/(hbarc**2)
      do 410 ibar=0,numbar
        Irigid(Zix,Nix,ibar)=(1+abs(beta2(Zix,Nix,ibar))/3.)*
     +    Irigid0(Zix,Nix)
  410 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
