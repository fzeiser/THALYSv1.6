      subroutine directecis
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 19, 2013
c | Task  : ECIS calculation of direct cross section
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      rotational,vibrational,jlmloc
      character*13 outfile
      integer      type,Zix,Nix,A,odd,i,l
c
c ********************** Set ECIS input parameters *********************
c
c rotational      : flag for rotational input
c vibrational     : flag for vibrational input
c jlmloc          : flag for JLM OMP
c legendre        : logical for output of Legendre coefficients
c title           : title of ECIS input file
c ecis1,ecis2     : 100 input flags ('T' or 'F') for ECIS
c flagrel         : flag for relativistic kinematics
c ncoll           : number of nuclear states
c iterm           : number of iterations
c flagstate       : flag for optical model potential for each excited
c                   state
c npp             : number of optical potentials
c hint            : integration step size h
c rmatch          : matching radius
c projmass,parmass: mass of projectile
c k0              : index of incident particle
c njmax           : maximal number of j-values in ECIS
c Atarget         : mass number of target nucleus
c onethird        : 1/3
c Einc            : incident energy in MeV
c numl            : maximum l-value (set in talys.cmb)
c spin,parspin    : spin of incident particle
c tarmass         : mass of target nucleus
c prodZ           : product of charges of projectile and target nucleus
c Ztarget         : charge number of target nucleus
c parZ            : charge number of particle
c Nband           : number of vibrational bands
c angbeg          : first angle
c angend          : last angle
c
c Specific ECIS flags:
c ecis1(30)=T : DWBA
c ecis2(9)=T  : output of total, reaction, elastic and inelastic c.s.
c ecis2(14)=T : output of inelastic angular distribution
c ecis2(15)=T : output of Legendre coefficients
c
      open (unit=9,status='unknown',file='ecisdisc.inp')
      rotational=.false.
      vibrational=.true.
      jlmloc=.false.
      legendre=.true.
      title='Direct discrete cross sections by DWBA            '
      ecis1='FFFFFTFFFFFFFFFFFFFFFFFFFFFTFTFFFFFFFFFFFFFFFFFFFF'
      ecis2='FFFFFFFFTFFFFTTFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
      if (flagrel) ecis1(8:8)='T'
      ncoll=2
      iterm=1
      if (flagstate) then
        npp=ncoll
      else
        npp=1
      endif
      hint=0.
      rmatch=0.
c
c We use a simple formula to estimate the required number of j-values:
c    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
c and we always take a minimum of njmax=20.
c
      projmass=parmass(k0)
      njmax=int(2.4*1.25*(real(Atarget)**onethird)*0.22*
     +  sqrt(projmass*Einc))
      njmax=max(njmax,20)
      njmax=min(njmax,numl)
      spin=parspin(k0)
      prodZ=real(Ztarget*parZ(k0))
      Nband=1
      if (k0.eq.1) then
        angbeg=0.
      else
        angbeg=0.00001
      endif
      angend=180.
c
c ******************* Write ECIS input files ***************************
c
c flagecisinp: flag for existence of ecis input file
c parskip   : logical to skip outgoing particle
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c AA,A      : mass number of residual nucleus
c deftype   : deformation length (D) or parameter (B)
c disp      : flag for dispersive optical model
c efer      : Fermi energy
c w2disp,...: constants for imaginary potentials
c odd       : odd (1) or even (0) nucleus
c resmass   : mass of residual nucleus
c nucmass   : mass of nucleus
c idvib     : identifier for existence of vibrational state inside
c             rotational model
c tarspin   : spin of target nucleus
c tarparity : parity of target nucleus
c anginc    : angle increment
c nangle    : number of angles
c
      flagecisinp=.false.
      do 10 type=k0,k0
        if (parskip(type)) goto 10
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        A=AA(0,0,type)
        if (deftype(Zix,Nix).eq.'B') ecis1(6:6)='F'
        if (disp(Zix,Nix,k0)) then
          ecis1(10:10)='T'
          efer=ef(Zix,Nix,k0)
          w2disp=w2(Zix,Nix,k0)
          d3disp=d3(Zix,Nix,k0)
          d2disp=d2(Zix,Nix,k0)
        endif
        odd=mod(A,2)
        resmass=nucmass(Zix,Nix)
        idvib(1)=0
        idvib(2)=0
        tarspin=0.
        tarparity='+'
        anginc=180./nangle
c
c 1. Direct collective states
c
c numlev2  : maximum number of levels
c deform   : deformation parameter
c eoutdis  : outgoing energy of discrete state reaction
c Elevel   : energy of level
c edis     : energy of level
c Q        : Q-value for target nucleus
c eninccm  : center-of-mass incident energy in MeV
c parA     : mass number of particle
c Jlevel   : spin of level
c jdis     : spin of level
c Plevel   : parity of level
c cparity  : parity (character)
c parlev   : parity of level
c jcore    : spin of level of core nucleus
c pcore    : parity of level of core nucleus
c nlev,nlev: number of excited levels for nucleus
c iband    : band number of level
c Jband    : angular momentum
c Kmag     : magnetic quantum number
c iph      : phonon (1 or 2)
c vibbeta  : vibrational deformation parameter
c ecisinput: subroutine to create ECIS input file
c
        do 20 i=0,numlev2
          if (i.eq.0.and.type.eq.k0) goto 20
          if (deform(Zix,Nix,i).eq.0.) goto 20
          if (eoutdis(type,i).le.0.) goto 20
          Elevel(2)=edis(Zix,Nix,i)-Q(type)
          if (eninccm.le.Elevel(2)+0.1*parA(type)) goto 20
          if (odd.eq.0) then
            Jlevel(2)=jdis(Zix,Nix,i)
            Plevel(2)=cparity(parlev(Zix,Nix,i))
          else
            if (nlev(Zix,Nix).eq.1) then
              Jlevel(2)=2.
              Plevel(2)='+'
            else
              Jlevel(2)=jcore(Zix,Nix,i)
              Plevel(2)=cparity(pcore(Zix,Nix,i))
            endif
          endif
          iband(2)=1
          Jband(1)=int(Jlevel(2))
          Kmag(1)=0
          iph(2)=1
          vibbeta(1)=deform(Zix,Nix,i)
          flagecisinp=.true.
          call ecisinput(Zix,Nix,type,Einc,rotational,vibrational,
     +      jlmloc)
   20   continue
c
c 2. Giant resonance states
c
c nanglecont: number of angles for continuum
c betagr    : deformation parameter for giant resonance
c Egrcoll   : energy of giant resonance
c sgn       : +1 for even argument, -1 for odd argument
c
        if (type.eq.k0) then
          anginc=180./nanglecont
          do 30 l=0,3
            do 30 i=1,2
              if (betagr(l,i).eq.0.) goto 30
              Elevel(2)=Egrcoll(l,i)
              if (eninccm.le.Elevel(2)+0.1*parA(type)) goto 30
              Jlevel(2)=real(l)
              iband(2)=1
              Jband(1)=int(Jlevel(2))
              Kmag(1)=0
              iph(2)=1
              Plevel(2)=cparity(int(sgn(l)))
              vibbeta(1)=betagr(l,i)
              if (deftype(Zix,Nix).eq.'D')
     +          vibbeta(1)=vibbeta(1)*rv*real(Atarget)**onethird
              flagecisinp=.true.
              call ecisinput(Zix,Nix,type,Einc,rotational,vibrational,
     +          jlmloc)
   30     continue
        endif
   10 continue
      if (.not.flagecisinp) return
      write(9,'("fin")')
      close (unit=9)
      legendre=.false.
c
c ************ ECIS calculation for discrete levels ********************
c
c flagoutecis: flag for output of ECIS results
c outfile    : output file
c nulldev    : null device
c ecist      : subroutine ecis, adapted for TALYS
c ecisstatus : status of ECIS file
c
      if (flagoutecis) then
        outfile='ecisdisc.out '
      else
        outfile=nulldev
      endif
      call ecist('ecisdisc.inp ',outfile,'ecis.dircs   ',
     +  'ecis.dirin   ','null         ','ecis.dirang  ','ecis.dirleg  ')
      open (unit=9,status='unknown',file='ecisdisc.inp')
      close (unit=9,status=ecisstatus)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
