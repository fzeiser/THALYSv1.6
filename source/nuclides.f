      subroutine nuclides
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 17, 2013
c | Task  : Properties of nuclides
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,type,Zix,Nix,A,odd,NL
c
c ************************ Properties of nuclides **********************
c
c Assignment of Z and N of all possible residual nuclei. The initial
c compound nucleus (created by projectile + target) has the indices
c (0,0). The first index represents the number of protons and
c the second index the number of neutrons away from the initial
c compound nucleus. Example: For the reaction p + 208Pb, the set
c (0,0) represents 209Bi and the set (1,2) represents 206Pb. Many arrays
c have Zindex and Nindex as their first two indices. At any point
c in the reaction calculation, given Zcomp, Ncomp and the particle type,
c these variables are directly known through the arrays we initialize
c here. ZZ is the charge number of the nucleus that is reached through
c emission of ejectile type from nucleus (Zcomp,Ncomp).
c Note that ZZ, NN, AA, Zinit and Ninit represent true values of the
c charge and neutron number and that Zindex, Nindex, Zcomp and Ncomp are
c indices relative to the initial compound nucleus.
c Example: 56Fe(n,p)56Mn reaction:
c
c Zcomp=0 (primary compound nucleus)
c Ncomp=0 (primary compound nucleus)
c Zindex=1
c Nindex=0
c Zinit=26
c Ninit=31
c ZZ=25
c NN=31
c
c Zcomp    : charge number index for compound nucleus
c maxZ,numZ: maximal number of protons away from initial compound
c            nucleus
c Ncomp    : neutron number index for compound nucleus
c maxN,numN: maximal number of neutrons away from initial compound
c            nucleus
c Zinit    : charge number of initial compound nucleus
c Ninit    : neutron number of initial compound nucleus
c Zindex   : charge number index for residual nucleus
c parZ     : charge number of particle
c Nindex   : neutron number index for residual nucleus
c parN     : neutron number of particle
c ZZ       : charge number of residual nucleus
c NN       : neutron number of residual nucleus
c AA       : mass number of residual nucleus
c
      do 10 Zcomp=0,numZ
        do 10 Ncomp=0,numN
          do 10 type=0,6
            Zindex(Zcomp,Ncomp,type)=Zcomp+parZ(type)
            Nindex(Zcomp,Ncomp,type)=Ncomp+parN(type)
            ZZ(Zcomp,Ncomp,type)=Zinit-Zindex(Zcomp,Ncomp,type)
            NN(Zcomp,Ncomp,type)=Ninit-Nindex(Zcomp,Ncomp,type)
            AA(Zcomp,Ncomp,type)=ZZ(Zcomp,Ncomp,type)+
     +        NN(Zcomp,Ncomp,type)
   10 continue
c
c We make sure the lightest possible residual nucleus is heavier than
c an alpha particle. Also, we set that at the beginning no structure
c information is available at all .
c
c strucexist: flag to determine whether structure info for this nucleus
c             already exists
c strucwrite: flag for output of nuclear structure info
c invexist  : logical to determine necessity of new inverse cross
c             section and transmission coefficients calculation
c
      maxZ=min(maxZ,Zinit-3)
      maxN=min(maxN,Ninit-3)
      do 20 Zix=0,numZ
        do 20 Nix=0,numN
          strucexist(Zix,Nix)=.false.
          strucwrite(Zix,Nix)=.false.
          invexist(Zix,Nix)=.false.
   20 continue
c
c Set the maximum number of discrete levels for each residual nucleus
c
c nlev      : number of excited levels for nucleus
c nlevbin   : number of excited levels for binary nucleus
c k0        : index of incident particle
c nlevmaxres: maximum number of included discrete levels for residual
c             nucleus
c nlevmax   : maximum number of included discrete levels for target
c
c Note the order of priority: nlevmax (keyword: maxlevels) overrules
c the value set for nlevbin (keyword: maxlevelsbin) for the inelastic
c channel. nlevbin overrules nlevmaxres (keyword: maxlevelsres) for
c particular binary channels.
c
      Zcomp=0
      Ncomp=0
      do 30 type=0,6
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        if (nlev(Zix,Nix).eq.0) nlev(Zix,Nix)=nlevbin(type)
   30 continue
      do 40 Zix=0,numZ
        do 40 Nix=0,numN
          if (Zix.eq.parZ(k0).and.Nix.eq.parN(k0)) goto 40
          if (nlev(Zix,Nix).eq.0) nlev(Zix,Nix)=nlevmaxres
   40 continue
      nlev(parZ(k0),parN(k0))=nlevmax
c
c ************ Nuclear structure for first set of nuclides *************
c
c strucinitial: subroutine for initialization of arrays for various
c               structure parameters
c masses      : subroutine for nuclear masses
c separation  : subroutine for separation energies and reduced and
c               specific masses
c primary     : flag to designate primary (binary) reaction
c flagpartable: flag for output of model parameters on separate file
c Einc        : incident energy
c enincmin    : minimum incident energy
c parskip     : logical to skip outgoing particle
c structure   : subroutine for nuclear structure parameters
c odd         : odd (1) or even (0) nucleus
c weakcoupling: subroutine for weak coupling model
c flagcomp    : flag for compound nucleus calculation
c parinclude  : logical to include outgoing particle
c radwidtheory: subroutine for theoretical calculation of total
c               radiative width
c flaggiant0  : flag for collective contribution from giant resonances
c sumrules    : subroutine for giant resonance sum rules
c targetspin  : spin of target
c jdis        : spin of level
c Ltarget     : excited level of target
c targetspin2 : 2 * spin of target
c targetP     : parity of target
c parlev      : parity of level
c targetE     : energy of target
c edis        : energy of level
c
c The nuclear masses and separation energies are read/calculated first,
c for all nuclides that can possibly be formed in multiple reactions.
c For all nuclides that can be formed by the first binary reaction, we
c determine the nuclear structure information (discrete levels,
c parameters for resonances, photons, fission and level densities). The
c logical strucexist ensures that we only need to do this once for each
c nucleus. (i.e., we do not waste time if the same nucleus is
c encountered later in the reaction chain).
c All parameters for each nucleus are written to file 'parameters.dat'.
c
      call strucinitial
      call masses
      call separation
      primary=.true.
      if (flagpartable) open (11,status='unknown',file='parameters.dat')
      Einc=enincmin
      do 110 type=0,6
        if (parskip(type).and.type.ne.0) goto 110
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        call structure(Zix,Nix)
        strucexist(Zix,Nix)=.true.
        A=AA(Zcomp,Ncomp,type)
        odd=mod(A,2)
        if (type.eq.k0.and.odd.eq.1) call weakcoupling(Zix,Nix,type)
  110 continue
      if (parinclude(0).or.flagcomp) call radwidtheory(Zcomp,Ncomp,0.)
      if (flaggiant0) call sumrules
      targetspin=jdis(parZ(k0),parN(k0),Ltarget)
      targetspin2=int(2.*targetspin)
      targetP=parlev(parZ(k0),parN(k0),Ltarget)
      targetE=edis(parZ(k0),parN(k0),Ltarget)
c
c ******************* Q-values and Coulomb barriers ********************
c
c tarmass : mass of target nucleus
c nucmass : mass of nucleus
c Q       : Q-value for target nucleus
c S       : separation energy per particle
c coulbar : Coulomb barrier
c Ztarget : charge number of target nucleus
c e2      : square of elementary charge in MeV.fm
c Atarget : mass number of target nucleus
c onethird: 1/3
c
      tarmass=nucmass(parZ(k0),parN(k0))
      do 210 type=0,6
        if (parskip(type)) goto 210
        Q(type)=S(0,0,k0)-S(0,0,type)
        coulbar(type)=Ztarget*parZ(type)*e2/(1.25*Atarget**onethird)
  210 continue
c
c * Off- and on-set energies for preequilibrium and width fluctuations *
c
c ewfc   : off-set incident energy for width fluctuation calculation
c eurr   : off-set incident energy for URR calculation
c flagurr: flag for output of unresolved resonance parameters
c
c Assignment of default values: Width fluctuation corrections and URR
c are included for incident energies up to the separation energy.
c
      if (k0.ge.1.and.ewfc.eq.-1.) ewfc=S(parZ(k0),parN(k0),k0)
      if (k0.eq.1.and.eurr.eq.-1..and.flagurr)
     +  eurr=S(parZ(k0),parN(k0),k0)
c
c Pre-equilibrium reactions are included for incident energies
c above the last discrete level.
c
c NLast     : last discrete level
c epreeq    : on-set incident energy for preequilibrium calculation
c kalbachsep: subroutine for separation energy for Kalbach systematics
c
      NL=Nlast(parZ(k0),parN(k0),0)
      if (epreeq.eq.-1.) epreeq=max(edis(parZ(k0),parN(k0),NL),1.)
      call kalbachsep
c
c For astrophysical calculations, the incident energy grid is determined
c
c flagastro : flag for calculation of astrophysics reaction rate
c egridastro: subroutine to calculate default incident energy grid
c             for astrophysical rate
c flagprod  : flag for isotope production
c decaydata : subroutine for decay data
c
      if (flagastro) call egridastro
      if (flagprod) call decaydata
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
