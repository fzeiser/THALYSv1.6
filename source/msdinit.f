      subroutine msdinit
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : Initialization of MSD model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nen
c
c ********* Set parameters and energy grid for MSD calculation *********
c
c maxmsd     : number of MSD steps
c maxpar     : maximal particle number
c maxJmsd    : maximal spin for MSD calculation
c Emsdmin    : minimal outgoing energy for MSD calculation
c Einc       : incident energy in MeV
c msdbins2   : number of energy points for MSD calculation
c msdbins    : number of energy points for DWBA calculation for MSD
c dEmsd      : energy bin for MSD
c Emsd       : MSD energy grid
c specmass   : specific mass for target nucleus
c parZ       : charge number of particle
c k0         : index of incident particle
c parN       : neutron number of particle
c flagonestep: flag for continuum one-step direct only
c flagddx    : flag for output of double-differential cross sections
c interangle : subroutine for intermediate angles by addition theorem
c              for MSD model
c
      maxmsd=maxpar-1
      maxJmsd=6
      if (Emsdmin.eq.0.or.Emsdmin.ge.Einc) Emsdmin=Einc/5.
      msdbins2=msdbins*2
      dEmsd=(Einc-Emsdmin)/msdbins2
      do 10 nen=0,msdbins2
        Emsd(nen)=Einc-nen*dEmsd
   10 continue
      Emsd(0)=real(Emsd(0)/specmass(parZ(k0),parN(k0),k0))
      if (.not.flagonestep.and.flagddx) call interangle
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
