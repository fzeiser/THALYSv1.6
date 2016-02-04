      subroutine inverse(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 19, 2013
c | Task  : Calculation of total, reaction and elastic cross sections
c |         and transmission coefficients for outgoing particles and
c |         energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical lexist
      integer Zcomp,Ncomp,Z,A
c
c ************************** ECIS calculation **************************
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c transfile: file with transmission coefficients
c csfile   : file with inverse reaction cross sections
c ZZ,Z     : charge number of residual nucleus
c AA,A     : mass number of residual nucleus
c
c All transmission coefficients calculated by ECIS will be written on a
c file trZZZAAA, where ZZZ and AAA are the charge and mass number in
c (i3.3) format. The reaction cross sections will be written to
c csZZZAAA.
c
      transfile='tr000000     '
      csfile='cs000000     '
      Z=ZZ(Zcomp,Ncomp,0)
      A=AA(Zcomp,Ncomp,0)
      if (Z.lt.10) then
        write(transfile(5:5),'(i1)') Z
        write(csfile(5:5),'(i1)') Z
      endif
      if (Z.ge.10.and.Z.lt.100) then
        write(transfile(4:5),'(i2)') Z
        write(csfile(4:5),'(i2)') Z
      endif
      if (Z.ge.100) then
        write(transfile(3:5),'(i3)') Z
        write(csfile(3:5),'(i3)') Z
      endif
      if (A.lt.10) then
        write(transfile(8:8),'(i1)') A
        write(csfile(8:8),'(i1)') A
      endif
      if (A.ge.10.and.A.lt.100) then
        write(transfile(7:8),'(i2)') A
        write(csfile(7:8),'(i2)') A
      endif
      if (A.ge.100) then
        write(transfile(6:8),'(i3)') A
        write(csfile(6:8),'(i3)') A
      endif
c
c Calculate transmission coefficients and inverse reaction cross
c sections.
c
c flagecisinp: flag for existence of ecis input file
c invexist   : logical to determine necessity of new inverse cross
c              section and transmission coefficients calculation
c inverseecis: subroutine for ECIS calculation for outgoing particles
c              and energy grid
c
      if (.not.invexist(Zcomp,Ncomp)) call inverseecis(Zcomp,Ncomp)
      invexist(Zcomp,Ncomp)=.true.
c
c Modification 5/5/11 by Kevin Kelley
c If the user has specified 'eciscalc n' and a cs/tr file pair is not
c present, go ahead and call the inverseecis subroutine to produce them.
c This is a little more user friendly than terminating 3/4 of the way
c through a calculation.
c This may save a lot of time in cases with 'optmodall y' when many
c cross sections and transmission coefficients are calculated, and
c which could then be transferred from directory to directory.
c
      flagecisinp=.true.
      if (.not.flageciscalc) then
        inquire(file=csfile,exist=lexist)
        if (lexist) inquire(file=transfile,exist=lexist)
        if (.not.lexist) then
          flageciscalc=.true.
          call inverseecis(Zcomp,Ncomp)
          flageciscalc=.false.
        endif
        invexist(Zcomp,Ncomp)=.true.
      endif
c End modification 5/5/11 by Kevin Kelley
c
c Read transmission coefficients and inverse reaction cross sections.
c
c inverseread: subroutine to read ECIS results for outgoing particles
c              and energy grid
c inversenorm: subroutine for normalization of reaction cross sections
c              and transmission coefficients
c flaginverse: flag for output of transmission coefficients and inverse
c              reaction cross sections
c inverseout : subroutine for reaction output for outgoing channels
c
      if (flagecisinp) call inverseread(Zcomp,Ncomp)
      call inversenorm(Zcomp,Ncomp)
      if (flaginverse) call inverseout(Zcomp,Ncomp)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
