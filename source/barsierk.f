      subroutine barsierk(iz,ia,il,bfis,egs,lbar0)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : July 2, 2004
c | Task  : Fission barrier heights, rotating gs energy and lbar0
c +---------------------------------------------------------------------
c
c *************************** Comments ********************************
c
c This subroutine returns the fission barrier height in MeV. It is
c based on calculations using yukawa-plus-exponential double folded
c nuclear energy, exact couloumb diffuseness corrections, and diffuse-
c matter moments of inertia (A.J. Sierk, Phys. Rev. C33, 2039 (1986)).
c The implementation is analogous to the subroutine "asierk" written
c by A.J. Sierk (LANL, 1984).
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer iz,ia,il,iloop,jloop,kloop
      real    bfis,egs,lbar0,zchar,amass,ll,amin,amax,a,z,l,plegendre,
     +        amin2,amax2,l80,l20,x,y,p1,p2,q1,q2,q3,r,a1,a2
c
c ******************** Rotating Finite Range Model *********************
c
c bfis     : barrier height
c egs      : rotating ground state energy
c lbar0    : l-value for which bfis becomes zero
c iz,zchar : charge number of residual nucleus
c ia,amass : mass number of residual nucleus
c il,ll    : angular momentum
c amax,amin: maximal and minimal mass number defining range of the fit
c a,z,l    : fit variables
c plegendre: function for calculation of Legendre polynomial
c
c Barrier height for l=0
c
      bfis=0.
      egs=0.
      lbar0=0.
      if (iz.lt.19.or.iz.gt.111) goto 900
      if (iz.gt.102.and.il.gt.0) goto 910
      zchar=real(iz)
      amass=real(ia)
      ll=real(il)
      amin=1.2*zchar+0.01*zchar*zchar
      amax=5.8*zchar-0.024*zchar*zchar
      if (amass.lt.amin.or.amass.gt.amax) goto 920
      a=amass/400.
      z=zchar/100.
      l=ll/100.
      do 10 iloop=1,7
        do 10 jloop=1,7
          bfis=bfis+barcof(jloop,iloop)*
     +      plegendre(jloop-1,z)*plegendre(iloop-1,a)
 10   continue
      if (il.lt.1) return
c
c L-values corresponding to fission barrier height which is 20% (80%)
c of L=0 fission barrier
c
c l20,l80: l-value for which bfis is 20% (80% resp.) of bfis(l=0)
c
      l80=0.
      l20=0.
      amin2=1.4*zchar+0.009*zchar*zchar
      amax2= 20.+3.0*zchar
      if ((amass.lt.amin2-5..or.amass.gt.amax2+10.).and.il.gt.0)
     +  goto 920
      do 20 iloop=1,4
        do 20 jloop=1,5
          l80=l80+l80cof(jloop,iloop)*
     +      plegendre(jloop-1,z)*plegendre(iloop-1,a)
          l20=l20+l20cof(jloop,iloop)*
     +      plegendre(jloop-1,z)*plegendre(iloop-1,a)
 20   continue
c
c L-value for which the fission barrier vanishes
c
      do 30 iloop=1,4
        do 30 jloop=1,6
          lbar0=lbar0+lmxcof(jloop,iloop)*
     +      plegendre(jloop-1,z)*plegendre(iloop-1,a)
 30   continue
c
c L-dependent fission barrier
c
c x,y,p1,p2,r  : help variables
c q1,q2,q3     : help variables
c barcof,l80cof: fitted parameters
c l20cof,lmxcof: fitted parameters
c egscof       : fitted parameters
c
      x=l20/lbar0
      y=l80/lbar0
      if (ll.gt.l20) then
         p1=(-(20.*x**5)+25.*x**4-4.)*(y-1.)**2*y*y
         p2=(-(20.*y**5)+25.*y**4-1.)*(x-1.)**2*x*x
         q1=0.2/((y-x)*((1.-x)*(1.-y)*x*y)**2)
         q2=q1*(p1*y-p2*x)
         q3=-(q1*(p1*(2.*y+1.)-p2*(2.*x+1.)))
         r=ll/lbar0
         a1=4.*r**5-5.*r**4+1.
         a2=q2*(2.*r+1.)
         bfis=bfis*(a1+(r-1.)*(a2+q3*r)*r*r*(r-1.))
      else
         q1=0.2/(l20**2*l80**2*(l20-l80))
         q2=q1*(4.*l80**3-l20**3)
         q3=-(q1*(4.*l80**2 -l20**2))
         bfis=bfis*(1.+q2*ll**2+q3*ll**3)
      endif
      if (bfis.le.0..or.ll.gt.lbar0) bfis=0.
c
c Rotating ground state energy
c
      if(ll.gt.lbar0)return
      do 40 iloop=1,4
        do 40 jloop=1,6
          do 40 kloop=1,5
            egs=egs + egscof(kloop,jloop,iloop)*
     +        plegendre(jloop-1,z)*plegendre(iloop-1,a)*
     +        plegendre(2*kloop-2,l)
 40   continue
      if (egs.lt.0.) egs = 0.
c
c warning messages for attempted use outside validity boundaries
c
      return
 900  write(*,100)
      return
 910  write(*,110)
      return
 920  write(*,120)
      return
 100  format(/,10x,'*  *  *  barfit called with  z  less than 19 or ',
     1 ' greater than 111.  bfis is set to 0.0  *  *  *')
 110  format(/,10x,'*  *  *  barfit called with  z  greater than 102',
     1 ' and  l  not equal to zero.  bfis is set to 0.0,  *  *  *')
 120  format(/,10x,'*  *  *  barfit called with a =',i3,', outside ',
     1 'the allowed values for z = ',i3,' *  *  *')
      end
