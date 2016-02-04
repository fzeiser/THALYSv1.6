      subroutine csunr2(ay,l,j)
c
c +---------------------------------------------------------------------
c | Author: Gilles Noguere
c | Date  : December 9, 2013
c | Task  : Unresolved resonance cross section from NJOY
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,j,mu,nu,lamda,odd
      real    gnox,gxx,ggx,gfx,dx,gnx,gam,galpha,gbeta,diff,den,amun,
     +        e,sqrte,k,const,temp,aj,ll,spi,gj,ac,ay,rho,
     +        rhoc,vl,ps,terg,ters,terf,gs,gc,gff,add
c
c      link with NJOY parameters
c
      spot(l)=0.
      gnox=urrwidth(3,l,j)
      gxx=urrwidth(1,l,j)
      ggx=urrwidth(0,l,j)
      gfx=urrwidth(-1,l,j)
      dx=Dlj(l,j)
      if (dx.eq.0.) return
      amun=real(nulj(0,l,j))
      mu=min(max(nulj(0,l,j),1),4)
      nu=min(max(nulj(-1,l,j),1),4)
      lamda=min(max(nulj(1,l,j),1),4)
      e=Einc*1.e6
      sqrte=sqrt(e)
      k=wavenum*10.
      const=2.*e*pi**2/k**2
      odd=mod(Atarget+1,2)
      aj=j+0.5*odd
      ll=real(l)
      spi=jdis(0,1,0)
      gj=(2.*aj+1.)/(4.*spi+2.)
      ac=0.123*tarmass**onethird+0.08
      rho=k*ac
      rhoc=k*ay
c
c      calculate penetrability (vl) and phase shift(ps)
c
      call unfac(l,rho,rhoc,amun,vl,ps)
      vl=vl*sqrte
c
c       calculate potential scattering
c
      spot(l)=2.*twopi*(2.*ll+1.)*(sin(ps)/k)**2
c
c      compute cross section contributions
c
      gnx=gnox*vl
      gam=ggx
      galpha=gnx
      gbeta=gfx
      diff=gxx
      den=e*dx
      temp=const*gj*gnx/den
      terg=temp*gam
      ters=temp*gnx
      terf=temp*gbeta
c
c      calculate fluctuation integrals
c
      call gnrl(galpha,gbeta,gam,mu,nu,lamda,gs ,diff,1)
      call gnrl(galpha,gbeta,gam,mu,nu,lamda,gc ,diff,2)
      call gnrl(galpha,gbeta,gam,mu,nu,lamda,gff,diff,3)
      gc=gc*terg
      gff=gff*terf
      gs=gs*ters
c
c      add interference correction term
c
      add=const*gj*2.*gnx*(sin(ps))**2
      add=add/(e*dx)
      gs=gs-add
c
c      cross sections
c
      sigurrs(l,j)=gs
      sigurrf(l,j)=gff
      sigurrc(l,j)=gc
      return
      end
