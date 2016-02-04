      subroutine foldalpha(Zix,Nix,E)
c
c Global alpha potential = real double folding potential + imaginary WS potential
c   dispersive relations also included if alphaomp=5
c Ref: Demetriou, Grama and Goriely (2002) NPA707, 253
c
      include "talys.cmb"
      character*80 key
      real rhomom(numjlm),radmom(numjlm)
      real va(1000),rva(1000)
      real z,a,ee,rb,e2exp,factor1,factor2
      real vopr(numjlm),wopr(numjlm)
      real alphav(10)
      real e,a2,a3,a13,as,es,aaw,raw,te,expj0,ww,dwv,dws,
     &  efermia,rww,aww,www,ws,rws,aws,wws,
     &  vdiv,vdis,am1,am2,am3,rr,c,hh,aa1,bb1,v5d,eref,dvolj
      integer Zix,Nix,i,k,nu,nradrho,klo,khi,kk
c
c ************************ Alpha Optical potential *************************
c
c Zix        : charge  number index of target nucleus
c Nix        : neutron number index of target nucleus
c E          : incident energy
c Z          : charge number of target nucleus
c A          : mass   number of target nucleus
c e2exp      : lowest inelastic threshold
c radmom     : radial grid for the potential, identical as the one used for JLM
c rhomom     : total density (neutron + proton) at a given radius radmom
c alphaomp   : model for the alpha OMP
c            : 3=OMP I, 4=OMP II, 5=OMP III of Demetriou et al. (2002)
c vopr       : Final real part of the alpha optical potential
c wopr       : Final imaginary part of the alpha optical potential
c
      a=AA(Zix,Nix,0)
      z=ZZ(Zix,Nix,0)
      ee=E
      do i=1,numjlm
        radmom(i)=radjlm(Zix,Nix,i)
        rhomom(i)=rhojlmn(Zix,Nix,i,1)+rhojlmp(Zix,Nix,i,1)
      enddo
      e2exp=edis(Zix,Nix,1)
      a2=a**2
      a3=a**3
      a13=a**(1./3.)
c
c
c--------------------------------------------------------------------------------------
c alphav =  1: reduced radius r0=r/a13 of real pot
c           2: diffuseness a of woods-saxon
c           3: reduced radius r0s=r/a13 of imag. pot.
c           4: diffuseness as of imag. pot.
c           5: fraction of surface-peaked im. pot. relative to volume part
c           6: depth vso of spin-orbit pot.
c           7: depth v of real pot.
c           8: depth w of imag. pot.
c rww        : radius of imaginary volume potential
c aww        : diffuseness of imaginary volume potential
c www        : depth of imaginary volume potential
c rws        : radius of imaginary surface potential
c aws        : diffuseness of imaginary surface potential
c wws        : depth of imaginary surface potential
c te         : energy-dependent term of the imaginary volume integral
c as         : diffuseness of the functional form of the imaginary volume integral
c es         : energy threshold of the functional form of the imaginary volume integral
c expj0      : saturation value of the imaginary volume integral
c--------------------------------------------------------------------------------------
c
      v5d=0.
      goto (10,20,30) alphaomp-2
   10 continue
c--------------------------------------------------------------------------------------
c Case alphaomp=3
c Imaginary potential OMP I of Demetriou et al. (2002)
c Woods-Saxon E-dependent volume imaginary potential
c--------------------------------------------------------------------------------------
      alphav(1)=1.25
      raw=0.85281+0.02202*a-2.14551e-4*a2+7.92942e-7*a3-
     1     9.94658e-10*a2*a2
      aaw=-0.13526+0.02029*a-1.98441e-4*a2+7.35104e-7*a3-
     2     9.15272e-10*a2*a2
      if (aaw.le.1.e-2) aaw=1.e-2
      as=7.65867-7.5669*e2exp+2.50486*e2exp**2
      es=0.1024*a+1.1307*as
      te=1./(1.+exp(-(ee-es)/as))
      expj0=77.
      if (a.le.90.) expj0=135.-0.644*a
      ww=3./pi/raw**3*expj0*te
      alphav(3)=raw
      alphav(4)=aaw
      alphav(8)=ww
      alphav(5)=0.
      dwv=0.
      dws=0.
      goto 50
c
   20 continue
c--------------------------------------------------------------------------------------
c Case alphaomp=4
c Imaginary potential OMP II of Demetriou et al. (2002)
c Woods-Saxon E-dependent volume+surface imaginary potential
c--------------------------------------------------------------------------------------
      alphav(1)=1.25
      raw=1.47385-0.00134615*a
      aaw=0.29
      alphav(5)=0.9
      as=7.65867-7.5669*e2exp+2.50486*e2exp**2
      es=0.1024*a+1.1307*as
      te=1./(1.+exp(-(ee-es)/as))
      expj0=77.
      if(a.le.90.) expj0=135.-0.644*a
      ww=expj0*te/(raw**3/3.+7.603*alphav(5)*aaw*raw**2/a13)/pi
      alphav(3)=raw
      alphav(4)=aaw
      alphav(8)=ww
      dwv=0.
      dws=0.
      goto 50
c
   30 continue
c--------------------------------------------------------------------------------------------------
c Case alphaomp=5
c Imaginary potential OMP III of Demetriou et al. (2002)
c Woods-Saxon E-dependent volume+surface imaginary potential plus Dispersion contribution
c  using prescription and functions of Capote et al. J. Phys. G 27 (2001) B15
c
c efermia    : Fermi energy
c am1        : mass excess of the target-projectile system
c am2        : mass excess of the target+projectile system
c am3        : mass excess of the projectile
c dwv        : dispersive contribution to the real potential arising from the volume imaginary pot
c dws        : dispersive contribution to the real potential arising from the surface imaginary pot
c--------------------------------------------------------------------------------------------------
      if (expmexc(Zix-2,Nix-2).ne.0..and.expmexc(Zix+2,Nix+2).ne.0.)
     &  then
        am1=expmexc(Zix+2,Nix+2)
        am2=expmexc(Zix-2,Nix-2)
      else
        am1=thmexc(Zix+2,Nix+2)
        am2=thmexc(Zix-2,Nix-2)
      endif
      am3=(parmass(6)-4.)*amu
      efermia=-0.5*(am1-am2+2.*am3)
c
      alphav(1)=1.25
      raw=1.47385-0.00134615*a
      aaw=0.30
      alphav(5)=0.35
      v5d=alphav(5)
      c=0.005
      if(ee.lt.13.) c=-0.165*ee+2.15
      alphav(5)=alphav(5)*exp(-c*abs(ee-efermia))
c
      as=7.65867-7.5669*e2exp+2.50486*e2exp**2
      es=0.0854*a+1.1307*as
      te=1./(1.+exp(-(ee-es)/as))
      expj0=75.
      if(a.le.96.) expj0=135.-0.644*a
c volume integral
      ww=expj0*te/((1.-alphav(5))*raw**3/3.+7.603*alphav(5)*aaw*
     &  raw**2/a13)/pi
      alphav(3)=raw
      alphav(4)=aaw
      alphav(8)=ww
c
c dispersive contributions: dwv = volume; dws = surface; dvolj = volume integral
c new method for dispersive relation following Mahaux, Ngo and Satchler, NPA 449 (1986) 354
      eref=150.
      call mahaux(a,ee,eref,expj0,as,es,v5d,raw,aaw,dwv,dws,dvolj,
     &  efermia)
c
   50 continue

c final radius, diffuseness and depth of the imaginary WS-type potential
c They are multiplied by the OMP adjustment keywords of TALYS
c (default 1.)
c   volume term
      rww=rvadjust(6)*alphav(3)*a13
      aww=avadjust(6)*alphav(4)
      www=w1adjust(6)*alphav(8)*(1.-v5d)
c   surface term
      rws=rwdadjust(6)*1.09*rww
      aws=awdadjust(6)*1.6*aww
      wws=d1adjust(6)*alphav(8)*alphav(5)

c--------------------------------------------------------------------------------------------------
c determination of the real folding potential through the product of the fourier transforms
c rb         : maximum radius value
c nu         : number of radial grid point
c rva        : radial coordinate for the double folding potential
c va         : double folding potential
c nradrho    : number of grid point used in radmom and rhomom
c vdiv       : dispersive contribution to the real potential arising from the volume imaginary pot
c vdis       : dispersive contribution to the real potential arising from the surface imaginary pot
c ws         : surface component of the imaginary potential
c rc         : coulomb radius based on elton's formula
c--------------------------------------------------------------------------------------------------
      rb=radjlm(Zix,Nix,numjlm)
      nu=1000
      nradrho=numjlm

      call afold(z,a,ee,rb,nu,rva,va,radmom,rhomom,nradrho)
c
c Depth and shape are multiplied by the OMP adjustment keywords of TALYS
c (default 1.)
c
c ompadjustp: flag for local optical model parameter adjustment
c adjust    : subroutine for energy-dependent parameter adjustment
c factor1,2 : multiplication factor
c
      if (ompadjustp(6)) then
        key='aradialcor'
        call adjust(E,key,Zix,Nix,0,0,factor1)
        key='adepthcor'
        call adjust(E,key,Zix,Nix,0,0,factor2)
      else
        factor1=1.
        factor2=1.
      endif
      do i=1,nu
        rva(i)=factor1*aradialcor*rva(i)
        va(i)=factor2*adepthcor*va(i)
      enddo
c
c interpolation of the alpha optical potential va(rva) on the given radial grid radjlm
c
      do 100 k=1,numjlm
        rr=radjlm(Zix,Nix,k)
        klo=1
        khi=nu
        if (rr.le.rva(klo)) then
           khi=klo+1
           goto 110
        endif
        if (rr.ge.rva(khi)) then
           klo=khi-1
           goto 110
        endif
  120   if (khi-klo.gt.1) then
          kk=(khi+klo)/2.
          if (rva(kk).gt.rr) then
            khi=kk
          else
          klo=kk
          endif
          goto 120
        endif
  110   hh=rva(khi)-rva(klo)
        aa1=(rva(khi)-rr)/hh
        bb1=(rr-rva(klo))/hh
        vopr(k)=aa1*va(klo)+bb1*va(khi)
c dispersive contributions to real potential for alphaomp=5  (OMP III)
        if (alphaomp.lt.5) then
          vdiv=0.
          vdis=0.
        else
          if (abs((rr-rww)/aww).lt.88.) then
            vdiv=-dwv/(1.+exp((rr-rww)/aww))
          else
            vdiv=0.
          endif
          if (abs((rr-rws)/aws).lt.88.) then
            vdis=-4.*dws*exp((rr-rws)/aws)/(1.+exp((rr-rws)/aws))**2
          else
            vdis=0.
          endif
        endif
        vopr(k)=vopr(k)+vdiv+vdis
c final imaginary potential
        wopr(k)=0.
        if (abs((rr-rww)/aww).lt.88.)
     &    wopr(k)=-www/(1.+exp((rr-rww)/aww))
        if (abs((rr-rws)/aws).lt.88.) then
          ws=-4.*wws*exp((rr-rws)/aws)/(1.+exp((rr-rws)/aws))**2
        else
          ws=0.
        endif
        wopr(k)=wopr(k)+ws
c-----------------------------------------------------------------
c potjlm   1 : final real central potential
c potjlm   2 : final imaginary central potential
c-----------------------------------------------------------------
        potjlm(Zix,Nix,k,1)=vopr(k)
        potjlm(Zix,Nix,k,2)=wopr(k)
        potjlm(Zix,Nix,k,3)=0.
        potjlm(Zix,Nix,k,4)=0.
        potjlm(Zix,Nix,k,5)=0.
        potjlm(Zix,Nix,k,6)=0.
  100 continue
c calculate the coulomb radius based on elton's formula
      rc=1.123+2.352*(a**(-.666666))-2.07*(a**(-1.333333))
      return
      end
c
c---------------------------------------------------------------------------
c subroutine for double folding potential for alpha projectiles
c---------------------------------------------------------------------------
c
      subroutine afold(z,a,e,rb,mp,rv,u,radmom,rhomom,nrad)
c
c determination of the real folding potential through the product of the fourier transforms
c Z          : charge number of target nucleus
c A          : mass   number of target nucleus
c E          : incident energy
c rb         : maximum radius value
c mp         : number of radial grid point
c rv         : radial coordinate for the double folding potential
c u          : double folding potential
c nrad       : number of grid point used in radmom and rhomom
c radmom     : radial grid for the potential, identical as the one used for JLM
c rhomom     : total density (neutron + proton) at a given radius radmom
c
c pass       : logical boolean to estimate alpha density distribution only once
c ce         : amplitude of the density-independent term of the DDM3Y interaction
c alpha      : amplitude of the density-dependent term of the DDM3Y interaction
c beta       : density-dependenceterm of the DDM3Y interaction
c am(1),amp  : mass number of projectile
c am(2),amt  : mass number of target
c am(3)      : charge number of projectile
c am(4)      : charge number of target
c r1,r2,ru   : maximum radius value
c n1,n2,nu   : number of coordinate grid points
c d1,d2,du   : equidistant radial increment
c rq         : maximum frequency value in the Fourier space
c nq         : number of frequency points in the Fourier space
c dq         : frequency increment in the Fourier space
c de         : single-nucleon exchange term J00
c ftrans     : subroutine performing the Fourier transform
c
      parameter(ngrid=1000)
      real rhomom(nrad),radmom(nrad)
      real*8 bcof,z0,d1,d2,dq,q,r,bint,b
      common /sp/ dq,d1,d2,nr(2),nq,nu,du,fro(ngrid,2),ro(ngrid,2),
     &  am(4),roe(ngrid,2),beta
      common/constf/ pi,sqrpi,pi4,pi32
      dimension bcof(ngrid),bint(9),b(1),ue(ngrid)
      dimension u(mp),rv(mp)
      dimension a1(6),a2(6),b1(6),b2(6)
      logical pass
      save
      equivalence (am(1),amp),(am(2),amt),(nr(1),n1),(nr(2),n2)
      data a1/1577.1071d0,34751.04d0,19839.156d0,25129.599d0,25129.599d0
     *,-15348.25d0/,b1/0.16,16.0,16.,16.,16.,16./,a2/-1239.9272d0,
     * -12754.865d0,-9857.0604d0,-10726.653d0,-10726.653d0,5908.71d0/,
     * b2/0.0625,6.2500,6.25,6.25,6.25,6.25/
      data pass/.false./
      pi=3.1415926536d0
      sqrpi=sqrt(pi)
      pi4=4.*3.1415926536d0
      pi32=sqrpi*pi
      nu=mp
      key=5
      ce=0.44073
      alpha=4.3529
      beta=10.639
      am(1)=4
      am(2)=a
      am(3)=2
      am(4)=z
      r1=rb
      d1=r1/nu
      r2=rb
      d2=r2/nu
c id=1 if d1=d2
      id=1
      rq=10.
      dq=rq/nu
      ru=rb
      nr(1)=nu
      nr(2)=nu
      nq=nu
      du=ru/nu
      do 34 i=1,nu
        u(i)=0.d0
        ue(i)=0.d0
   34 continue
      akoff=0.5/(pi*pi)
      is1=1
      is2=1
      do 33 iso=is1,is2
        if(iso.eq.2) key=6
        de=0.
        goto (20,20,20,20,15,16),key
   15   de=-276.*(1.-0.005*e/amp)
        goto 20
   16   de=227*(1.-0.004*e/amp)
   20   continue
        if (pass) go to 25
        no=1
        if(amp.gt.1.) call dense(no)
        pass=.true.
   25   continue
        no=2
        r=0.
        do 42 i=1,nu
c analytic density distribution as a fermi function
c         ro(i,2)=(cn/(1.+exp((r-rn)/bn))+cz/(1.+exp((r-rz)/bz)))*r*r*pi4
c grid density distribution
c  --> determination of the density at radius r
          klo=1
          khi=nrad
          if (r.le.radmom(klo)) then
            aa1=1.
            bb1=0.
            goto 465
          endif
          if (r.ge.radmom(khi)) then
            aa1=0.
            bb1=1.
            goto 465
          endif
  470     if (khi-klo.gt.1) then
            kk=(khi+klo)/2.
            if (radmom(kk).gt.r) then
              khi=kk
            else
              klo=kk
            endif
            goto 470
          endif
          hh=radmom(khi)-radmom(klo)
          aa1=(radmom(khi)-r)/hh
          bb1=(r-radmom(klo))/hh
  465     rho=aa1*rhomom(klo)+bb1*rhomom(khi)
c
          ro(i,2)=rho*r*r*pi4
          roe(i,2)=ro(i,no)*exp(-ro(i,no)*beta)
          r=r+d2
   42   continue
        nmint=max0(nr(1),nr(2))
        if(amp.eq.1.) nmint=nr(2)
        iflg=0
        q=0.
        do 200 k=1,nq
          qq=q*q
          sum1=0.
          sum2=0.
          s1=0.
          s2=0.
          i1=1
          z0=0.
          ndd=nu
          if(id.eq.1.and.iflg.eq.0)
     &      call ftrans(bcof,i1,nmint,z0,d1,q,bint,b,ndd)
           if(id.ne.1.and.iflg.eq.0)
     &      call ftrans(bcof,i1,n1,z0,d1,q,bint,b,ndd)
          if(key.ne.5) goto 101
          do 100 j=1,n1
  100     s1=s1+bcof(j)*roe(j,1)
  101     do 102 j=1,n1
  102     sum1=sum1+bcof(j)*ro(j,1)
  126     if(id.ne.1.and.iflg.eq.0)
     &      call ftrans(bcof,i1,n2,z0,d2,q,bint,b,ndd)
          if(key.ne.5) goto 131
          do 130 j=1,n2
  130     s2=s2+bcof(j)*roe(j,2)
  131     do 132 j=1,n2
  132     sum2=sum2+bcof(j)*ro(j,2)
  172     c=a1(key)/(qq+b1(key))+a2(key)/(qq+b2(key))+de
          fro(k,1)=(sum1*sum2+s1*s2*alpha)*c*qq
  200   q=q+dq
        r=0.
        do 300 i=1,nu
          sum=0.
          call ftrans(bcof,i1,nq,z0,dq,r,bint,b,ndd)
          do 400 k=1,nq
  400     sum=sum+bcof(k)*fro(k,1)
          u(i)=akoff*sum*ce+u(i)
          rv(i)=r
  300   r=r+du
   33 continue
      r=0
      do 520 i=1,nu
        ud=u(i)
        u(i)=ud+ue(i)
  520 r=r+du
      call intu(u,nu,du,volj,rms,vmi33)
      volj=volj/(amp*amt)
      rms=sqrt(rms)
c Renormalization of the potential depth through the volume integral
c  cf demetriou et al. (2002) NPA707, page 258
      cee=(e-19.)**2/19000.
      if (abs(cee).lt.88.) then
        cej=337.2/exp(cee)
        rap=-cej/volj
      else
        rap=1.
      endif
      do 521 i=1,nu
        u(i)=rap*u(i)
  521 continue
      return
      end
c
      subroutine intu(u,nr,du,volj,rms,rin)
c----------------------------------------------------------------------------------
c determination of the Volume Integral and root-mean-square radius of the potential
c----------------------------------------------------------------------------------
      dimension u(nr)
c----------------------------------------------------------------
c u          : potential
c nr         : number of radial grid point
c du         : equidistant radial increment
c volj       : volume integral of the potential
c rms        : root-mean-square radius of the potential
c rin        : mean radius of the potential
c----------------------------------------------------------------
c
      fact=0.375*du*12.5663706144d0
      volj=0.
      rms=0.
      rin=0.
      sum0=0.
      sum1=0.
      sum2=0.
      a1=0.
      a2=0.
      r=0.
      nr2=nr-1
      nr1=(nr2/3)*3
      do 400 kr=3,nr1,3
        r1=r+du
        r2=r1+du
        r3=r2+du
        rr1=r1*r1
        rr2=r2*r2
        rr3=r3*r3
        b1=u(kr-1)*rr1+u(kr)*rr2
        c1=u(kr+1)*rr3
        sum1=sum1+a1+3.*b1+c1
        a1=c1
        rrr1=rr1*rr1
        rrr2=rr2*rr2
        rrr3=rr3*rr3
        b2=u(kr-1)*rrr1+u(kr)*rrr2
        c2=u(kr+1)*rrr3
        sum2=sum2+a2+3.*b2+c2
        a2=c2
        r=r3
  400 continue
      if (nr2.le.nr1) goto 3
      nr1=nr1+1
      c5=du*4.d0/3.d0
      do 4 kr=nr1,nr2
        sum0=sum0+c5*(u(kr)*(kr-1)+u(kr+1)*kr)
        sum1=sum1+c5*du*(u(kr)*(kr-1)**2+u(kr+1)*kr**2)
        sum2=sum2+c5*du**3*(u(kr)*(kr-1)**4+u(kr+1)*kr**4)
    4 continue
    3 rin=fact*sum0
      volj=fact*sum1
      if (sum1.ne.0.) rms=sum2/sum1
      return
      end
c
      subroutine dense(n0)
c
c determination of radial density distribution of the alpha system
c  as a sum of Gaussians
      parameter(ngrid=1000)
      real*8 dq,d1,d2
      common /sp/ dq,d1,d2,nr(2),nq,nu,du,fro(ngrid,2),ro(ngrid,2),
     &  am(4),roe(ngrid,2),beta
      common /constf/ pi,sqrpi,pi4,pi32
      dimension a(12),ri(12),sq(12)
      dimension fr(2*ngrid),f(2*ngrid)
      equivalence (fro(1,1),f(1)),(ro(1,1),fr(1))
      no=n0
      do 1 i=1,ngrid
        ro(i,no)=0.
        roe(i,no)=0.
   1    fro(i,no)=0.
      n=nr(no)
  700 if(no.eq.2)go to 50
      uu=2./3.
      g=sqrt(uu)
      ri(1)=0.2
      ri(2)=0.6
      ri(3)=0.9
      ri(4)=1.4
      ri(5)=1.9
      ri(6)=2.3
      ri(7)=2.6
      ri(8)=3.1
      ri(9)=3.5
      ri(10)=4.2
      ri(11)=4.9
      ri(12)=5.2
      sq(1)=0.034724
      sq(2)=0.430761
      sq(3)=0.203166
      sq(4)=0.192986
      sq(5)=0.083866
      sq(6)=0.033007
      sq(7)=0.014201
      sq(8)=0.
      sq(9)=0.006860
      sq(10)=0.
      sq(11)=0.000438
      sq(12)=0.
      do 11 i=1,12
   11 a(i)=4.*sq(i)/(2.*pi32*g**3*(1.+2.*(ri(i)/g)**2))
      r=0.
      do 3 j=1,nu
      s=0.
      do 22 i=1,12
   22 s=s+a(i)*(exp(-((r-ri(i))/g)**2)+exp(-((r+ri(i))/g)**2))
      ro(j,1)=s
      r=r+d1
    3 continue
      r=0.
      roe(1,1)=0.
      ro(1,1)=0.
      do 80 i=2,n
        rr=r*r
        roe(i,1)=pi4*rr*ro(i,no)*exp(-ro(i,no)*beta)
        ro(i,1)=ro(i,no)*rr*pi4
   80 r=r+d1
      return
   50 continue
      return
      end
c
      subroutine ftrans(bcof,lm,npoint,r,h,p,bint,b,ndim)
c performs Fourier transform
      implicit double precision (a-h,o-z)
      dimension bcof(ndim,lm),bint(lm,3,3),f(2),bfac(2),ll(2),
     1 fi(3),b(lm),bdif(3)
      data eps,eps1,eps4 /1.d-3,1.d-13,1.d-1/
      l1=lm-1
      l2=lm-2
      d=h*p
      dflm=dfac(l1)
      if(lm.gt.1) dflm1=dfac(l2)
      xeps=0.
      if(l1.gt.0) xeps=(eps*dflm)**(1./real(l1))
      if(abs(d).lt.eps4) goto 300
      d2=d+d
      dd=d*d
      xfac=1./(p*dd)
      x=r*p
      n=(npoint+1)/2
      do 5 in=1,n
      xx=x*x
      m=mod(in,3)+1
      m1=mod(in-1,3)+1
      m2=mod(in-2,3)+1
      if(x.eq.0) goto 8
      if(lm.eq.1.or.abs(x).ge.xeps) goto 100
      if(abs(x).gt.eps1) goto 11
    8 do 9 j=1,3
      do 9 i=1,lm
    9 bint(i,j,m)=0.
      goto 200
   11 do 10 i=1,2
        ll(i)=lm-i
        f(i)=1.
   10 bfac(i)=1.
      do 15 i2=1,3
   15 fi(i2)=1./(ll(1)+i2)
      ind=0
   20 ind=ind+1
      do 30 i=1,2
        bfac(i)=-bfac(i)/(4.*ind*(ll(i)+ind+0.5))*xx
   30 f(i)=f(i)+bfac(i)
      do 35 i2=1,3
   35 fi(i2)=fi(i2)+bfac(1)/(ll(1)+i2+ind+ind)
      if(abs(bfac(2)).gt.eps1) goto 20
      xex=x**l2*xfac
      b(lm)=xex*x*f(1)/dflm
      b(l1)=xex  *f(2)/dflm1
      do 60 i2=1,3
   60 bint(lm,i2,m)=xex*x**(i2+1)*fi(i2)/dflm
      bint(l1, 1,m)=(bint(lm,2,m)+x*b(l1))/l1
      bint(l1, 2,m)=l1*bint(lm,1,m)+x*b(lm)
      bint(l1, 3,m)=l2*bint(lm,2,m)+xx*b(lm)
      if(lm.le.2) goto 200
      li=l1
      do 40 l=1,l2
        li=li-1
        b(li)=(li+li+1)/x*b(li+1)-b(li+2)
        bint(li,1,m)=((li+1)*bint(li+2,1,m)+(li+li+1)
     &   *b(li+1))/li
        bint(li, 2,m)=li*bint(li+1,1,m)+x*b(li+1)
   40 bint(li,3,m)= (li-1)*bint(li+1,2,m)+xx*b(li+1)
      goto 200
  100 s=sin(x)
      c=cos(x)
      b(1)=s/x*xfac
      bint(1,1,m)=sinint(x)*xfac
      bint(1,2,m)=(-c+1.)*xfac
      bint(1,3,m)=(s-x*c)*xfac
      if(lm.eq.1) goto 200
      b(2)=(b(1)-c*xfac)/x
      bint(2,1,m)=-b(1)+xfac
      bint(2,2,m)=bint(1,1,m)-s*xfac
      bint(2,3,m)=(-2.*c-x*s+2.)*xfac
      if(lm.eq.2) goto 200
      do 110 l=3,lm
        b(l)=(l+l-3)/x*b(l-1)-b(l-2)
        bint(l,1,m)=((l-2)*bint(l-2,1,m)-(l+l-3)*b(l-1))/(l-1)
        bint(l,2,m)= (l-1)*bint(l-1,1,m)-x*b(l-1)
  110 bint(l,3,m)= l*bint(l-1,2,m)-xx*b(l-1)
  200 continue
      if(in.eq.1) goto 5
      x1=x-d
      x2=x-d2
      do 230 l=1,lm
        do 210 i2=1, 3
  210   bdif(i2)=bint(l,i2,m)-bint(l,i2,m1)
        bcof(2*in-2,l)=-bdif(3)+2.*x1*bdif(2)
     +    -(x1*x1-dd)*bdif(1)
        if(in.gt.2) goto 240
        bcof(1,l)=0.5*bdif(3)-(x2+1.5*d)*bdif(2)
     +    +((.5*x2+1.5*d)*x2+dd)*bdif(1)
        goto 230
  240   bcof(2*in-3,l)=.5*(bint(l,3,m)-bint(l,3,m2))
     1    -x2*(bint(l,2,m)-bint(l,2,m2))
     2    -1.5*d*(bint(l,2,m)-2.*bint(l,2,m1)+bint(l,2,m2))
     3    +(.5*x2*x2+dd)*(bint(l,1,m)-bint(l,1,m2))
     4    +1.5*x2*d*(bint(l,1,m)-2.*bint(l,1,m1)+bint(l,1,m2))
        if(in.eq.n) bcof(npoint,l)=.5*bdif(3)-(x-1.5*d)
     1    *bdif(2)+((.5*x-1.5*d)*x+dd)*bdif(1)
  230 continue
    5 x=x+d2
      return
  300 continue
      xfac=h/3.
      x=p*r
      ivf=4
      do 340 in=1,npoint
        if(x.eq.0.0) goto 301
        if(lm.eq.1.or.abs(x).ge.xeps) goto 400
  301   ivf=6-ivf
        if(abs(x).ge.eps1) goto 360
        bcof(in,1)=ivf*xfac
        if(lm.eq.1) goto 340
        do 370 i1=2,lm
  370   bcof(in,i1)=0.
        goto 340
  360   xx=x*x
        do 310 i=1,2
          ll(i)=lm-i
          f(i)=1.
  310   bfac(i)=1.
        ind=0
  320   ind=ind+1
        do 330 i=1,2
          bfac(i)=-bfac(i)/(4.*ind*(ll(i)+ind+.5))*xx
  330   f(i)=f(i)+bfac(i)
        if(abs(bfac(2)).gt.eps1) goto 320
        xex=x**l2*xfac
        bcof(in,lm)=xex*x*f(1)*ivf/dflm
        bcof(in,l1)=xex  *f(2)*ivf/dflm1
        if(lm.le.2) goto  340
        li=l1
        do 345 l=1,l2
          li=li-1
  345   bcof(in,li)=(li+li+1)/x*bcof(in,li+1)-bcof(in,li+2)
  340 x=x+d
      goto 410
  400 n1=in
      do 351 in=n1,npoint
        ivf=6-ivf
        s=sin(x)*xfac
        c=cos(x)*xfac
        bcof(in,1)=s/x*ivf
        if(lm.eq.1) goto 351
        bcof(in,2)=(bcof(in,1)-c*ivf)/x
        if(lm.eq.2) goto 351
        do 350 l=3,lm
  350   bcof(in,l)=(l+l-3)/x*bcof(in,l-1)-bcof(in,l-2)
  351 x=x+d
  410 do 420 l=1,lm
        bcof(1,l)=.5*bcof(1,l)
  420 bcof(npoint,l)=.5*bcof(npoint,l)
      return
      end
c
      function dfac(l)
      implicit double precision (a-h,o-z)
      dfac=1.
      if(l.eq.0) return
      x=1.
      do 10 i=1,l
        x=x+2.
   10 dfac=dfac*x
      return
      end
c
      function sinint(x)
      implicit double precision (a-h,o-z)
      real*8 x,sinint
      if(abs(x).ge.16) goto 1
      z=x/16.
      y=4.*z**2-2.
      b=       0.00000 00000 00007d0
      a= y*b  -0.00000 00000 00185d0
      b= y*a-b+0.00000 00000 04185d0
      a= y*b-a-0.00000 00000 84710d0
      b= y*a-b+0.00000 00015 22370d0
      a= y*b-a-0.00000 00241 00076d0
      b= y*a-b+0.00000 03329 88589d0
      a= y*b-a-0.00000 39729 08746d0
      b= y*a-b+0.00004 04202 71419d0
      a= y*b-a-0.00034 54691 55569d0
      b= y*a-b+0.00243 62214 04749d0
      a= y*b-a-0.01386 74455 89417d0
      b= y*a-b+0.06203 36794 32003d0
      a= y*b-a-0.21126 37809 76555d0
      b= y*a-b+0.53014 88479 16522d0
      a= y*b-a-0.96832 22369 87086d0
      b= y*a-b+1.38930 87711 71888d0
      a= y*b-a-1.92656 50911 50656d0
      b= y*a-b+2.77875 63817 42663d0
      a= y*b-a-4.06398 08449 11986d0
      a= y*a-b+8.10585 29553 61245d0
      sinint=z*.5*(a-b)
      return
    1 z=16./x
      g=z**2
      y=4.*g-2.
      b=       0.00000 00000 00002d0
      a= y*b  -0.00000 00000 00014d0
      b= y*a-b+0.00000 00000 00107d0
      a= y*b-a-0.00000 00000 00964d0
      b= y*a-b+0.00000 00000 10308d0
      a= y*b-a-0.00000 00001 36096d0
      b= y*a-b+0.00000 00023 56196d0
      a= y*b-a-0.00000 00586 70317d0
      b= y*a-b+0.00000 24537 55677d0
      a= y*b-a-0.00023 37560 41393d0
      a= y*a-b+0.12452 74580 57854d0
      f= z*.5*(a-b)
      b=       0.00000 00000 00002d0
      a= y*b  -0.00000 00000 00012d0
      b= y*a-b+0.00000 00000 00087d0
      a= y*b-a-0.00000 00000 00717d0
      b= y*a-b+0.00000 00000 06875d0
      a= y*b-a-0.00000 00000 79604d0
      b= y*a-b+0.00000 00011 69202d0
      a= y*b-a-0.00000 00234 68225d0
      b= y*a-b+0.00000 07249 95950d0
      a= y*b-a-0.00004 26441 82622d0
      a= y*a-b+0.00772 57121 93407d0
      g= g*.5*(a-b)
      b= 1.57079 63267 94897d0
      if (x.lt.0.) b=-b
      sinint=b-f*cos(x)-g*sin(x)
      return
      end
c
      subroutine mahaux(at,elab,eref,expj0,as,es,v5,rw,aw,dwv,dws,dvolj,
     &  efermi)
      implicit double precision(a-h,o-z)
      real at,elab,eref,expj0,as,es,v5,rw,aw,dwv,dws,dvolj,efermi
      dimension ee(400),aji(400),ww(400),wd(400)
c
      pi=4.*atan(1.)
      a13=at**0.3333333
      ebeg=0.1
c calculate  volume and surface imaginary pots on energy grid ee
      do i=1,150
       ee(i)=ebeg+dble(i-1)
c damping term
       C=0.005
       if(ee(i).lt.13.) C=-0.165*ee(i)+2.15
c add damping term to surface potential
       v5f=v5*exp(-C*abs(ee(i)-efermi))
       te=1./(1.+exp(-(ee(i)-es)/as))
       aji(i)=-expj0*te
       wamp=-expj0*te/((1.-v5f)*rw**3/3.+7.603*v5f*aw*rw**2/a13)/pi
       ww(i)=wamp*(1.-v5f)
       wd(i)=-wamp*v5f
      enddo
c calculate integrals for dispersive corrections according to Mahaux, Ngo and Satchler, NPA 449 (1986) 354
c and volume integral
      dvef=0.
      dverf=0.
      dvefd=0.
      dverfd=0.
      dvvol=0.
      do j=1,149
       dj = ee(j + 1) - ee(j)
       an1 = (elab - ee(j))*log(abs((elab-ee(j))/dj))
       an2 = (elab - ee(j + 1))*log(abs((elab-ee(j+1))/dj))
       dve = ((ww(j+1) - ww(j))/dj)*(an1 - an2)/pi
       dved = ((wd(j+1) - wd(j))/dj)*(an1 - an2)/pi
       dva = ((aji(j+1) - aji(j))/dj)*(an1 - an2)/pi
       anr1 = (eref - ee(j))*log(abs((eref-ee(j))/dj))
       anr2 = (eref - ee(j + 1))*log(abs((eref-ee(j+1))/dj))
       dver = ((ww(j+1) - ww(j))/dj)*(anr1 - anr2)/pi
       dverd = ((wd(j+1) - wd(j))/dj)*(anr1 - anr2)/pi
       dvef = dvef + dve
       dverf = dverf + dver
       dvefd= dvefd + dved
       dverfd= dverfd + dverd
       dvvol =  dvvol + dva
      enddo
      dwv = -dvef
      dws = -dvefd
      dvolj = dvvol
c
      return
      end

