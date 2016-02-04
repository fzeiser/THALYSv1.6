      subroutine neck(Z,A,fmass,fmasscor,fmz,fmzcor,ap,edefo,elt,crel)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : June 11, 2007
c | Task  : Fission fragment mass yields per fission mode based on RNRM
c +---------------------------------------------------------------------
c
c *************************** Comments ********************************
c
c This subroutine is based on the subroutine neck originally developed
c by U. Brosa (version 28.7.89).
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Z,A,Zix,Nix,izmax,jimax,izloop,mcount,jmx,imax,i,k,
     +  izstepnum,iaf1,izf1,iaf2,izf2,irn1,irn2
      real sumw,expo,rhodi,eob,ezdisnorm,ezdis,coulel,vnel,rnma,rnmi,
     +     es1,es2,bind01,bind02,es,bind1,bind2,ze1,ze2,am1,am2,coul12,
     +     sform,s12,x1,x2,b1,b2,de,vr3,ve1,ve2,a1,a2,zriss,zo,zu,
     +     vr2,vr1,amin,fmin,fimin,delt,tmp,d,dum,
     +     bcom,at,gam,zda,zstepsize,astepsize,astepnum,r0,xnu,
     +     rayl,psh(7),pd(7),pa(7),pe(7),af1(4000),rn1(4000,numelem),
     +     rn2(4000,numelem),zdis(4000,numelem),af2(4000),
     +     zf2(4000,numelem),sumtmp(4000),wgt(4000),wlog(4000),
     +     fmass(nummass),rtbis,fmasscor(nummass),crel,aloop,elt,edefo,
     +     ap,ztot,atot,fmzcor(nummass,numelem),fmz(nummass,numelem),
     +     zf1(4000,numelem),ald,ignatyuk,dumm
      external fidi,rpoint,evap
      data r0,xnu,rayl/1.15,1.0,11.00/
c
c determine mass and charge grid (depending whether evaporation
c correction is required)
c
      Zix=Zinit-Z
      Nix=Ninit-(A-Z)
      atot=real(A)
      ztot=real(Z)
      if(flagffevap)then
         astepnum=10.
         astepsize=0.1
         izstepnum=30
         zstepsize=0.1
      else
         astepnum=1.
         astepsize=1.
         izstepnum=3
         zstepsize=1.
      endif
c
c di   : nucleon number density
c vtot : total volume of the complex
c rp,rt: radii of projectile and target
c c    : curvature of neck
c r    : ratio of neck contribution
c tmp  : temperature
c zda  : charge over mass ratio
c gam  : surface tension
c
      totl = (1.15/1.2249)*2*elt
      rayl =  11*elt/( 2.4*1.2249*((atot)**.333333) )
      di=3./(4.*pi*r0**3)
      vtot=atot/di
      zda=ztot/atot
      gam=.9517*(1.-1.7826*((atot-2.*ztot)/atot)**2)
      at=atot-ap
      rt=r0*at**(1./3.)
      rp=r0*ap**(1./3.)
c
c fission option
c
      call bdef(atot,ztot,0.,dum,dumm,bcom)
      d=totl-rt-rp
      r2=totl/rayl
      c0=crel*2.*(rt+rp+2.*(SQRT((rt-r2)*(rp-r2))-r2))/d**2
      rpt=(rp/rt)**xnu
      ald=ignatyuk(Zix,Nix,edefo,0)
      tmp=sqrt(edefo/ald)
c
c Starting values for fmin to compute the shape of the dinuclear system
c psh(1) for a, psh(2) for z2, psh(3,4) for z1,z3, psh(5,6) for r1,r3.
c psh(7) curvature at the smallest cross section of the neck.
c
      psh(1)=.5
      pa(1)=.0
      pe(1)=1.
      pd(1)=.1
      psh(2)=0.
      pa(2)=-.2
      pe(2)=.2
      pd(2)=.05
      do 30 i=3,4
      psh(i)=.65
      pa(i)=.3
      pe(i)=1.
 30   pd(i)=.35
      do 31 i=5,6
      psh(i)=.8
      pa(i)=.6
      pe(i)=1.
 31   pd(i)=.2
      psh(7)=1.
      pa(7)=.5
      pe(7)=10.
      pd(7)=.5
      r1=rp
      r3=rt
      z1=r1*psh(3)
      z3=r3*psh(4)
c
c fmin needs a few starts to find the right values.
c
      do 102 k=1,15
      if (k.LT.2) go to 100
      delt=.2**((k+1)/2)
      do 101 i=1,7
      pa(i)=psh(i)-delt
      pe(i)=psh(i)+delt
 101  pd(i)=delt
 100  fimin=fmin(fidi,7,psh,pd,pa,pe,200,-1.E-4)
 102  if (fimin.lt.1.E-3) goto 103
 103  d=totl-r1-r3
c
c graphical discussion of the rupture shape.
c
      amin=di*vr1(z1)
      imax=di*vr2(z1,z2,z3)*astepnum
c
c In this loop the properties of the different fragmentations are
c calculated, as tke(a,z), neutron number rn(a,z), and yield wgt(a),
c zdis(a,z).
c
      jmx=0
      do 3 i=0,imax
      aloop=amin+i*astepsize
c
c calculate the rupture cut at zriss depending on the mass number.
c
      rest=amin-aloop
      zu=z1
      zo=z3
      zriss=rtbis(rpoint,zu,zo,2.E-4)
c
c calculation of the equivalent ellipsoidal shapes and the
c Coulomb and nuclear proximity repulsion energies
c
      a1=.5*(zriss+r1)
      a2=.5*(d+r3-zriss)
      ve1=vr1(z1)+vr2(z1,z2,zriss)
      ve2=vr3(z3)+vr2(zriss,z2,z3)
      b1=sqrt(3.*ve1/(4.*pi*a1))
      b2=sqrt(3.*ve2/(4.*pi*a2))
      de=a1+a2
      x1=a1**2-b1**2
      if(x1.ge.0.)then
         x1=sqrt(x1)/de
      else
         x1=-sqrt(-x1)/de
      endif
      x2=a2**2-b2**2
      if(x2.ge.0.)then
         x2=sqrt(x2)/de
      else
         x2=-sqrt(-x2)/de
      endif
      s12=sform(x1,x2)
      coul12=1.44*s12/de
      x1=x1*de/a1
      x2=x2*de/a2
      am1=aloop
      am2=atot-am1
      jmx=jmx+1
      mcount=1
      do 3, izloop=1,2*izstepnum+1
      ze1=zda*am1+(izloop-izstepnum-1)*zstepsize
      ze2=ztot-ze1
      if(ze1.lt.0..or.ze2.lt.0.)goto 3
c
c excess internal energy of ruptures nucleus es. es1 and es2
c excitation energies of separated fragments, rn1 and rn2 number of
c evaporated neutrons from light and heavy fragments.
c
      call bdef(am1,ze1,x1,dum,dumm,bind1)
      call bdef(am2,ze2,x2,dum,dumm,bind2)
      es=edefo
      call bdef(am1,ze1,0.,dum,dumm,bind01)
      call bdef(am2,ze2,0.,dum,dumm,bind02)
      es1=es/atot*am1
      es2=es-es1
      es1=bind1-bind01+es1
      es2=bind2-bind02+es2
      wlog(jmx)=es1+es2
      amm=am1
      zee=ze1
      ess=es1
      rnmi=-1.
      rnma=wlog(jmx)*.167
      rn1(jmx,izloop)=0.
      if(rnma.gt.0.7) rn1(jmx,izloop)=rtbis(evap,rnmi,rnma,1.E-2)
      if(rn1(jmx,izloop).lt.0.)rn1(jmx,izloop)=0.
      amm=am2
      zee=ze2
      ess=es2
      rnmi=-1.
      rnma=wlog(jmx)*.167
      rn2(jmx,izloop)=0.
      if(rnma.gt.0.7) rn2(jmx,izloop)=rtbis(evap,rnmi,rnma,1.E-2)
      if(rn2(jmx,izloop).lt.0.)rn2(jmx,izloop)=0.
c
c charge distribution
c
      zdis(jmx,izloop)=0.
      vnel=4.*pi*gam*(b1*b2)**2/(a1*b2*b2+a2*b1*b1)*(-1.7817)
      coulel=coul12*ze1*ze2
      ezdis=(bind1+bind2+coulel+vnel)/tmp
      if(mcount.eq.1)ezdisnorm=ezdis
      mcount=mcount+1
      ezdis=ezdis-ezdisnorm
      if (abs(ezdis).le.80.) zdis(jmx,izloop)=exp(-ezdis)
c
c total kinetic energy tke, mass probability distribution wgt(a).
c
      eob=2.*pi*gam*(rhodi(zriss)**2-rhodi(z2)**2)
      wgt(jmx)=0.
      expo=eob/tmp
      if(expo.lt.80.) wgt(jmx)=exp(-expo)
      af1(jmx)=am1
      zf1(jmx,izloop)=ze1
      af2(jmx)=am2
      zf2(jmx,izloop)=ze2
 3    continue
      jimax=jmx
c
c mass and charge distributions, with(out) corrections for evaporated
c neutrons
c
      sumw=0.
      izmax=izstepnum*2+1
      do 71 k=1,jimax
         sumw=sumw+wgt(k)
         sumtmp(k)=0.
         do 72 i=1,izmax
            sumtmp(k)=sumtmp(k)+zdis(k,i)
 72      continue
 71   continue
      do 74 k=1, jimax
         wgt(k)=wgt(k)/sumw
         do 73, i=1,izmax
            zdis(k,i)=zdis(k,i)/sumtmp(k)
            iaf1=max(int(af1(k)+0.5),1)
            izf1=max(int(zf1(k,i)+0.5),1)
            iaf2=max(int(af2(k)+0.5),1)
            izf2=max(int(zf2(k,i)+0.5),1)
            fmz(iaf1,izf1)=wgt(k)*zdis(k,i)+fmz(iaf1,izf1)
            fmz(iaf2,izf2)=wgt(k)*zdis(k,i)+fmz(iaf2,izf2)
c
c calculate corrections for neutron evaporation if evapcor=true
c
            if(flagffevap)then
              irn1=max(int(af1(k)-rn1(k,i)+0.5),1)
              irn2=max(int(af2(k)-rn2(k,i)+0.5),1)
              fmzcor(irn1,izf1)=fmzcor(irn1,izf1)+wgt(k)*zdis(k,i)
              fmzcor(irn2,izf2)=fmzcor(irn2,izf2)+wgt(k)*zdis(k,i)
            endif
 73      continue
         fmass(int(af1(k)+0.5))=fmass(int(af1(k)+0.5))+wgt(k)
         fmass(int(af2(k)+0.5))=fmass(int(af2(k)+0.5))+wgt(k)
 74   continue
      do 76 k=1,nummass
         do 75 i=1,numelem
            fmasscor(k)=fmasscor(k)+fmzcor(k,i)
 75      continue
 76   continue
      end
