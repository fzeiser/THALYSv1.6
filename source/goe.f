      subroutine goe(tjl,numtjl,na,nb,st,nweip,nweis,nweit,tav,x02,x102,
     +  x202,x2i,x12i,x22i,fpst1,fpst2,s1,s2,s3,s4,s5,res,numtr,ielas)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : February 6, 2007
c | Task  : GOE triple integral width fluctuation correction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numtjl,na,nb,nweip,nweis,nweit,numtr,ielas,i,j,k
      real             x02(nweip,nweis,nweit),x102(nweip,nweis,nweit),
     +                 x202(nweip,nweis,nweit),x2i(nweip,nweis,nweit),
     +                 x12i(nweip,nweis,nweit),x22i(nweip,nweis,nweit),
     +                 fpst1(nweip,nweis,nweit),
     +                 fpst2(nweip,nweis,nweit),s1,s2,s3,s4,s5,res,
     +                 dab,y1,y2,y3,y4,y5,y6,y7,y8,res02,func1,res2i
      double precision tjl(0:5,numtr),st,tav(numtr),res1,ta,tb,
     +                 rhob,s12,s13,s14,s15,ta2,ta3,ta4,ta5,tb2,tb3,tb4,
     +                 tb5,tb6,taptb,taptb2,taptb3,taptb4,taptb5,c1s1,
     +                 c1s2,c1s3,c1s4,c2s1,c2s2,c2s3,c2s4,c2s5,tt
c
c *************************** GOE calculation **************************
c
c ielas   : designator for elastic channel
c tjl     : transmission coefficients
c numtjl  : number of transmission coefficients
c na,nb   : counters for width fluctuation calculation
c tav,ta..: average transmission coefficients
c rhob    : help variable
c st      : denominator of compound nucleus formula
c nweip,..: variables for GOE triple integral calculation
c res     : width fluctuation factor
c res1,dab: help variables
c st      : denominator of compound nucleus formula
c
c Initialization
c
      res1=0.
      dab=0.
      if (ielas.eq.1) dab=1.
c
c ****** Numerical calculation of triple integral for few channels *****
c
      if (st.lt.20.) then
        ta=tav(na)
        tb=tav(nb)
        rhob=tjl(1,nb)
        if (nb.eq.numtjl+1) tb=dble(0.)
c
c Loop over p
c
        do 10 i=1,nweip
c
c Loop over s
c
          do 20 j=1,nweis
c
c Loop over t
c
            do 30 k=1,nweit
c
c y1.......res2i: help variables
c fpst1,fpst2   : variables for final GOE calculation
c x02,x102,.....: variables for final GOE calculation
c func1         : function for GOE
c
              y1=fpst1(i,j,k)
              y2=fpst2(i,j,k)
              y3=x102(i,j,k)
              y4=x202(i,j,k)
              y5=x02(i,j,k)
              y6=x12i(i,j,k)
              y7=x22i(i,j,k)
              y8=x2i(i,j,k)
              res02=y1*func1(y3,y4,y5,ta,tb,dab)
              res2i=y2*func1(y6,y7,y8,ta,tb,dab)
              res1=res1+real(ta*(res02+res2i)*rhob)
   30       continue
   20     continue
   10   continue
      else
c
c ***** Numerical calculation of triple integral for many channels *****
c
c Compute sum of tc,tc**2,tc**3 ....
c
c s12,.......,tt: help variables
c s1,s2,s3,s4,s5: variables for final GOE calculation
c
        s12=s1*s1
        s13=s12*s1
        s14=s13*s1
        s15=s14*s1
        ta=tjl(1,na)
        ta2=ta*ta
        ta3=ta2*ta
        ta4=ta3*ta
        ta5=ta4*ta
        tb=tjl(1,nb)
        tb2=tjl(2,nb)
        tb3=tjl(3,nb)
        tb4=tjl(4,nb)
        tb5=tjl(5,nb)
        tb6=tb*tjl(5,nb)/max(tjl(0,nb),1.d0)
        taptb=ta*tb+tb2
        taptb2=ta2*tb+ta*tb2+tb3
        taptb3=ta3*tb+ta2*tb2+ta*tb3+tb4
        taptb4=ta4*tb+ta3*tb2+ta2*tb3+ta*tb4+tb5
        taptb5=ta5*tb+ta4*tb2+ta3*tb3+ta2*tb4+ta*tb5+tb6
        c1s1=(-2.-4.*ta)/s1
        c1s2=(6.+3.*s2+12.*ta+34.*ta2)/s12
        c1s3=(-32.-12.*s2-16.*s3-64.*ta-40.*ta*s2-136.*ta2-304.*ta3)/s13
        c1s4=(240.+80.*s2+25*s2*s2+80*s3+100.*s4+480.*ta+200.*ta*s2+
     +    240.*ta*s3+864.*ta2+524.*ta2*s2+1520.*ta3+3508.*ta4)/s14
        c2s1=-taptb/s1
        c2s2=(s2*tb-2.*taptb+4.*taptb2)/s12
        c2s3=(s2*(2.*tb-5.*taptb)-4.*s3*tb+4.*taptb+8.*taptb2-20.*
     +    taptb3)/s13
        c2s4=(s2*(-4.*tb-16.*taptb+36.*taptb2)+s3*(-8.*tb+24.*taptb)+
     +    20.*s4*tb-12.*taptb+5.*s2*s2*tb-20.*taptb2-68.*taptb3+148.*
     +    taptb4)/s14
        c2s5=(s2*(12.*tb+48.*taptb+156.*taptb2-336.*taptb3)-60.*s2*s3*
     +    tb-148.*s5*tb+s3*(20.*tb+108.*taptb-228.*taptb2)+s4*(68.*tb-
     +    168.*taptb)+s2*s2*(16.*tb-41.*taptb)+64.*taptb+88.*taptb2+
     +    204.*taptb3+608.*taptb4-1348.*taptb5)/s15
        res1=dab*2.*ta2*(1.-ta)/s12*(1.+c1s1+c1s2+c1s3+c1s4)+
     +    (1.+dab)*ta/s1*(tb+c2s1+c2s2+c2s3+c2s4+c2s5)
      endif
      tt=tjl(1,na)*tjl(1,nb)/tjl(0,na)
      if (tt.ne.0.) then
        res=real(res1*st/tt)
      else
        res=0.
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
