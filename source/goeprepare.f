      subroutine goeprepare(tjl,numtjl,st,nweip,nweis,nweit,p,s,
     +  t,wp,ws,wt,tav,x02,x102,x202,x2i,x12i,x22i,fpst1,fpst2,s1,s2,
     +  s3,s4,s5,numtr,numinc)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : January 24, 2007
c | Task  : Preparation of GOE triple integral width fluctuation
c |         correction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numtjl,nweip,nweis,nweit,numtr,numinc,i,j,k
      real             p(nweip),s(nweis),t(nweit),wp(nweip),ws(nweis),
     +                 wt(nweit),x02(nweip,nweis,nweit),
     +                 x102(nweip,nweis,nweit),x202(nweip,nweis,nweit),
     +                 x2i(nweip,nweis,nweit),x12i(nweip,nweis,nweit),
     +                 x22i(nweip,nweis,nweit),fpst1(nweip,nweis,nweit),
     +                 fpst2(nweip,nweis,nweit),s1,s2,s3,s4,s5,p1,pp,p2,
     +                 d,ds2,wp1,wpspp2,ex02,ex2i,s21,s22,ums21,ums22,
     +                 ps21,umps21,pm2ps21,um2s21,pums212,uppm2ps2,
     +                 pm2s22,uppm2s2,pms222,ws2,ws1,e,ume,umes2,upes2,
     +                 wpws1,wpws2,t1,t2,wt2,t21,t22,s2t22,s2t21,umt21,
     +                 umt22,ps2t21,x1rat,x2rat
      double precision tjl(0:5,numtr),st,tav(numtr),tgamma,prodm,prodp
c
c *************************** GOE calculation **************************
c
c Initialisation and average transmission coefficients
c
c tgamma: gamma transmission coefficient
c numtjl: number of transmission coefficients
c tjl   : transmission coefficients
c tav   : average transmission coefficients
c
      tgamma=0.
      do 10 i=1,numtr
        tav(i)=0.
   10 continue
      do 20 i=1,numtjl
        if (tjl(0,i).gt.0.) tav(i)=tjl(1,i)/tjl(0,i)
   20 continue
      tgamma=tjl(1,numtjl+1)
c
c ****** Numerical calculation of triple integral for few channels *****
c
c st: denominator of compound nucleus formula
c
      if (st.lt.20.) then
c
c Loop over p
c
c p1,...: help variables
c
        do 30 i=1,nweip
          p1=p(i)+1.
          pp=0.5*p1
          p2=(3.-pp)/pp
          d=sqrt(0.5*p2)
          ds2=0.5*d
          wp1=2.*wp(i)
          wpspp2=3.*wp(i)/(pp*pp)
c
c Special case for capture
c
c exc02,...: help variables
c
          if (tgamma.eq.0.) then
            ex02=1.
            ex2i=1.
          else
            ex02=real(exp(-p1*tgamma*0.5))
            ex2i=real(exp(-p2*tgamma*0.5))
          endif
c
c Loop over s
c
c s1,...: help variables
c
          do 40 j=1,nweis
            s1=0.25*sqrt(2.)*(s(j)+1.)
            s2=ds2*s(j)+ds2
            s21=real(s1*s1)
            s22=real(s2*s2)
            ums21=1.-s21
            ums22=1.-s22
            ps21=p1*s21
            umps21=1.-ps21
            pm2ps21=p1-ps21-ps21
            um2s21=ums21-s21
            pums212=p1*ums21*ums21
            uppm2ps2=pm2ps21+1.
            pm2s22=p2-s22-s22
            uppm2s2=1.+pm2s22
            pms222=(p2-s22)*(p2-s22)
            ws2=d*ws(j)
            ws1=0.5*sqrt(2.)*ws(j)
            if (s2.gt.1.) then
              e=sqrt(1.-1./s22)
            else
              e=0.
            endif
            ume=1.-e
            umes2=0.5*ume
            upes2=0.5+0.5*e
            wpws1=wp1*ws1
            wpws2=wpspp2*ws2
c
c Loop over t
c
c t1,t2,...x2rat: help variables
c x02,x102,.....: variables for final GOE calculation
c fpst1,fpst2   : variables for final GOE calculation
c prodm,prodp   : product functions for GOE
c s1,s2,s3,s4,s5: variables for final GOE calculation
c
            do 50 k=1,nweit
              t1=0.5*t(k)+0.5
              t2=umes2*t(k)+upes2
              wt2=ume*wt(k)
              t21=t1*t1
              t22=t2*t2
              s2t22=s22*t22
              s2t21=s21*t21
              umt21=1.-t21
              umt22=1.-t22
              ps2t21=ps21*t21
              x02(i,j,k)=ps21-ps2t21
              x102(i,j,k)=ps2t21
              x202(i,j,k)=pm2ps21+ps2t21
              x2i(i,j,k)=s22-s2t22
              x12i(i,j,k)=s2t22
              x22i(i,j,k)=pm2s22+s2t22
              x1rat=real(prodm(x02(i,j,k),tjl,numtjl,numtr,numinc)/
     +          prodp(x102(i,j,k),tjl,numtjl,numtr,numinc)/
     +          prodp(x202(i,j,k),tjl,numtjl,numtr,numinc))
              fpst1(i,j,k)=wpws1*wt(k)*(ps2t21+umps21)*um2s21*umt21*
     +          ex02*x1rat**10/pums212/
     +          sqrt((um2s21+s2t21)*(1.+ps2t21)*(ps2t21+uppm2ps2))
              x2rat=real(prodm(x2i(i,j,k),tjl,numtjl,numtr,numinc)/
     +          prodp(x12i(i,j,k),tjl,numtjl,numtr,numinc)/
     +          prodp(x22i(i,j,k),tjl,numtjl,numtr,numinc))
              fpst2(i,j,k)=wpws2*wt2*umt22*(ums22+s2t22)*pm2s22*
     +          ex2i*x2rat**10/pms222/
     +          sqrt((1.+s2t22)*(pm2s22+s2t22)*(uppm2s2+s2t22))
   50       continue
   40     continue
   30   continue
      else
c
c ***** Numerical calculation of triple integral for many channels *****
c
        s1=0.
        s2=0.
        s3=0.
        s4=0.
        s5=0.
        do 110 i=numinc+1,numtjl+1
          s1=s1+real(tjl(1,i))
          s2=s2+real(tjl(2,i))
          s3=s3+real(tjl(3,i))
          s4=s4+real(tjl(4,i))
          s5=s5+real(tjl(5,i))
  110   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
