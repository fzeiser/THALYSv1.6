      function fmin(extern,k,psh,pd,pa,pe,ncall,eps)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : February 26, 2005
c | Task  : Searches the minimal value fmin of the function extern.
c |         If eps positive the Tchebychev accuracy within the p is
c |         controlled, if it is negative, the standard deviation of
c |         the function values is examined.
c +---------------------------------------------------------------------
c
c *************************** Comments ********************************
c
c This function is based on the function fmin originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      integer   mmax,nmaxp1,icall,ncall,i,n1,jk,k,minf,maxf
      parameter (mmax=15,nmaxp1=16)
      real      eps,fmin,extern,psh(k),pd(k),pa(k),pe(k),pr(mmax),
     +          pex(mmax),pc(mmax),gam,betap,alpha,fn,fmax,fmax2,pmax,
     +          sum,sum2,fr,fex,std,f(nmaxp1),ssh(mmax,nmaxp1)
      data alpha/1./,betap/.5/,gam/2./
c
c **********************************************************************
c
      if (k.le.mmax) goto 4
      write (*,5) k,mmax
 5    format(/1x,"from fmin: introduced number of parameters",i4,
     + "    larger than allowed number",i3/)
      stop
 4    icall=0
      do 2 i=1,k
        if (psh(i).lt.pa(i)) psh(i)=pa(i)
        if (psh(i).gt.pe(i)) psh(i)=pe(i)
 2    continue
      fn=float(k)
      n1=k+1

c construct initial simplex.
c
      do 10  jk=1,n1
        do 10  i=1,k
          ssh(i,jk)=psh(i)
 10   continue
      do 11  i=1,k
        ssh(i,i+1)=ssh(i,i+1) + pd(i)
 11   continue
      do 16 jk=1,n1
        if (jk.eq.1) go to 14
        do 12 i=1,k
          psh(i)=ssh(i,jk)
 12     continue
 14     icall=icall+1
        if (icall.gt.ncall) go to 80
        f(jk)=extern(psh,k)
 16   continue
c
c super loop.
c find largest, second largest, and smallest function value.
c
 17   fmax=f(1)
      maxf=1
      do 18  i=2,n1
        if (f(i).lt.fmax)  goto 18
        maxf=i
        fmax=f(i)
 18   continue
      fmin=f(1)
      minf=1
      do 19  i=2,n1
        if (f(i).gt.fmin)  goto 19
        minf=i
        fmin=f(i)
 19   continue
      fmax2=f(minf)
      do 20  i=1,n1
        if (i.eq.maxf)  goto 20
        if (f(i).lt.fmax2)  goto 20
        fmax2=f(i)
 20   continue
c
c determine centroid.
c
      do 22  i=1,k
        pc(i)=0.
 22   continue
      do 24  jk=1,n1
        if (jk.eq.maxf)  goto 24
        do 26  i=1,k
          pc(i)=pc(i)+ssh(i,jk)
 26     continue
 24   continue
      do 25  i=1,k
        pc(i)=pc(i)/fn
 25   continue
c
c     check accuracy.
c
      if (eps.lt.0.) go to 74
      std=0.
      do 72 i=1,k
        pmax=abs(pc(i)-ssh(i,maxf))
 72   continue
      if (pmax.gt.std) std=pmax
      go to 80
 74   sum=0.
      do 76 jk=1,n1
        sum=sum+f(jk)
 76   continue
      sum=sum/float(n1)
      sum2=0.
      do 78 jk=1,n1
        sum2=sum2+(f(jk)-sum)**2
 78   continue
      std=sqrt(sum2/fn)
 80   if (icall.lt.ncall.and.std.gt.abs(eps)) go to 27
      do 82 i=1,k
        psh(i)=ssh(i,minf)
 82   continue
      return
c
c reflect.
c
 27   do 28 i=1,k
        pr(i)=(1.+alpha)*pc(i)-alpha*ssh(i,maxf)
 28   continue
      do 29 i=1,k
        if (pr(i).lt.pa(i).or.pr(i).gt.pe(i)) go to 30
 29   continue
      icall=icall+1
      if (icall.gt.ncall) go to 80
      fr=extern(pr,k)
      go to 31
 30   icall=icall+1
      if (icall.gt.ncall) go to 80
      fr=abs(fmax)*(1.+abs((pr(i)-pc(i))/pd(i)))
 31   if (fr.lt.fmin)  goto 40
      if (fr.gt.fmax2)  goto 50
      goto 44
c
c expand.
c
 40   do 42 i=1,k
        pex(i)=gam*pr(i)+(1.-gam)*pc(i)
 42   continue
      do 32 i=1,k
        if (pex(i).lt.pa(i).or.pex(i).gt.pe(i)) go to 33
 32   continue
      icall=icall+1
      if (icall.gt.ncall) go to 80
      fex=extern(pex,k)
      go to 34
 33   icall=icall+1
      if (icall.gt.ncall) go to 80
      fex=abs(fmax)*(1.+abs((pex(i)-pc(i))/pd(i)))
 34   if (fex.gt.fmin)  goto 44
 45   do 46  i=1,k
        ssh(i,maxf)=pex(i)
 46   continue
      f(maxf)=fex
      goto 17
 44   do 48  i=1,k
        ssh(i,maxf)=pr(i)
 48   continue
      f(maxf)=fr
      goto 17
 50   if(fr.gt.fmax)  goto 54
      do 52  i=1,k
        ssh(i,maxf)=pr(i)
 52   continue
      f(maxf)=fr
c
c contract.
c
 54   do 56  i=1,k
        pex(i)=betap*ssh(i,maxf)+(1.-betap)*pc(i)
 56   continue
      icall=icall+1
      if (icall.gt.ncall) go to 80
      fex=extern(pex,k)
      if (fex.lt.fmax)  goto 45
c
c halve.
c
      do 62  jk=1,n1
        if (jk.eq.minf)  goto 62
        do 64  i=1,k
          ssh(i,jk)=.5*(ssh(i,jk)+ssh(i,minf))
          psh(i)=ssh(i,jk)
 64     continue
        icall=icall+1
        if (icall.gt.ncall) go to 80
        f(jk)=extern(psh,k)
 62   continue
      go to 17
      end
