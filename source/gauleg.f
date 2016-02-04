      subroutine gauleg(ngl,tgl,wgl)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : July 9, 2004
c | Task  : Calculation of Gauss-Legendre arrays
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer ngl,ns2,i,it,k,j
      real    tgl(ngl),wgl(ngl),pi,ti,oi,p2,p1,p
c
c ****************** Calculation of Gauss-Legendre arrays **************
c
c ngl: number of points for Gauss-Legendre integration
c tgl: points for Gauss-Legendre integration
c wgl: weights for Gauss-Legendre integration
c
      data pi /3.14159265358979323/
      ns2=ngl/2
      do 10 i=1,ns2
        ti=(4*i-1)*pi/(4*ngl+2.)
        oi=cos(ti+1./(8.*ngl*ngl*tan(ti)))
        do 20 it=1,10
          p2=1.
          p1=oi
          do 30 k=2,ngl
            p=((2*k-1)*oi*p1-(k-1)*p2)/float(k)
            p2=p1
            p1=p
   30     continue
          oi=oi-p*(1.-oi*oi)/(ngl*(p2-oi*p))
   20   continue
        tgl(i)=oi
        wgl(i)=(1.-oi*oi)/(ngl*p2*ngl*p2)
   10 continue
      do 40 j=ns2+1,ngl
        tgl(j)=-tgl(j-ns2)
        wgl(j)=wgl(j-ns2)
   40 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
