      function racah(a,b,c,d,e,f,g,numfac)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire (adapted from O. Bersillon)
c | Date  : November 22, 2006
c | Task  : Calculation of Racah coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numfac,ja,jb,jc,jd,je,jf,i(16),n,j,k,il,j1,j2,
     +                 j3,j4,j5,j6,j7
      real             racah,a,b,c,d,e,f,g(numfac),eps,q,p,r,o,v,w,x,y,
     +                 z,t
      double precision h,s
c
c *********** Calculation of Racah coefficients   w(a,b,c,d;e,f) *******
c
c Attention:   w(a,b,c,d;e,f) = (-1)**(a+b+c+d)*6-J(a,b,e;d,c,f)
c
c from John.G.Wills  ORNL-TM-1949 (August 1967)
c                    and Comp.Phys.Comm. 2(1971)381
c     O.Bersillon    August 1977
c
      eps  = 1.0e-03
      racah = 0.
c
c  Convert arguments to integer and make useful combinations
c
      ja  = int(2.*a + eps)
      jb  = int(2.*b + eps)
      jc  = int(2.*c + eps)
      jd  = int(2.*d + eps)
      je  = int(2.*e + eps)
      jf  = int(2.*f + eps)
      i(1)  = ja + jb - je
      i(2)  = jb + je - ja
      i(3)  = je + ja - jb
      i(4)  = jc + jd - je
      i(5)  = jd + je - jc
      i(6)  = je + jc - jd
      i(7)  = ja + jc - jf
      i(8)  = jc + jf - ja
      i(9)  = jf + ja - jc
      i(10) = jb + jd - jf
      i(11) = jd + jf - jb
      i(12) = jf + jb - jd
      i(13) = ja + jb + je
      i(14) = jc + jd + je
      i(15) = ja + jc + jf
      i(16) = jb + jd + jf
c
c Check triangular inequalities, find no. of terms in sum,
c divide I's by 2
c
      n = i(16)
      do j=1,12
        k = i(j)/2
        if(i(j) .ne. 2*k) return
        if(k .lt. 0) return
        if(k .lt. n) then
          n = k
        endif
        i(j) = k + 1
      end do
c
c Find minimum value of summation index
c
      il = 0
      do j=13,16
        i(j) = i(j)/2
        if(il .lt. i(j)) then
          il = i(j)
        endif
      end do
      j1 = il  - i(13) + 1
      j2 = il  - i(14) + 1
      j3 = il  - i(15) + 1
      j4 = il  - i(16) + 1
      j5 = i(13) + i(4)  - il
      j6 = i(15) + i(5)  - il
      j7 = i(16) + i(6)  - il
      h  = - exp(0.5*(g(i(1))+g(i(2))+g(i(3))-g(i(13)+2)+g(i(4))+g(i(5))
     ++g(i(6))-g(i(14)+2)+g(i(7))+g(i(8))+g(i(9))-g(i(15)+2)+g(i(10))+
     +g(i(11))+g(i(12))-g(i(16)+2))+g(il+2)-g(j1)-g(j2)-g(j3)-g(j4)-
     +g(j5)-g(j6)-g(j7))
      if((j5 - 2*(j5/2)) .ne. 0) h = -h
      if(n .lt. 0) return
      if(n .eq. 0) then
        racah = real(h)
        return
      else
        s = 1.
        q = n  - 1
        p = il + 2
        r = j1
        o = j2
        v = j3
        w = j4
        x = j5 - 1
        y = j6 - 1
        z = j7 - 1
        do j=1,n
          t = (p+q)/(r+q)*(x-q)/(o+q)*(y-q)/(v+q)*(z-q)/(w+q)
          s = 1. - s*t
          q = q   - 1.
        end do
        racah = real(h*s)
      endif
      return
      end
