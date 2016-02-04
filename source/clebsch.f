      function clebsch(aj1,aj2,aj3,am1,am2,am3,g,numfac)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire (adapted from O. Bersillon)
c | Date  : November 22, 2006
c | Task  : Calculation of Clebsch-Gordan coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numfac,j1,j2,j3,m1,m2,m3,i(11),n,j,k,il,la,lb,
     +                 l,m
      real             clebsch,aj1,aj2,aj3,am1,am2,am3,g(numfac),eps,a,
     +                 b,h,d,e,f,q,t,x
      double precision c,s
c
c ************** Calculation of Clebsch-Gordan coefficients ************
c
c Attention:  cg(j1,j2,j3,;m1,m2,m3) = (-1)**(j1+j2-m3)*
c                                      3-J(j1,j2,j3;m1,m2,-m3)
c
c from John.G. Wills ORNL-TM-1949 (August 1967)
c                    and Comp.Phys.Comm. 2(1971)381
c     O.Bersillon    August 1977
c
      eps   = 1.0e-03
      clebsch = 0.
c
c Convert the arguments to integer
c
      j1 = int(2.*aj1 + eps)
      j2 = int(2.*aj2 + eps)
      j3 = int(2.*aj3 + eps)
      m1 = int(2.*am1 + sign(eps,am1))
      m2 = int(2.*am2 + sign(eps,am2))
      m3 = int(2.*am3 + sign(eps,am3))
c
c Test m1 + m2 = m3
c
      if(m1 + m2 - m3 .ne. 0) return
c
c Test table size
c
      i(10) = (j1 + j2 + j3)/2 + 2
      n     = i(10)
      i(11) = j3 + 2
      if(i(10) .gt. numfac) return
      i(1) = j1 + j2 - j3
      i(2) = j2 + j3 - j1
      i(3) = j3 + j1 - j2
      i(4) = j1 - m1
      i(5) = j1 + m1
      i(6) = j2 - m2
      i(7) = j2 + m2
      i(8) = j3 - m3
      i(9) = j3 + m3
c
c Check i(j) = even, triangular inequality, m less than j,
c find number of terms
c
      do j=1,9
        k = i(j)/2
        if(i(j) .ne. 2*k) return
        if(k .lt. 0) return
        if(k .lt. n) then
          n = k
        endif
        i(j) = k + 1
      end do
      if(m3 .ne. 0 .or. m1 .ne. 0 .or. m1 .ne. 1) then
        il = 0
        la = i(1) - i(5)
        lb = i(1) - i(6)
        if(il .lt. la) then
          il = la
        endif
        if(il .lt. lb) then
          il = lb
        endif
c
c Form coefficients of sum
c
        c  = (g(i(11)) - g(i(11)-1) + g(i(1)) + g(i(2)) + g(i(3)) -
     +        g(i(10)) + g(i(4)) + g(i(5)) + g(i(6)) + g(i(7)) +
     +        g(i(8)) + g(i(9)))*0.5
        j1 = i(1) - il
        j2 = i(4) - il
        j3 = i(7) - il
        m1 = il + 1
        m2 = il - la + 1
        m3 = il - lb + 1
        c  = c - g(j1) - g(j2) - g(j3) - g(m1) - g(m2) - g(m3)
        c  = exp(c)
        if((il - 2*(il/2)) .ne. 0) c = -c
        if(n .lt. 0) return
        if(n .eq. 0) then
          clebsch = real(c)
          return
        else
c
c Form sum
c
          a = j1 - 1
          b = j2 - 1
          h = j3 - 1
          d = m1
          e = m2
          f = m3
          s = 1.
          q = n - 1
          do j=1,n
            t = (a-q)/(d+q)*(b-q)/(e+q)*(h-q)/(f+q)
            s = 1. - s*t
            q = q - 1.
          end do
          clebsch = real(c*s)
          return
        endif
      else
c
c Special formula for m3 = 0 and m1 = 0 or 1/2
c
        k = i(10)/2
        if(i(1) .eq. 2*k) then
          k = 0
        else
          k = 1
        endif
        if(m1 .eq. 0) then
          l = 0
          if(k .ne. 0) return
        else if(m1 .eq. 1) then
          l = 1
        endif
        x  = l
        m  = i(3) + (i(1) + k + 1)/2 - l
        m1 = i(10)/2 + k
        m2 = i(4) + i(5)
        m3 = i(6) + i(7)
        j1 = (i(1) + 1 - k    )/2
        j2 = (i(2) + 1 + k - l)/2
        j3 = (i(3) + 1 + k - l)/2
        clebsch = exp((g(i(11))-g(i(11)-1)+g(i(1))+g(i(2))+g(i(3))-
     +          g(i(10)))/2. + g(m1)-g(j1)-g(j2)-g(j3)+
     +          x*(g(3)-(g(m2)-g(m2-1)+g(m3)-g(m3-1))/2.))
        if((m - 2*(m/2)) .ne. 0) clebsch = -clebsch
        return
      endif
      end
