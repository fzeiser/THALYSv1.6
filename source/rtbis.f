      function rtbis(func,x1,x2,xacc)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning (adapted from Numerical Recipes)
c | Date  : July 5, 2004
c | Task  : Search for zero crossings of the function
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit  none
      integer   jmax
      parameter (jmax=40)
      integer   j
      real      rtbis,func,x1,x2,xacc,fmid,f,dx,xmid
c
c ************************** Search ************************************
c
c func: function for which zero crossing is searched
c
      fmid=func(x2)
      f=func(x1)
      if (f.lt.0.) then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 10 j=1,jmax
        dx=dx*0.5
        xmid=rtbis+dx
        fmid=func(xmid)
        if (fmid.le.0.) rtbis=xmid
        if (abs(dx).lt.xacc.or.fmid.eq.0.) return
   10 continue
      end
