      function sform(x,y)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Form factor for the Coulomb interaction energy between two
c |         spheroids.
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function sform originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      integer k,m
      real    x,y,sform,z1loc,z2loc,zn,sn,zm,cnm,sm
c
c **********************************************************************
c
      if (abs(x).gt.abs(y)) then
c
c ph. quentin, j. de physique 30 (1969) 497.
c
        z1loc=x
        z2loc=y
      else
        z1loc=y
        z2loc=x
      endif
      if (z1loc.eq.0.) then
        sform=1.
        return
      endif
      z1loc=sign(z1loc**2,z1loc)
      z2loc=sign(z2loc**2,z2loc)
      sform=0.
      zn=1./z1loc
      do 20 k=0,20
        zn=zn*z1loc
        sn=0.75/((k+0.5)*(k+1.5))*zn
        zm=1.
        cnm=1.
        do 10 m=1,15
          zm=zm*z2loc
          cnm=cnm*(k+m)*(k+m-0.5)/(m*(m-0.5))
          sm=0.5625/((k+0.5)*(k+1.5)*(m+0.5)*(m+1.5))*cnm*zn*zm
          sn=sn+sm
 10     continue
        if (abs(sm).lt.1.e-6) go to 15
 15     sform=sform+sn
 20   continue
      if (abs(sn).lt.1.e-5) go to 25
 25   return
      end
