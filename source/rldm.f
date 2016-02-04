      subroutine rldm(iz,ia,il,egs,esp)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : December 13, 2013
c | Task  : Saddle point energies, rotating gs energy
c +---------------------------------------------------------------------
c
c *************************** Comments ********************************
c
c This subroutine returns the fission barrier height in MeV. It is
c based on calculations using the rotating liquid drop model (S. Cohen,
c F. Plasil, and W.J. Swiatecki, Ann. of Phys. 82, 577 (1974)). This
c subroutine is based on the subroutine FISROT incorporated in ALICE-91
c (M. Blann, presented at the "Workshop on Computation and Analysis of
c Nuclear Data Relevant to Nuclear Energy and Safety", Trieste, Italy,
c 1992.)
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ia,iz,il,ix,iy
      real    amass,zchar,neut,ll,paren,eso,ero,x,y,egs,esp,cx,bx,dx,
     +        by,cy,dy,h1,h2,hf,b1,b2,bfrldm
c
c ******************** Rotating Liquid Drop Model *********************
c
c iz,zchar: charge number of residual nucleus
c ia,amass: mass number of residual nucleus
c neut    : neutron number of residual nucleus
c il,ll   : angular momentum
c cx,bx...: help variables
c x,y,....: help variables
c paren   : help variables
c x1b,x2b.: fitted parameters
c eso     : surface energy
c egs     : rotating ground state energy
c esp     : saddle point energy
c
      zchar=real(iz)
      amass=real(ia)
      ll=real(il)
      neut=amass-zchar
      paren=1.-1.7826*((neut-zchar)/amass)**2
      eso=17.9439*paren*amass**0.666666
      x=0.019655*zchar*(zchar/amass)/paren
      ero=34.548*ll*ll/amass**1.666666
      y=1.9254*ll*ll/(paren*amass**2.3333333)
      ix=int(20.*x+1.)
      cx=real(ix)
      bx=int(20.*x+1.)
      dx=real(bx)-cx
      if(x-.25)5,5,30
    5 by=10.*y+1.
      if(by-9.)15,15,10
   10 by=9.
   15 if(by-1.)20,20,25
   20 by=1.
   25 iy=int(by)
      cy=iy
      dy=by-cy
      h1=(x1h(ix+1,iy)-x1h(ix,iy))*dx+x1h(ix,iy)
      h2=(x1h(ix+1,iy+1)-x1h(ix,iy+1))*dx+x1h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      b2=(x1b(ix+1,iy+1 )-x1b(ix,iy+1))*dx+x1b(ix,iy+1)
      b1=(x1b(ix+1,iy)-x1b(ix,iy))*dx+x1b(ix,iy)
      bfrldm=(b2-b1)*dy+b1
      go to 95
   30 if(x-.5)35,35,60
   35 by=20.*y+1.
      if(by-10.)45,45,40
   40 by=10.
   45 if(by-1.)50,50,55
   50 by=1.
   55 ix=ix-5
      iy=int(by)
      cy=iy
      dy=by-cy
      h1=(x2h(ix+1,iy)-x2h(ix,iy))*dx+x2h(ix,iy)
      h2=(x2h(ix+1,iy+1)-x2h(ix,iy+1))*dx+x2h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      b1=(x2b(ix+1,iy)-x2b(ix,iy))*dx+x2b(ix,iy)
      b2=(x2b(ix+1,iy+1 )-x2b(ix,iy+1))*dx+x2b(ix,iy+1)
      bfrldm=(b2-b1)*dy+b1
      go to 95
   60 if(x-.95)70,70,65
   65 x=.95
   70 ix=int(20.*x+1.)
      ix=ix-10
      by=100.*y+1.
      if(by-19.)80,80,75
   75 by=19.
   80 if(by-1.)85,85,90
   85 by=1.
   90 iy=int(by)
      cy=iy
      dy=by-cy
      ix=min(ix,9)
      iy=min(iy,19)
      h1=(x3h(ix+1,iy)-x3h(ix,iy))*dx+x3h(ix,iy)
      h2=(x3h(ix+1,iy+1)-x3h(ix,iy+1))*dx+x3h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      b1=(x3b(ix+1,iy)-x3b(ix,iy))*dx+x3b(ix,iy)
      b2=(x3b(ix+1,iy+1 )-x3b(ix,iy+1))*dx+x3b(ix,iy+1)
      bfrldm=(b2-b1)*dy+b1
   95 egs=ero+hf*eso
      esp=egs+bfrldm*eso
      return
      end
