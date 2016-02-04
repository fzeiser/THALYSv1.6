      subroutine ompadjust(E,k)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 15, 2012
c | Task  : Local optical model parameter adjustment
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      integer      k
      real         E,factor
c
c ******************* Create adjustable OMP factors *********************
c
c E          : incident energy
c flagincadj : flag for OMP adjustment on incident channel also
c flaginvecis: logical for calculating inverse channel OMP
c ompadjustp : flag for local optical model parameter adjustment
c adjust     : subroutine for energy-dependent parameter adjustment
c factor     : multiplication factor
c Fv1,...    : factor
c v1adjust.. : adjustable factors for OMP (default 1.)
c
      if (ompadjustp(k).and.(flagincadj.or.flaginvecis)) then
        key='v1adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fv1=factor*v1adjust(k)
        key='v2adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fv2=factor*v2adjust(k)
        key='v3adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fv3=factor*v3adjust(k)
        key='v4adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fv4=factor*v4adjust(k)
        key='rvadjust'
        call adjust(E,key,0,0,k,0,factor)
        Frv=factor*rvadjust(k)
        key='avadjust'
        call adjust(E,key,0,0,k,0,factor)
        Fav=factor*avadjust(k)
        key='w1adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fw1=factor*w1adjust(k)
        key='w2adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fw2=factor*w2adjust(k)
        key='w3adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fw3=factor*w3adjust(k)
        key='w4adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fw4=factor*w4adjust(k)
        key='rwadjust'
        call adjust(E,key,0,0,k,0,factor)
        Frw=factor*rwadjust(k)
        key='awadjust'
        call adjust(E,key,0,0,k,0,factor)
        Faw=factor*awadjust(k)
        key='rvdadjust'
        call adjust(E,key,0,0,k,0,factor)
        Frvd=factor*rvdadjust(k)
        key='avdadjust'
        call adjust(E,key,0,0,k,0,factor)
        Favd=factor*avdadjust(k)
        key='d1adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fd1=factor*d1adjust(k)
        key='d2adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fd2=factor*d2adjust(k)
        key='d3adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fd3=factor*d3adjust(k)
        key='rwdadjust'
        call adjust(E,key,0,0,k,0,factor)
        Frwd=factor*rwdadjust(k)
        key='awdadjust'
        call adjust(E,key,0,0,k,0,factor)
        Fawd=factor*awdadjust(k)
        key='vso1adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fvso1=factor*vso1adjust(k)
        key='vso2adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fvso2=factor*vso2adjust(k)
        key='rvsoadjust'
        call adjust(E,key,0,0,k,0,factor)
        Frvso=factor*rvsoadjust(k)
        key='avsoadjust'
        call adjust(E,key,0,0,k,0,factor)
        Favso=factor*avsoadjust(k)
        key='wso1adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fwso1=factor*wso1adjust(k)
        key='wso2adjust'
        call adjust(E,key,0,0,k,0,factor)
        Fwso2=factor*wso2adjust(k)
        key='rwsoadjust'
        call adjust(E,key,0,0,k,0,factor)
        Frwso=factor*rwsoadjust(k)
        key='awsoadjust'
        call adjust(E,key,0,0,k,0,factor)
        Fawso=factor*awsoadjust(k)
        key='rcadjust'
        call adjust(E,key,0,0,k,0,factor)
        Frc=factor*rcadjust(k)
      else
        Fv1=1.
        Fv2=1.
        Fv3=1.
        Fv4=1.
        Frv=1.
        Fav=1.
        Fw1=1.
        Fw2=1.
        Fw3=1.
        Fw4=1.
        Frw=1.
        Faw=1.
        Frvd=1.
        Favd=1.
        Fd1=1.
        Fd2=1.
        Fd3=1.
        Frwd=1.
        Fawd=1.
        Fvso1=1.
        Fvso2=1.
        Frvso=1.
        Favso=1.
        Fwso1=1.
        Fwso2=1.
        Frwso=1.
        Fawso=1.
        Frc=1.
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
