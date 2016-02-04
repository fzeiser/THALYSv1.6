      subroutine timer
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 25, 2006
c | Task  : Output of execution time
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      integer hour,minute,second,hundred
      real    time,etime,tarray(2)
c
c ****** Get elapsed time in seconds from beginning of execution *******
c
c etime  : time function
c tarray : help variable
c hour   : number of hours
c minute : number of minutes
c second : number of seconds
c hundred: number of 1/100th of seconds
c
c The returned time should be "charge time" (e.g., cp+pp+sys). This
c could be machine dependent.
c
      time=etime(tarray)
      hour=int(time/3600.)
      minute=int((time-hour*3600)/60.)
      second=int(time-hour*3600-minute*60)
      hundred=int(100*(time-int(time)))
      write(*,'(/" Execution time:",i3," hours ",i2," minutes ",i2,".",
     +  i2.2," seconds ")') hour,minute,second,hundred
      write(*,'(/" The TALYS team congratulates you with this ",
     +  "successful calculation.")')
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
