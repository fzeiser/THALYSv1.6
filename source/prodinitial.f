      subroutine prodinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 5, 2013
c | Task  : Initialization of isotope production info
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer k
      real    rhotable(numelem),TT
c
c Default material densities
c
c rhotable,rhotarget: target material density
c Ztarget           : charge number of target nucleus
c
      data (rhotable(k),k=1,numelem) /8.99e-5, 1.79e-4,
     +  0.534, 1.85, 2.34, 2.267, 1.25e-3, 1.43e-3, 1.70e-3, 9.0e-3,
     +  0.97, 1.74, 2.7, 2.33, 1.82, 2.07, 3.21e-3, 1.78e-3, 0.86, 1.55,
     +  2.99, 4.54, 6.11, 7.19, 7.43, 7.87, 8.9, 8.9, 8.96, 7.13,
     +  5.91, 5.32, 5.72, 4.79, 3.12, 3.75e-3, 1.53, 2.64, 4.47, 6.51,
     +  8.57, 10.22, 11.5, 12.37, 12.41, 12.02, 10.5, 8.65, 7.31, 7.31,
     +  6.68, 6.24, 4.93, 5.9e-3, 1.87, 3.59, 6.15, 6.77, 6.77, 7.01,
     +  7.3, 7.52, 5.24, 7.9, 8.23, 8.55, 8.8, 9.07, 9.32, 6.9,
     +  9.84, 13.31, 16.65, 19.35, 21.04, 22.6, 22.4, 21.45,19.32,13.55,
     +  11.85, 11.35, 9.75, 9.3, 7.0, 9.73e-3, 1.87, 5.5, 10.07, 11.72,
     +  15.4, 18.95, 20.2, 19.84, 13.67, 13.5, 14.78, 15.1, 13.5, 13.5,
     +  13.5, 13.5, 13.5, 13.5, 13.5, 13.5, 13.5, 13.5, 13.5, 13.5/
      if (rhotarget.eq.-1.) rhotarget=rhotable(Ztarget)
c
c ******************* Irradiation and cooling times ********************
c
c Determine irradiation time in seconds
c
c minutesec: number of seconds in a minute
c hoursec  : number of seconds in an hour
c daysec   : number of seconds in a day
c yearsec  : number of seconds in a year
c Tirrad   : irradiation time
c unitT    : irradiation time unit (y,d,h,m,s)
c T        : irradiation time per unit
c
      Tir=0.
      do k=1,5
        if (unitTirrad(k).eq.'y') Tir=Tir+yearsec*Tirrad(k)
        if (unitTirrad(k).eq.'d') Tir=Tir+daysec*Tirrad(k)
        if (unitTirrad(k).eq.'h') Tir=Tir+hoursec*Tirrad(k)
        if (unitTirrad(k).eq.'m') Tir=Tir+minutesec*Tirrad(k)
        if (unitTirrad(k).eq.'s') Tir=Tir+Tirrad(k)
      enddo
c
c Transform back to year, day, hour, minute, second
c (enabling general input)
c
c TT: help variable
c
      TT=Tir
      Tirrad(1)=int(TT/yearsec)
      TT=TT-Tirrad(1)*yearsec
      Tirrad(2)=int(TT/daysec)
      TT=TT-Tirrad(2)*daysec
      Tirrad(3)=int(TT/hoursec)
      TT=TT-Tirrad(3)*hoursec
      Tirrad(4)=int(TT/minutesec)
      TT=TT-Tirrad(4)*minutesec
      Tirrad(5)=int(TT)
c
c Determine cool time in seconds
c
c Tcool,Tcoo: cooling time
c unitTcool : cooling time unit (y,d,h,m,s)
c Tcool     : cooling time per unit
c
      Tco=0.
      do k=1,5
        if (unitTcool(k).eq.'y') Tco=Tco+yearsec*Tcool(k)
        if (unitTcool(k).eq.'d') Tco=Tco+daysec*Tcool(k)
        if (unitTcool(k).eq.'h') Tco=Tco+hoursec*Tcool(k)
        if (unitTcool(k).eq.'m') Tco=Tco+minutesec*Tcool(k)
        if (unitTcool(k).eq.'s') Tco=Tco+Tcool(k)
      enddo
c
c Transform back to year, day, hour, minute, second
c (enabling general input)
c
      TT=Tco
      Tcool(1)=int(TT/yearsec)
      TT=TT-Tcool(1)*yearsec
      Tcool(2)=int(TT/daysec)
      TT=TT-Tcool(2)*daysec
      Tcool(3)=int(TT/hoursec)
      TT=TT-Tcool(3)*hoursec
      Tcool(4)=int(TT/minutesec)
      TT=TT-Tcool(4)*minutesec
      Tcool(5)=int(TT)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
