      function rpoint(z)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : rpoint=0 gives for a given nucleon number hidden in the
c |         left-hand side of the dinuclear complex the rupture point
c |         zriss
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This function is based on the function rpoint originally developed by
c U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real rpoint,z,vr2
c
c **********************************************************************
c
      rpoint=di*vr2(z1,z2,z)+rest
      return
      end
