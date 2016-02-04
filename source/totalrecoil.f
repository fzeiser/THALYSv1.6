      subroutine totalrecoil
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : July 25, 2005
c | Task  : Recoil results
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen,iang,Zcomp,Ncomp,ixl,nexrec,iyl
      real    sumen,sum,dErecoil
c
c ******************** Angle-integrated spectra ************************
c
c flaglabddx : flag for calculation of DDX in LAB system
c sumen,sum  : help variables
c iejlab     : number of ejectile lab bins
c parskip    : logical to skip outgoing particle
c nanglecont : number of angles for continuum
c ddxejlab   : LAB ddx array
c dcosangcont: width of cosine bin
c xsejlab    : LAB ejectile spectrum
c twopi      : 2.*pi
c dEejlab    : width of ejectile lab bin
c xsejlabint : LAB energy-integrated spectrum
c
      if (flaglabddx) then
        do 10 type=0,6
          if (parskip(type)) goto 10
          sumen=0.
          do 20 nen=1,iejlab(type)
            sum=0.
            do 30 iang=0,nanglecont
              sum=sum+ddxejlab(type,nen,iang)*dcosangcont(iang)
   30       continue
            xsejlab(type,nen)=sum*twopi
            sumen=sumen+sum*twopi*dEejlab(type,nen)
   20     continue
          xsejlabint(type)=sumen
   10   continue
      endif
c
c ************************* Recoils in LAB system **********************
c
c Zcomp     : charge number index for compound nucleus
c maxZ      : maximal number of protons away from the initial compound
c             nucleus
c Ncomp     : neutron number index for compound nucleus
c maxN      : maximal number of neutrons away from the initial compound
c             nucleus
c xspopnuc  : population cross section per nucleus
c xseps     : limit for cross sections
c maxenrec  : number of recoil energies
c dErecoil  : width of recoil bin
c Erecmax   : maximal energy limit of recoil bin
c Erecmin   : minimal energy limit of recoil bin
c Nlast     : last discrete level
c ddxrec    : array containing the lab double differential xs of the
c             residual nucleus with (Zinit-i) protons and (Ninit-j)
c             neutrons, with excitation energy given by Ex(i,j,k) in
c             the lab bin with energy limits Erecmin(i,j,l),
c             Erecmax(i,j,l) and angular limits angrecmin(m),
c             angrecmax(m)
c areareclab: Total surface of LAB ddx bins
c specrecoil: recoil spectrum
c recoilint : total recoil integrated over spectrum
c
      do 110 Zcomp=0,maxZ
        do 110 Ncomp=0,maxN
          if (xspopnuc(Zcomp,Ncomp).lt.xseps) goto 110
          sumen=0.
          do 120 ixl=0,maxenrec
            dErecoil=abs(Erecmax(Zcomp,Ncomp,ixl)-
     +        Erecmin(Zcomp,Ncomp,ixl))
            if (dErecoil.le.1.0e-14) dErecoil=1.
            sum=0.
            do 130 nexrec=0,max(Nlast(Zcomp,Ncomp,0),1)
              do 130 iyl=0,nanglerec
                sum=sum+ddxrec(Zcomp,Ncomp,nexrec,ixl,iyl)*
     +            areareclab(Zcomp,Ncomp,ixl,iyl)
  130       continue
            specrecoil(Zcomp,Ncomp,ixl)=sum*twopi/dErecoil
            sumen=sumen+sum*twopi
  120     continue
          recoilint(Zcomp,Ncomp)=sumen
  110 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
