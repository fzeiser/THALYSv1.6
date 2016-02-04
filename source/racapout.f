      subroutine racapout
c
c +---------------------------------------------------------------------
c | Author: Yi Xu, Stephane Goriely, Arjan Koning
c | Date  : December 22, 2012
c | Task  : Output for racap into racap.out and racap.tot files
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nracapout,i,nen,JJJ
      integer Acpracap,Zcpracap,partar
      real    ecm,spintar
      character*50 racapfile
      character*1 cpar
c
c************************************************************************
c
      nracapout=nin
      Acpracap=Atarget+parA(k0)
      Zcpracap=ZZ(0,0,0)
      ecm=real(einc*specmass(parZ(k0),parN(k0),k0))

ctest racapgsspin=jdis(0,0,0)
ctest racapgsprty=parlev(0,0,0)
      spintar=jdis(parZ(k0),parN(k0),0)
      partar=parlev(parZ(k0),parN(k0),0)
      cpar='+'
      if (partar.eq.-1) cpar='-'

      if (nin.eq.1) then
        open (unit=2,status='unknown',file='racap.tot')
        write(2,'("Direct capture reaction on target ",i3,a2," (Z =",
     +    i3,") with projectile ",a1, " (Qvalue=",f7.3," MeV)"/)')
     +    Atarget,Starget,Ztarget,parsym(k0),Q(0)

        write(2,'("Characteristics of the Direct Radiative Capture",
     +   " calculation:")')
        if (ldmodelracap.eq.1) write(2,'("  ldmodelracap=1: spin-, ",
     +    "parity-dependent pphh NLDs used for racap calculation")')
        if (ldmodelracap.eq.2) write(2,'("  ldmodelracap=2: spin-, ",
     +    "parity-independent pphh NLDs used for racap calculation")')
        if (ldmodelracap.eq.3) write(2,'("  ldmodelracap=3: spin-, ",
     +    "parity-dependent total NLDs used for racap calculation ",
     +    "(ldmodel=5)")')
        if (racopt.eq.1) write(2,'("  racopt=1: Woods-Saxon potential",
     +    " KD used for racap calculation")')
        if (racopt.eq.3) write(2,'("  racopt=3: JLMB potential",
     +    " used for racap calculation")')

        write(2,'("  Total number of transitions from",f5.1,a1,
     +    " GS to ",i3,a2," experimental levels:",i3)') spintar,
     +    cpar,Acpracap,nuc(Zcpracap),nlevexpracap
        write(2,'("  Total number of transitions from",f5.1,a1,
     +    " GS to ",i3,a2," levels:",i3)') spintar,cpar,Acpracap,
     +    nuc(Zcpracap),maxex(0,0)-nlevexpracap+1
        write(2,'(/"  Spectroscopic factors ",/)')
        do i=0,numex
          if (i.eq.0.or.edis(0,0,i).gt.0.) then
            write(2,'(1p,g12.4,0p,f5.1,1x,a2,1p,g12.4,0p)')
     +        edis(0,0,i),jdis(0,0,i),cparity(parlev(0,0,i)),
     +        spectfac(0,0,i)
          endif
        enddo
        write(2,*)
      else
        open (unit=2,status='unknown',file='racap.tot')
   10   read(2,*,end=20,err=20)
        goto 10
      endif

   20 backspace 2
      write(2,*)
      write(2,'("==========  Direct capture at Elab=",1p,e10.3,
     +  "  ==========")') einc
      write(2,'("   Ecm =",1p,e10.3," MeV: ",/,
     +  "   Direct   radiative capture xs =",e12.5," mb",
     +  "  (Discrete=",e10.3," - Continuum=",e10.3,")",/,
     +  "   HF+Preeq radiative capture xs =",e12.5," mb",/,
     +  "   Total    radiative capture xs =",e12.5," mb")')
     +  ecm,xsracape,xsracapedisc,xsracapecont,
     +  xspopnuc(0,0)-xsracape,xspopnuc(0,0)

      write(2,*)

      if (mod(Acpracap,2).eq.0) then
        write(2,1222) (real(JJJ),JJJ=0,10),(real(JJJ),JJJ=0,10)
 1222   format('Nlvl',2x,' Ex  ',2x,' J ','   Pi ',2x,' Sp ',2x,
     +  'J/p=tot',3x,11('Jp=',f4.1,'+',1x),11('Jp=',f4.1,'-',1x))
      else
        write(2,1223) (real(JJJ+0.5),JJJ=0,10),(real(JJJ+0.5),JJJ=0,10)
 1223   format('Nlvl',2x,' Ex  ',2x,' J ','   Pi ',2x,' Sp ',2x,
     +  'J/p=tot',3x,11('Jp=',f4.1,'+',1x),11('Jp=',f4.1,'-',1x))
      endif

      do i=0,maxex(0,0)
           if (i.le.nlevexpracap-1) then
            write(2,1225) i,edis(0,0,i),jdis(0,0,i),parlev(0,0,i),
     &      spectfac(0,0,i),xsracappopex(i)
          else
            write(2,1226) i,Ex(0,0,i),spectfac(0,0,i),
     &        xsracappopex(i),(xsracappop(i,JJJ,1),JJJ=0,10),
     &        (xsracappop(i,JJJ,-1),JJJ=0,10)
          endif
      enddo
 1225 format(1x,i3,1x,f6.3,1x,f5.1,1x,i3,1x,f6.2,1x,1p,e10.3,32e9.2)
 1226 format(1x,i3,1x,f6.3,1x,' 99.9',1x,'  0',1x,f6.2,1x,1p,e10.3,
     &  32e9.2)
      close(2)
c
c write output racap.out with summary of reaction cross section
c
      racapfile='racap.out'
      if (nin.eq.1) then
        open(unit=1,status='unknown',file=racapfile)
        write(1,'("# ",a1," + ",i3,a2,": Direct Capture to ",i3,a2)')
     +    parsym(k0),Atarget,Starget,Acpracap,nuc(Zcpracap)
        write(1,'("# Q-value    =",1p,e12.5,0p," mass=",f11.6,
     +    " Emax=",f11.6)') Qres(0,0,0),nucmass(0,0),
     +    min(S(0,0,1),S(0,0,2))
        write(1,'("# # transitions from ",f5.1,a1,
     +    " GS to ",i3,a2," levels:",i3)') spintar,cpar,Acpracap,
     +    nuc(Zcpracap),nlevracap(0,0)
        write(1,'("# # energies =",i3)') nracapout
        write(1,'("#    E         xs                    ",
     +    " xs(E1)      xs(E2)      xs(M1)     xs(tot)")')
        close(1)
      endif
      open(unit=1,status='unknown',file=racapfile)
      do 30 nen=1,nin+4
        read(1,*,end=40,err=40)
   30 continue
      write(1,'(1p,e10.3,e12.5,12x,4e12.5)') eninc(nin),xsracap(nin),
     +  xsracapEM(nin,1,1), xsracapEM(nin,1,2),
     +  xsracapEM(nin,0,1),xspopnuc(0,0)
   40 close(1)
c
      return
      end
