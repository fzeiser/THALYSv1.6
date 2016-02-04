      subroutine adjust(E,key,Zix,Nix,type,ibar,factor)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 2, 2012
c | Task  : Energy-dependent parameter adjustment
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      integer      Zix,Nix,type,ibar,i,k,lenkey,N,j,npow,nen
      real         E,factor,Ea,Eb,Em,D,Da,Db,C,Eadj(0:numenadj),
     +             Dadj(0:numenadj)
c
c ************************* OMP Wv function ****************************
c
c E         : incident energy
c key       : keyword
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c type      : particle type
c ibar      : fission barrier
c factor    : multiplication factor
c Nadjust   : number of adjustable parameters
c adjustkey : keyword for local adjustment
c adjustpar : local adjustment parameters
c adjustix  : local adjustment index
c adjustfile: file for local adjustment
c nenadjust : number of tabulated energies of local adjustment
c Eadjust   : tabulated energy of local adjustment
c Dadjust   : tabulated depth of local adjustment
c Ea        : start energy of local adjustment
c Eb        : end energy of local adjustment
c Em        : intermediate energy of local adjustment
c D         : depth of local adjustment
c
      factor=1.
      Eadj(0)=0.
      Dadj(0)=1.
      do 10 i=1,Nadjust
        do 20 k=1,80
          if (adjustkey(i)(k:k).eq.' ') then
            lenkey=k-1
            goto 30
          endif
   20   continue
   30   if (key(1:lenkey).eq.adjustkey(i)(1:lenkey).and.
     +    Zix.eq.adjustix(i,1).and.Nix.eq.adjustix(i,2).and.
     +    type.eq.adjustix(i,3).and.ibar.eq.adjustix(i,4)) then
          if (adjustfile(i)(1:1).ne.' ') then
            N=nenadjust(i)
            do 40 j=1,nenadjust(i)
              Eadj(j)=Eadjust(i,j)
              Dadj(j)=Dadjust(i,j)
   40       continue
            if (E.ge.Eadj(1).and.E.le.Eadj(N)) then
              call locate(Eadj,1,N,E,nen)
              Ea=Eadj(nen)
              Eb=Eadj(nen+1)
              Da=Dadj(nen)
              Db=Dadj(nen+1)
              call pol1(Ea,Eb,Da,Db,E,factor)
            endif
          else
            Ea=adjustpar(i,1)
            Eb=adjustpar(i,2)
            Em=adjustpar(i,3)
            D=adjustpar(i,4)
            if (E.le.Ea.or.E.ge.Eb) goto 10
            npow=4
            Db=(D-1.)*(1.+0.5**npow)
            if (E.le.Em) then
              C=Ea+0.5*(Em-Ea)
              factor=1.+Db*((E-Ea)**npow)/((E-Ea)**npow+(C-Ea)**npow)
            else
              C=Em+0.5*(Eb-Em)
              factor=1.+Db*((Eb-E)**npow)/((Eb-E)**npow+(C-Eb)**npow)
            endif
          endif
        endif
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
