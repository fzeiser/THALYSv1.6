      subroutine convert(i)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 4, 2011
c | Task  : Convert input line from upper case to lowercase
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 str
      integer      i,k
c
c ************** Convert uppercase to lowercase characters *************
c
c inline,str: input line
c
c For easy handling of all the input parameters, the whole input is
c converted to lowercase characters, with the exception of filenames or
c other character strings.
c
      str(1:80)=inline(i)(1:80)
      do 10 k=1,80
        if (inline(i)(k:k).ge.'A'.and.inline(i)(k:k).le.'Z')
     +    inline(i)(k:k)=char(ichar(inline(i)(k:k))+32)
  10  continue
      do 20 k=0,60
        if (inline(i)(k+1:k+7).eq.'energy ') then
          inline(i)(k+8:80)=str(k+8:80)
          return
        endif
        if (inline(i)(k+1:k+7).eq.'optmod ') then
          inline(i)(k+8:80)=str(k+8:80)
          return
        endif
        if (inline(i)(k+1:k+7).eq.'nulldev') then
          inline(i)(k+8:80)=str(k+8:80)
          return
        endif
        if (inline(i)(k+1:k+8).eq.'bestpath') then
          inline(i)(k+9:80)=str(k+9:80)
          return
        endif
        if (inline(i)(k+1:k+8).eq.'integral') then
          inline(i)(k+9:80)=str(k+9:80)
          return
        endif
        if (inline(i)(k+1:k+9).eq.'abundance') then
          inline(i)(k+10:80)=str(k+10:80)
          return
        endif
        if (inline(i)(k+1:k+9).eq.'levelfile') then
          inline(i)(k+10:80)=str(k+10:80)
          return
        endif
        if (inline(i)(k+1:k+9).eq.'strucpath') then
          inline(i)(k+10:80)=str(k+10:80)
          return
        endif
        if (inline(i)(k+1:k+10).eq.'class2file') then
          inline(i)(k+11:80)=str(k+11:80)
          return
        endif
        if (inline(i)(k+1:k+10).eq.'deformfile') then
          inline(i)(k+11:80)=str(k+11:80)
          return
        endif
        if (inline(i)(k+1:k+10).eq.'radialfile') then
          inline(i)(k+11:80)=str(k+11:80)
          return
        endif
        if (inline(i)(k+1:k+10).eq.'rescuefile') then
          inline(i)(k+11:80)=str(k+11:80)
          return
        endif
        if (inline(i)(k+1:k+11).eq.'hbtransfile') then
          inline(i)(k+12:80)=str(k+12:80)
          return
        endif
        if (inline(i)(k+1:k+11).eq.'optmodfilen') then
          inline(i)(k+12:80)=str(k+12:80)
          return
        endif
        if (inline(i)(k+1:k+11).eq.'optmodfilep') then
          inline(i)(k+12:80)=str(k+12:80)
          return
        endif
        if (inline(i)(k+1:k+13).eq.'ompenergyfile') then
          inline(i)(k+14:80)=str(k+14:80)
          return
        endif
   20 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
