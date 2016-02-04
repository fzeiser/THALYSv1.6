      subroutine getkeywords(line,word)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 30, 2006
c | Task  : Retrieve keywords and values from input line
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      character*1  chprev,ch
      character*80 line,word(40)
      integer      i,nkey,ibeg,iend
c
c **************** Read keywords and values from input line ************
c
c line     : input line
c word     : words on input line
c chprev,ch: character
c nkey     : word counter
c ibeg     : index to mark begin of word
c iend     : index to mark end of word
c
c From each input line we retrieve the keyword, the number of
c values, and the values themselves. These are all stored in strings
c (the word array).
c
      do 10 i=1,40
        word(i)='                                                  '
   10 continue
      chprev=' '
      nkey=0
      do 20 i=1,80
        if (i.gt.1) chprev=line(i-1:i-1)
        ch=line(i:i)
        if (ch.ne.' '.and.chprev.eq.' ') then
          nkey=nkey+1
          ibeg=i
        endif
        if (ch.eq.' '.and.chprev.ne.' ') then
          iend=i-1
          word(nkey)=line(ibeg:iend)
        endif
   20 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
