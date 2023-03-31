      subroutine dsacbnd(excl)
      implicit none
      include 'dsaccom.h'
      integer excl
      if (aclabel.eq.'default') then
         call dsacbnd10(excl)
      elseif (aclabel.eq.'pdg2009') then
         call dsacbnd9(excl)
      elseif (aclabel.eq.'pdg2002e') then
         call dsacbnd8(excl)
      elseif (aclabel.eq.'pdg2002d') then
         call dsacbnd7(excl)
      elseif (aclabel.eq.'pdg2002c') then
         call dsacbnd6(excl)
      elseif (aclabel.eq.'pdg2002b') then
         call dsacbnd5(excl)
      else if (aclabel.eq.'pdg2002') then
         call dsacbnd4(excl)
      else if (aclabel.eq.'pdg2000') then
         call dsacbnd3(excl)
      else if (aclabel.eq.'mar2000') then
         call dsacbnd2(excl)
      else if (aclabel.eq.'pdg1999') then
         call dsacbnd1(excl)
      else
         write(*,*) 'Error in dsacbnd: invalid option: ',aclabel
         stop
      endif
      return
      end
