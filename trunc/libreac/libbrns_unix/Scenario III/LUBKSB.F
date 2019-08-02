c    $Id: LUBKSB.F 16 2007-10-18 12:36:47Z centler $
c  ****************LUBKSB.F********************      
       subroutine LUBKSB(aam,n,n0,indxx,bbm)
!      implicit double precision(a-h,o-z)

       implicit real*8 (a-h,o-z)
       dimension aam(n0,n0),bbm(n0),indxx(n)
c
       II = 0
       do 12 i = 1,n
          ll = indxx(i)
          sum = bbm(ll)
          bbm(ll) = bbm(i)      
          if (II .ne. 0) then
             do 11 j = II,i-1
                sum = sum - aam(i,j)*bbm(j)
   11       continue
          else if (sum .ne. 0.0d0) then
             II = i
          endif
          bbm(i) = sum
   12 continue
c
       do 14 i = n,1,-1
          sum = bbm(i)
          if (i .lt. n) then 
             do 13 j = i+1,n
                sum = sum - aam(i,j)*bbm(j)
   13       continue
          endif
          bbm(i) = sum/aam(i,i)
   14 continue
c
       return 
       end
c***********************************************************
