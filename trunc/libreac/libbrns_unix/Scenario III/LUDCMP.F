c    $Id: LUDCMP.F 16 2007-10-18 12:36:47Z centler $
c ******************LUDCMP.F**************************
c      
       subroutine LUDCMP(aam,n,n0,indxx,d)
!      implicit double precision(a-h,o-z)
!     include 'common_geo.inc'
!     include 'common.inc'

       implicit real*8 (a-h,o-z)
       parameter (tiny=1.d-20)
       dimension aam(n0,n0),indxx(n0),vv(n0)
!     real*8 d,aam(n0,n0),indxx(n0),vv(n0)
!     integer n,n0


       D = 1.0D0
       do 12 i = 1,n
          aammax = 0.0d0
         do 11 j = 1,n
           If (dabs(aam(i,j)) .gt. aammax) aammax=dabs(aam(i,j))

!     if(aam(i,j).ne.mat(i,j))then
!     write(*,*) 'aam.ne.mat', i,j, aam(i,j),mat(i,j)
!     endif

   11   continue
         if (aammax .eq. 0.0d0) then
              write(*,*) 'Singular matrix', i, n, n0
!            write(*,*) (aam(i,j),j=i,i+ncomp)
!                write(*,*) (mat(i,j),j=i,i+ncomp)
!              write(*,*) (aam(i,j),j=(i-1)*17+1,i*17)
             pause ! stop
              
         endif
         vv(i) = 1.0d0/aammax
   12 continue
c
       do 19 j = 1,n
         do 14 i = 1,j-1
           sum = aam(i,j)
           do 13 k = 1,i-1
             sum = sum - aam(i,k)*aam(k,j)
   13     continue
           aam(i,j) = sum
   14   continue
         aammax = 0.0d0
         do 16 i = j,n
           sum = aam(i,j)
           do 15 k = 1,j-1
             sum = sum - aam(i,k)*aam(k,j)
   15     continue
           aam(i,j) = sum
           dum = vv(i)*dabs(sum)
           if (dum .ge. aammax) then
             imax = i
             aammax = dum
           endif
   16   continue
         if (j .ne. imax) then
           do 17 k = 1,n
              dum = aam(imax,k)
              aam(imax,k) = aam(j,k)
              aam(j,k) = dum
   17     continue
           d = -d
           vv(imax) = vv(j)
         endif
         indxx(j) = imax
         if (aam(j,j) .eq. 0.0d0) aam(j,j) = tiny
         if (j .ne. n) then
            dum = 1.0d0/aam(j,j)
            do 18 i = j+1,n
              aam(i,j) = aam(i,j)*dum
   18      continue
         endif
   19 continue
c
       return 
       end
c ********************************************************************
