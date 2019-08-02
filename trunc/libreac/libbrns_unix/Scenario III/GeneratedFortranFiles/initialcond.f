c      
c     SUBROUTINE initialcond
c      
      subroutine initialcond()
          include 'common_geo.inc'
          include 'common.inc'
          include 'common_drive.inc'
        real*8 spi(ncomp)
          if (ic.eq.1) then
          open(unit=85,file='initialconc.txt',status='old')
c           
          do 1000, j=1,nx,2
            read(85,*) de,(sp(i,j),i=1,ncomp) 
 1000     continue
c         
          close(85)
          endif
          if (ic.eq.2) then
              spi(1) = 0.1D1
              spi(2) = 0.1D1
              spi(3) = 0.1D1
              spi(4) = 0.1D1
              spi(5) = 0.1D1
c             
            do 1001, i=1,ncomp
c               
              do 1002, j=1,nx,2
                sp(i,j) = spi(i)
 1002         continue
c             
 1001       continue
c           
          endif
          if (ic.eq.3) then
              open(unit=101,file='dummy.inp',status='old')
              open(unit=102,file='dummy.inp',status='old')
              open(unit=103,file='dummy.inp',status='old')
              open(unit=104,file='dummy.inp',status='old')
              open(unit=105,file='dummy.inp',status='old')
c             
            do 1003, j=1,nx,2
                read(101,2000) sp(1,j),depth
 2000           format(1x,e14.7,2x,f8.4)
                read(102,2001) sp(2,j),depth
 2001           format(1x,e14.7,2x,f8.4)
                read(103,2002) sp(3,j),depth
 2002           format(1x,e14.7,2x,f8.4)
                read(104,2003) sp(4,j),depth
 2003           format(1x,e14.7,2x,f8.4)
                read(105,2004) sp(5,j),depth
 2004           format(1x,e14.7,2x,f8.4)
 1003       continue
c           
              close(101)
              close(102)
              close(103)
              close(104)
              close(105)
          endif
      end
