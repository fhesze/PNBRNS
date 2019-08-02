c    $Id: printdepth.f 16 2007-10-18 12:36:47Z centler $
       subroutine printdepth()

        include 'common_geo.inc'
        include 'common.inc'

        open(unit=9,file='echogrid.dat',status='unknown')

        write(9,*) 'echo of grid & depth dependent properties'
       write(9,*) 'note that first and last depth are ghost points'
       write(9,*) 'x, por, A, vd, vs, disp(species)'
        do j=0,nx+1,1
        write(9,8) x(j),por(j),area(j),vd(j),vs(j),(disp(k,j),k=1,ncomp)
       end do

        close(9)

 8     format(200(e14.7,2x))

      end
