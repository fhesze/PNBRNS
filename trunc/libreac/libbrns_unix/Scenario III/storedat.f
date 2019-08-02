c    $Id: storedat.f 16 2007-10-18 12:36:47Z centler $
c
c     SUBROUTINE storedat
c
       subroutine storedat()
         include 'common_geo.inc'
         include 'common.inc'
         include 'common_opt.inc'
         include 'common_meas.inc'
         integer t,j,i
         character*30 file_opt_name(ntopt)

c***********************************************************************
c      STOREDAT: read the measured data                                *
c                                                                      *
c      the filenames are obtained from the subroutine getdat           *
c      these files must have the following structure:                  *
c      1st line: time when measurements were taken                     *
c      2nd line: number of species measured at this time               *
c      3rd line: index number of the species                           *
c      4th line: number of measurements for the species                *
c      5th ff. : depth of measurements, measured concentration         *
c                                                 & standard deviation *
c      note that if std.dev<=0, a relative error is assumed (relsig)   *
c       together with a minimum error (stdmin, both defined in main.f) *
c      continue as 3rd ff.                                             *
c                                                                      *
c      example:                                                        *
c      2 species at time 1.5, nr.3 with 6 and nr.7 with 3 measurements *
c      1.5                                                             *
c      2                                                               *
c      3                                                               *
c      6                                                               *
c      0    0.01  0.                                                   *
c      0.1  0.02  0.                                                   *
c      0.2  0.04  0.                                                   *
c      0.3  0.06  0.                                                   *
c      0.45 0.07  0.                                                   *
c      0.6  0.05  0.                                                   *
c      7                                                               *
c      3                                                               *
c      0.1  0.08  0.                                                   *
c      0.15 0.1   0.                                                   *
c      0.2  0.024 0.                                                   *
c                                                                      *
c      ntopt: number of times measurements are availabe                *
c      timemeas: time of measurements                                  *
c      nrspmeas: number of species measured at a given time            *
c      idspmeas: index of what was measured                            *
c      nrxmeas: number of points in a given profile                    *
c      xmeas: depth of measurements                                    *
c      spmeas: concentration measured                                  *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************
c
c***********************************************************************
c      GETDATA: communication to maple, called from storedat.f         *
c      contains                                                        *
c      - index of parameters to be optimized                           *
c        (i.e. the identification number in the entire parameter list  *
c      - filenames from where data is read                             *
c        for required file content see storedat.f                      *
c     idpar(...)= xxx index of parameters to optimize                  *
c     e.g. idpar(1) = 12: parameter 12 gets optimized                  *
c          idpar(2) = 15: parameter 15 gets optimized                  *
c          for sequence of parameters see transferXX.f                 *
c     file_opt_name(...) = 'some name'                                 *
c                                                                      *
c     note:                                                            *
c     - list of files to read from include the extension in the name   *
c     - one file per time of measurement, containing all the species   *
c     e.g. file_opt_name(1) = 'meas_time1.dat'                         *
c                                                                      *
c     CM, Dec 2001                                                     *
c***********************************************************************
         call getdat(file_opt_name)

         do t=1,ntopt
           open(unit=2,file=file_opt_name(t),status='old')
           read(2,*) timemeas(t)
           read(2,*) nrspmeas(t)

           do j=1,nrspmeas(t)
               read(2,*) idspmeas(t,j)
               read(2,*) nrxmeas(t,j)
               do i=1,nrxmeas(t,j)
                 read(2,*) xmeas(t,j,i),spmeas(t,j,i),sigmeas(t,j,i)
c      write(*,*) xmeas(t,j,i),spmeas(t,j,i),sigmeas(t,j,i)
                end do
           end do

           close(2)
        end do
c     pause

       end
