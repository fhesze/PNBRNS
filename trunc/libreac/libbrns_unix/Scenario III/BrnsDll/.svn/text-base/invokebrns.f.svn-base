#include<defines.inc>
      
      subroutine brnsIsAlive
      
cDEC$ ATTRIBUTES ALIAS:'_brnsIsAlive', DLLEXPORT::BRNSISALIVE
      
      implicit none
      write(*,*) "BRNS is Alive!"
      
      end subroutine brnsIsAlive
      
      subroutine invokebrns(concAfterTransport,
     +                           concBeforeTransport,
     +                        outputConcentrations,
     +                           numberOfSpecies,
     +                           timeStep,
     +                        fixedConcentrationBoundary,
     +                           returnValue,
     +                        pos_x, pos_y, pos_z,
     +                        porosity, waterSaturation,
#ifdef RETURNRATES
     +                        parameterVector,
     +                        rateVector)
#else
     +                        parameterVector)
#endif

c     Format of Arrays: conc_species1, conc_species2, ...
      
cExpose subroutine invokebrns to users of this DLL
c
cDEC$ ATTRIBUTES ALIAS:'_invokebrns', DLLEXPORT::INVOKEBRNS

c     implicit none
      
      include '../common_geo.inc'
      include '../common.inc'
      include '../common_drive.inc'
      
      
c Variables
      
      real*8 x_pos
      real*8 y_pos
      real*8 z_pos

      real*8 porosity
      real*8 waterSaturation
      
      real*8 concAfterTransport(1)
      real*8 concBeforeTransport(1)
      real*8 outputConcentrations(1)
      integer numberOfSpecies
      real*8 timeStep
      integer fixedConcentrationBoundary(1)
      integer returnValue
      integer k
      integer l
      real*8 parameterVector(1)
#ifdef RETURNRATES
      real*8 rateVector(1)
#endif
c     integer ibct(ncomp,2)

      integer doRateOutput
      
      real*8 pos_x
      real*8 pos_y
      real*8 pos_z

	  skipChemistry = 0
      

c Body of invokebrns
c TODO: also implement newt() (?)
      
      
#ifdef DLLWRITELOG
c     open(unit=25, file='brnsdll.log', access='append', action='write')
      open(unit=25, file='brnsdll.log', access='append')
      write(25,*) 'Inside invokebrns(): numberOfSpecies=',
     +             numberOfSpecies,
     +             'timestepsize=', timeStep,
     +             'returnValue=', returnValue,
     +             'pos_x=', pos_x,
     +             'pos_y=', pos_y,
     +             'pos_z=', pos_z
      write(25,*) 'concAfterTransport:'
      do k = 1,ncomp
      write(25,*) concAfterTransport(k)
      end do
      write(25,*) 'concBeforeTransport:'
      do k = 1,ncomp
      write(25,*) concBeforeTransport(k)
      end do
#endif
      
      doRateOutput = returnValue

      returnValue = -1
      
      x_pos = pos_x;
      y_pos = pos_y;
      z_pos = pos_z;
      
      CALL drivervalues()  ! Setting up ilud, ftol, etol, newton, isol
c                          (for choosing between netwonsub and newt)
      
      call basic()
c     call molecular()      

c     in switches.f, the x coordinate will be fetched from array x:
      x(1) = x_pos
      call switches(1)
      call parameters(parameterVector)
      
      call biogeo()
      
c     call gridsetup()
c     call porarea()
c     call advdiffcoeff()
c     call printdepth()
      
      call boundaries()     
c     call initialcond()
      
c     call transferfw()
c     call storedat()
      
c Setting up Everything so that newtonsub() can be called

	  por0 = porosity

      IF (ncomp /= numberOfSpecies) THEN
      write (*,*) "In BrnsDll: Number of species mismatch! Expected ",
     +      ncomp, " but received ", numberOfSpecies, " species!"
      returnValue=111
      return
      END IF
      
      delt = timeStep

      do k = 1,ncomp
      sp(k,1)    = concBeforeTransport(k)
      if (sp(k,1).lt.1.0d-20) then
      sp(k,1)=1.0d-20
#ifdef DLLWRITELOG
      write(25,*) 'concBeforeTransport: species',k,'<1e-20, setting to ',
     +  'this value'
#endif
      endif
      spold(k,1) = concAfterTransport(k)
      if (spold(k,1).lt.1.0d-20) then
      spold(k,1)=1.0d-20
#ifdef DLLWRITELOG
      write(25,*) 'concAfterTransport: species',k,'<1e-20, setting to ',
     +       'this value'
#endif
      endif
c     set ibc != 0 so that boundaries are not considered in newtonsub.f
c        ibct(k,1) = ibc(k,1)    
c        ibct(k,2) = ibc(k,2)
      if ( fixedConcentrationBoundary(k) == 0 ) then
      ibc(k,1) = 5
      else
      ibc(k,1) = 0
      spb(k,1) = concAfterTransport(k)
      end if 
      ibc(k,2)  = 5
      end do 
      
c Calling newtonsub()
    
	  if ( skipChemistry == 0 ) then
      	CALL newtonsub(1,returnValue,ftol,etol,newton,ilud)
	  else
	    do k = 1,ncomp
		  sp(k,1) = spold(k,1)
		end do
		call rates(1)
		do k = 1,ncomp
			if (spb(k,2).ne.0) then
			 sp(k,1) = spold(k,1) + r(spb(k,2),1) * delt
			end if
		end do
	  end if

c Returning new concentration values
      
      do k = 1,ncomp
      call limits(k,1)
      outputConcentrations(k) = sp(k,1)
      end do 

	  call updateporosity(1)
	  porosity = por0

c     output rates
      if (doRateOutput < 0) then
c     this is the last timestep of the simulation
c     rates. calls switches.f; in switches.f, the x coordinate will be fetched from array x:
      x(1) = x_pos
      call rates(1)
      if (doRateOutput == -1) then
c     replace existing file
      open(unit=11,file='ratesAtFinish.dat',access='append')
      close(11,status='delete')
      endif
      open(unit=11,file='ratesAtFinish.dat',access='append')
      write(11,2000) x_pos, y_pos, z_pos, (r(l,1), l = 1, nreac)
2000  format(3f12.4, 512e14.7)      
      close(11)
      endif

#ifdef RETURNRATES
      x(1) = x_pos
      call rates(1)
      do k = 1, nreac
      rateVector(k) = r(k,1)
      end do
#endif
      
#ifdef DLLWRITELOG
      write(25,*) 'returning from invokebrns(): returnValue=',
     +             returnValue      
      write(25,*) 'concAfterTransport:'
      do k = 1,ncomp
      write(25,*) concAfterTransport(k)
      end do
      write(25,*) 'concBeforeTransport:'
      do k = 1,ncomp
      write(25,*) concBeforeTransport(k)
      end do
      write(25,*) 'outputConcentrations:'
      do k = 1,ncomp
      write(25,*) outputConcentrations(k)
      end do
      close(unit=25)
#endif
      
      end subroutine invokebrns

