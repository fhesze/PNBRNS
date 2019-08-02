c    $Id: optimde.f 16 2007-10-18 12:36:47Z centler $
       subroutine optimde()
!  *********************************************************************
!  Differential Evolution for Optimal Control Problems                 *
!                                                                      *
!  This code is based on the DE code outlined in                       *
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of *
!  the ICEC'96 contest by differential evolution. IEEE Conf. on        *
!  Evolutionary Computation, 842-844.                                  *
!  for original source codes see                                       *
!  http://www.icsi.berkeley.edu/~storn/code.html                       *
!  The Fortran 90 program on which this Fortran 77 version is based    *
!  was written by Dr. Feng-Sheng Wang, Department of Chemical          *
!  Engineering, National Chung Cheng University, Chia-Yi 621, Taiwan   *
!  e-mail: chmfsw@ccunix.ccu.edu.tw                                    *
! **********************************************************************

       include 'common_geo.inc'
       include 'common.inc'
        include 'common_opt.inc'
       include 'common_meas.inc'
       include 'common_drive.inc'

      parameter (np=10*nopt)
       integer strategy, refresh ! intent in
       integer nfeval ! intent out
       real*8 VTR, CR_XC, F_XC, XCmin(nopt), XCmax(nopt) ! intent in
       real*8 bestmem_XC(nopt) ! intent inout
       real*8 bestval ! intent out
       real*8 pop_XC(np,nopt), bm_XC(np,nopt), mui_XC(np,nopt)
       real*8 mpo_XC,popold_XC(np,nopt), rand_XC(np,nopt),ui_XC(np,nopt)
       integer i, ibest, iter, k, ind(4), t
       integer rot(np), a1(np), a2(np), a3(np), a4(np), a5(np), rt(np)
       real*8 tempval, rand_C1(nopt), val(np), bestmemit_XC(nopt)
      real*8 popvec(nopt), rand2(np), funk
      integer num, jj, number, ilin
        real*8  rhelp(nopt)

       external funk
!!      intrinsic random_number, mod
       intrinsic mod

!      funk : The user provided file for evaluating the objective
!             function. originally subroutine obj(xc,fitness), where
!             "xc" is the real decision parameter vector.(input) and
!             "fitness" is the fitness value.(output)
!             changed to the function funk(xc), fitness = funkvalue
!      Dim_XC : Dimension of the real decision parameters, number of
!               parameters of the objective function = nopt
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
! this range should contain the global
! minimum (else algorithm doesn't work so well)

!   there is some argument if XCmin/max should simply be range of search
!      VTR : The expected fitness value to reach. should be related to
!      the error of measurements. here set to 0 but problim-stopping
!      criteria added (iof=1!)
!      NP : Population size - not so critical, must be >= 5
!           10*#parameters is good initial guess, defined in parameter
!           statement above
!      itermax : The maximum number of iterations (generations). set
!      to "itmax/np"
!      F_XC : Mutation scaling factor for real decision parameters [0,2]
!             A good initial guess is 0.8 (0.5-1)
!      CR_XC : Crossover factor for real decision parameters. [0,1].
!              For correlated paramters, high values of CR work better.
!              start at 0.5
!     strategy : The strategy of the mutation operations is used in HDE.
!           1=DE/best/1/bin: best & mutation from 2
!           2=DE/rand/1/bin: a node & mutation from 2
!           3=DE/rand/2/bin: a node & mutation from 2 plus best
!           4=DE/best/2/bin: best & mutation from 4
!           5=DE/rand/2/bin: a node & mutation from 4
!     refresh : The intermediate output will be produced after "refresh"
!     iterations. No intermediate output will be produced if refresh < 1
!      bestmen_XC(Dim_XC) : The best parameter set
!      bestval : The best objective function
!      nfeval : The number of function call
!
!      ilin: 1: uniform random selection, else random out of log(range)
!      omagn: orders of magnitude for range around the initial profile


       ! driver variables, adapt if necessary
       strategy = 1
      ilin = 1 ! if 1 then uniform random selection, else loguniform
       omagn = 10.0d0 ! order of magnitude for initial range

       ! driver variables. may be left the way they are
       Dim_XC = nopt
       itermax = itmax/np
       F_XC = 0.8d0
       CR_XC = 0.5d0
       refresh = -1
      VTR = 0.d0

       ! define range where to look for minimum
       do i=1,nopt
       XCmin(i) = par(idpar(i))/omagn!lower bound or use XCmin(i) = 0.d0
      XCmax(i) = par(idpar(i))*omagn!upper bound x times > initial guess
         if (XCmin(i).eq.XCmax(i)) then  ! sort
           XCmin(i)=-1.d0
           XCmax(i)=1.d0
           write(*,*) 'DE: p=0, range set to -1..1; not very smart...'
        end if
        if (XCmin(i).gt.XCmax(i)) then
          temp=XCmax(i)
          XCmax(i)=XCmin(i)
          XCmin(i)=temp
        end if
       end do

!      initialize probability and get degrees of freedom
      if (iof.ne.1) write(*,*) 'iof.ne.1 -> but stopping depends on it!'
        if (iof.eq.1) then
          idegf = 0
          do t=1,ntopt
           do k=1,nrspmeas(t)
              do i=1,nrxmeas(t,k)
               idegf=idegf+1
             end do
           end do
         end do
          idegf = idegf - nopt
        else
         prob = 2.d0 ! value greater than 1, so no interference
        end if

!!c     CALL RANDOM_SEED()
           seed = 310952
           call rnfstr(seed)
!!-----Initialize a population --------------------------------------!!
       do i=1,NP
          do j=1,nopt
             pop_XC(i,j)=0.0d0
          enddo
       enddo
       do i=1,NP
!!c        call random_number(rand_C1)
          call rnfarr(rand_C1,nopt)
        do j=1,nopt
           if (ilin.eq.1) then
             ! __linear random initialization
              pop_XC(i,j)=XCmin(j)+rand_C1(j)*(XCmax(j)-XCmin(j))
           else
             ! __log distributed starting points
      pop_XC(i,j)=XCmin(j)+10.0d0**(rand_C1(j)*log10(XCmax(j)-XCmin(j)))
            endif
        end do
       end do

!!--------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------!!
       do i=1,np
        val(i)=0.0d0
       end do
       do j=1,nopt
        popvec(j) = pop_XC(1,j)
       end do
      nfeval=0

       val(1) = funk(popvec)
       bestval=val(1)

       ibest=1

       nfeval=nfeval+1
       do i=2,NP
         do j=1,nopt
          popvec(j) = pop_XC(i,j)
        end do
         val(i) = funk(popvec)
         nfeval=nfeval+1
         if (val(i) < bestval) then
           ibest=i
          bestval=val(i)
         end if
       end do

      do j=1,nopt
         bestmemit_XC(j)=pop_XC(ibest,j)
         bestmem_XC(j)=bestmemit_XC(j)
      end do
!!--------------------------------------------------------------------!!
       do i=1,np
        do j=1,nopt
          bm_XC(i,j)=0.0d0
        end do
      end do
      rot(1) = 0
      do i=2,np
        rot(i) = rot(i-1)+1
      end do

       iter=1
!!------Perform evolutionary computation------------------------------!!
       do while (iter .le. itermax)
        do i=1,np
          do j=1,nopt
             popold_XC(i,j)=pop_XC(i,j)
          end do
        end do

!!------Mutation operation--------------------------------------------!!
        do jj=1,2
         if (jj.eq.1) num = 4  ! ind=randperm(4)
         if (jj.eq.2) num = np ! a1=randperm(NP)
!!c        call random_number(rand2)
            call rnfarr(rand2,np)
          do i=1,num
            number=1
            do j=1,num
              if (rand2(i) .gt. rand2(j)) then
             number=number+1
              end if
            end do
            do k=1,i-1
             if (rand2(i) .eq. rand2(k)) then
               number=number+1
             end if
           end do
           if (jj.eq.1) ind(i)=number
           if (jj.eq.2) a1(i)=number
          end do
         end do

         do i=1,np
           rt(i)=mod(rot(i)+ind(1),NP)
           a2(i)=a1(rt(i)+1)
        end do
         do i=1,np
           rt(i)=mod(rot(i)+ind(2),NP)
           a3(i)=a2(rt(i)+1)
        end do
         do i=1,np
           rt(i)=mod(rot(i)+ind(3),NP)
           a4(i)=a3(rt(i)+1)
        end do
         do i=1,np
           rt(i)=mod(rot(i)+ind(4),NP)
           a5(i)=a4(rt(i)+1)
        end do
         ! bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)
         do i=1,np
          do j=1,nopt
            bm_XC(i,j) = bestmemit_XC(j)
          end do
        end do

! ---- select a mutation strategy--------------------------------------!
         select case (strategy)
           case (1)
          do i=1,np
             do j=i,nopt
               ui_XC(i,j)=bm_XC(i,j)+F_XC*(popold_XC(a1(i),j)
     +                                     -popold_XC(a2(i),j))
            end do
          end do
           case (2)
          do i=1,np
             do j=i,nopt
               ui_XC(i,j)=popold_XC(a3(i),j)+F_XC*(popold_XC(a1(i),j)
     +                                     -popold_XC(a2(i),j))
            end do
          end do
           case (3)
          do i=1,np
             do j=i,nopt
               ui_XC(i,j)=popold_XC(i,j)+F_XC*(bm_XC(i,j)-popold_XC(i,j)
     +                           +popold_XC(a1(i),j)-popold_XC(a2(i),j))
            end do
          end do
           case (4)
           do i=1,np
             do j=i,nopt
               ui_XC(i,j)=bm_XC(i,j)+F_XC*(popold_XC(a1(i),j)
     +        -popold_XC(a2(i),j)+popold_XC(a3(i),j)-popold_XC(a4(i),j))
            end do
          end do
           case (5)
          do i=1,np
             do j=i,nopt
               ui_XC(i,j)=popold_XC(a5(i),j)+F_XC*(popold_XC(a1(i),j)
     +        -popold_XC(a2(i),j)+popold_XC(a3(i),j)-popold_XC(a4(i),j))
            end do
          end do
         end select

!!--------------------------------------------------------------------!!
!!------Crossover operation-------------------------------------------!!
!!c        call random_number(rand_XC)
         do i=1,np
            call rnfarr(rhelp,nopt)
            do jj=1,nopt
               rand_XC(i,jj)=rhelp(jj)
            enddo
         end do
         do i=1,np
          do j=1,nopt
             mui_XC(i,j)=0.0d0
             mpo_XC=0.0d0
             if (rand_XC(i,j).lt.CR_XC) then
               mui_XC(i,j)=1.0d0
               mpo_XC=0.0d0
            else
               mui_XC(i,j)=0.0d0
               mpo_XC=1.0d0
             end if
            ui_XC(i,j)=popold_XC(i,j)*mpo_XC+ui_XC(i,j)*mui_XC(i,j)
          end do
        end do

!!--------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------!!
         do i=1,NP
!!------Confine each of feasible individuals in the lower-upper bound-!!
           do j=1,nopt
            if (ui_XC(i,j).lt.XCmin(j)) ui_XC(i,j) = XCmin(j)
            if (ui_XC(i,j).gt.XCmax(j)) ui_XC(i,j) = XCmax(j)
            popvec(j) = ui_XC(i,j)
          end do

           tempval = funk(popvec)
          nfeval=nfeval+1
           if (tempval .lt. val(i)) then
            do j=1,nopt
               pop_XC(i,j)=ui_XC(i,j)
            end do
             val(i)=tempval
             if (tempval .lt. bestval) then
               bestval=tempval

              do j=1,nopt
                bestmem_XC(j)=ui_XC(i,j)
              end do
             end if
           end if
         end do

        if( (refresh .gt. 0) .and. (mod(iter,refresh).eq.0)) then
             write(3,FMT=203) iter
          write(*, FMT=203) iter
           do i=1,Dim_XC
            write(3, FMT=202) i, bestmem_XC(i)
            write(*,FMT=202) i,bestmem_XC(i)
           end do
           write(3, FMT=201) bestval
           write(*, FMT=201) bestval
         end if
         iter=iter+1
         if ( bestval .le. VTR .and. refresh .gt. 0) then
           write(3, FMT=*) 'The best fitness is smaller than VTR'
           write(*, FMT=*) 'The best fitness is smaller than VTR'
           exit
         endif

        if (iof.eq.1) prob = gammp(dfloat(idegf)/2.d0,bestval/2.d0)

         if(prob.lt.problim) then

          write(3,*) 'prob<problim: ', prob, problim
          write(*,*) 'prob<problim: ', prob, problim
          iter = itermax+1
        end if
       end do
!!------end the evolutionary computation------------------------------!!
201   format(2x, 'bestval=', E14.7, /)
202   format(5x, 'bestmem_XC(', I3, ')=', E12.5)
203   format(2x, 'No. of iteration =', I8)

! output________________

       write(3,*) 'np&nr.function calls: ', NP, nfeval
      write(3,*) 'F and CR: ', F_XC, CR_XC
       write(3,*) 'best objective function: ',bestval
      write(3,*) 'best par: ', (bestmem_XC(i),i=1,Dim_XC)
       write(*,*) 'best OF&par', bestval,(bestmem_XC(i),i=1,Dim_XC)

       return
       end
