c    $Id: optimlm.f 16 2007-10-18 12:36:47Z centler $
      subroutine optimlm()

c***********************************************************************
c      OPTIMLM: optimization using Levenberg-Marquardt                 *
c      brief theory:                                                   *
c      chi^2 (=of) = sum_observations((cmeas-ccalc)/sigma)^2           *
c      minimize of:  d(of)/dpar(j) = 0                                 *
c      iterate, newton type (cut off series expansion)                 *
c      0 = d(of)/d[par(j)+delp]                                        *
c        = d(of)/dpar(j) + sum_parameters[d2(of)/(d(par(j)*d(par(k))]  *
c      d(of)/dpar(j)                                                   *
c      =-2*sum_observations{((cmeas-ccalc)/sigma)^2 * d(ccalc)/dpar(j)}*
c      d2(of)/(dpar(j)*dpar(k))                                        *
c      = 2*sum_obs{ 1/sigma^2 * [d(ccalc)/dpar(j)*d(ccalc)/dpar(k)     *
c                      - (cmeas-ccalc) * d2(ccalc)/(dpar(j)*dpar(k))]} *
c      now introduce matrix alpha and vector beta                      *
c      alpha is 1/2*d2(of)/(dpar(j)*dpar(k)), without 2nd derivative ! *
c      beta is  -1/2*d(of)/dpar(j)                                     *
c      -> solve alpha*delta = beta, where delta is the update of par   *
c      in Levenberg-Marquardt, the diagonal elements of alpha are      *
c      augmented by a factor (1+lambda), and lambda (>0) is adapted    *
c      in the algorithm. Thus, the method switches between deepest     *
c      descent (large lambda) and newton type (small lambda)           *
c                                                                      *
c      algorithm:                                                      *
c      1.choose initial guess for parameters ("a") (from maple)        *
c      2.compute chi^2(a)                                              *
c      3.start with modest lambda (0.001)                              *
c      4.solve alpha*delta = beta, determine chi^2(a+delta)            *
c      5.if chi^2(a+delta)>chi^2(a) increase lambda (*10) and go to 4  *
c      6.if chi^2(a+delta)<chi^2(a) decrease lambda (/10), a=a+delta   *
c      7.check stopping criteria: stop or go back to 4                 *
c      - note that computing chi^2(a) involves a perturbation of each  *
c        individual parameter being optimized...not so cheap           *
c      - stopping: either related to very small change in parameters   *
c        or a low value of chi^2 (for the given degrees of freedom)    *
c                                                                      *
c      subroutines used                                                *
c      olmtime: drives forward calculations                            *
c      olmcalcco: calculates alpha&beta based on concentration profiles*
c      olmstop: contains stopping criteria                             *
c      olmcov: covariance matrix -  uncertainty of optimized parameters*
c                                                                      *
c      CM, Jan 2002                                                    *
c***********************************************************************

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

      integer t,i,j,k,indmeas,jj,iter,idegf,indxx(nopt)
      real*8 bestpar(nopt), optpar(nopt),delp(nopt)
      real*8 spint(maxxmeas), yid(nopt,nopt), beta(nopt)
      real*8 alpha(nopt,nopt),alphasave(nopt,nopt),betasave(nopt)
      real*8 yint(nopt+1,maxxmeas*maxspmeas*ntopt), yof(nopt+1)
      real*8 yofopt, of,lambda
      real*8 sig(maxxmeas*maxspmeas*ntopt*nopt)

       external gammp

c***********************************************************************
c      INITIALIZE & BASE RUN                                           *
c      par: current parameter set being tested                         *
c      bestpar: currently best parameters                              *
c      optpar: best parameters ever encountered                        *
c      sig: standard deviation (sigmeas or relsig*spmeas)              *
c      indmeas: total number of measurements                           *
c      ideg: degrees of freedom for chi^2 distribution                 *
c      iter: iteration counter                                         *
c      spmeas: measured concentrations                                 *
c      yof(nopt+1): objective function for the parameter vertices      *
c      yint: calculated conc., interpolated to measured locations      *
c
c      stdmin: minimum standard deviation                              *
c      relsig&stdmin MUST have same values as in objf.f ..common block?*
c      common_meas.inc, as well as sigmeas                             *
c        (read from file in getdat/storedat.f)                         *
c      NOTE that in objf.f, iof MUST be set to 1 when using optimlm.f  *
c***********************************************************************
        ! determine standard deviation (once, in the beginning)
        indmeas = 0
        do t=1,ntopt
         do k=1,nrspmeas(t)
            do i=1,nrxmeas(t,k)
             indmeas=indmeas+1
             sig(indmeas) = relsig*spmeas(t,k,i)
              if (sigmeas(t,k,i).gt.0.) sig(indmeas) = sigmeas(t,k,i)
             if (sig(indmeas).lt.stdmin) sig(indmeas) = stdmin
           end do
         end do
       end do
        idegf = indmeas - nopt

       ! evaluate the original profile & adjust the settings
        iter = 1
        call olmtime(yof(1),1,yint)
        lambda = 0.001d0
       ofold = yof(1)
        yofopt = yof(1)
       do i=1,nopt
         bestpar(i) = par(idpar(i))
         optpar(i)  = bestpar(i)
       end do

       write(*,*) 'first ', ofold,(par(idpar(i)),i=1,nopt)
       write(3,*) ofold,(par(idpar(i)),i=1,nopt)

c***********************************************************************
c      evaluate the effect of the parameters on the concentrations     *
c      each parameter to be optimized is perturbed one by one          *
c      and the sensitivity of the concentration towards this change is *
c      evaluated (1st derivative, done in olmcalcco                    *
c      idpar: index of parameter that gets optimized                   *
c      delp(nopt) is the change in the parameter value (pertubation)   *
c***********************************************************************

20     do jj=2, nopt+1
          delp(jj-1) = perturb*par(idpar(jj-1))
          par(idpar(jj-1)) = par(idpar(jj-1)) + delp(jj-1)
22        call olmtime(yof(jj),jj,yint)

      ! ___addition to avoid singular matrix
         if (yof(jj).eq.yof(1)) then
            par(idpar(jj-1)) = par(idpar(jj-1)) - delp(jj-1)
           delp(jj-1) = 10.0d0 * delp(jj-1)
            par(idpar(jj-1)) = par(idpar(jj-1)) + delp(jj-1)
           write(*,*) 'increase dp',jj,iter
            write(*,*) yof(1), yof(jj), par(idpar(jj-1)),delp(jj-1)
           go to 22
         end if
      ! ___end addition

         iter = iter+1
         write(*,*) 'vertex', yof(jj),(par(idpar(i)),i=1,nopt)
         write(3,*)  yof(jj),(par(idpar(i)),i=1,nopt)
         ! check if a new overall best profile was found
         if (yof(jj).lt.yofopt) then
           yofopt = yof(jj)
           do k=1,nopt
             optpar(k) = par(idpar(k))
              bestpar(k) = par(idpar(k))  ! new
            end do
         end if
        end do ! end parameter perturbation loop

c***********************************************************************
c      calculate alpha&beta                                            *
c***********************************************************************

        call olmcalcco(alpha,beta,yint,delp,sig)

c***********************************************************************
c      update parameters by solving alpha*delta_par = beta             *
c      delta_par = inv(alpha)*beta                                     *
c      ludcmp: on output, alpha is its LU decomposition ("LUd")        *
c      lubksb: beta becomes the suggested change for the parameters,   *
c      thus save both alpha&beta and restored them afterwards again    *
c***********************************************************************

        ! conserve (alpha is conserved without the lambda term!)
       do i=1,nopt
          betasave(i) = beta(i)
         do k=1,nopt
           alphasave(i,k) = alpha(i,k)
         end do
       end do

! evaluate change in parameters (after adding lambda term to alpha)
10      do i=1,nopt
          alpha(i,i) = alpha(i,i) * (1.d0+lambda)
        end do
      ! ____start addition avoiding singular matrix
         do 12 i = 1,nopt
           aamax = 0.0d0
           do 11 j = 1,nopt
             if (dabs(alpha(i,j)) .gt. aamax) aamax=dabs(alpha(i,j))
   11     continue
           if (aamax .eq. 0.0d0) then
              write(*,*) 'Singular matrix despite large increase in p'
             do k=1,nopt
               write(*,*) (alpha(k,j),j=1,nopt)
              end do
              stop
              ! rather than stop go up and modify the parameter profile
c             pause 'if you want to restart hit enter'
c            perturb = perturb*100.0d0
c             iter = 2
c            go to 20
           endif
   12   continue
      ! ____end addition

c     xxx lud
        call ludcmp(alpha,nopt,nopt,indxx,det)
        call lubksb(alpha,nopt,nopt,indxx,beta)

      ! update parameters and restore alpha (without lambda) and beta
        do i=1,nopt
          par(idpar(i)) = bestpar(i) + beta(i) ! beta is change in par
         if (ngtzero.eq.1) then
          if (par(idpar(i)).lt.0.d0) then ! avoid negative parameters
            if (bestpar(i).lt.0.d0) par(idpar(i)) = -bestpar(i)*0.1d0
            if (bestpar(i).gt.0.d0) par(idpar(i)) = bestpar(i)*0.1d0
            if (bestpar(i).eq.0.d0) par(idpar(i)) = 0.1d0
          end if
         endif
         beta(i) = betasave(i)                ! restore beta
          do k=1,nopt
           alpha(i,k) = alphasave(i,k)! restore alpha (w.o. lambda term)
         end do
       end do

c***********************************************************************
c      evaluate obj.function at this new trial parameter set           *
c***********************************************************************

        call olmtime(yof(1),1,yint)
       iter = iter+1
       write(*,*) 'trial ', yof(1),(par(idpar(i)),i=1,nopt)
       write(3,*) yof(1),(par(idpar(i)),i=1,nopt)

c***********************************************************************
c      compare the new unperturbed profile to the previous             *
c      unperturbed parameter profile                                   *
c      if new one is better, check if stopping criteria are met        *
c      -> stop or proceed with new one as new starting point           *
c      if old one is better, go back and solve for new delta_parameter *
c       using a larger value for lambda                                *
c      note that if stopping criterium is met, iter = itmax in olmstop *
c***********************************************************************

        if (iter.lt.itmax) then
         if (yof(1).gt.ofold) then
            lambda = lambda*10.d0
            go to 10 !solve for new dpar with larger lambda,same profile
          end if
          if (yof(1).le.ofold) then
           lambda = lambda*0.1d0
           do k=1,nopt
             bestpar(k) = par(idpar(k)) ! updata best parameter set
           end do
           call olmstop(ofold,yof,yofopt,optpar,idegf,iter)
         end if
       else
         write(*,*) 'exceeding itmax'
         call olmstop(ofold,yof,yofopt,optpar,idegf,iter)
        end if
        if (iter.lt.itmax) go to 20 ! new node, full perturbation

c***********************************************************************
c      calculate covariance matrix - measure of parameter uncertainty  *
c      at this point, alpha is not including lambda (i.e. lambda"="0)  *
c***********************************************************************

        call olmcov(alpha)

      return
      end

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       subroutine olmtime(of,jj,yint)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

       real*8 of,tstart,tend
      real*8 spint(maxxmeas), yint(nopt+1,maxxmeas*maxspmeas*ntopt)
      integer t,j,k,jj,imeas

        call transferback() ! assign to their real maple names

c***********************************************************************
c      TIME LOOP                                                       *
c     tstart  = starting time of simulation. either 0 or last time     *
c               measurements have already been compared to calculated  *
c               values                                                 *
c     tend    = end time of next simulation. time of next measurement  *
c     ntopt   = number of different times measurements are available   *
c***********************************************************************
        imeas = 0
         of = 0.d0
        do t=1,ntopt
          if (t.eq.1) then
            tstart = 0.d0
          else
            tstart = timemeas(t-1)
          end if
          tend = timemeas(t)

c***********************************************************************
ccalculated concentration profiles at the time of measurements   
c***********************************************************************
          call diagenesis(tstart,tend)

c***********************************************************************
c     GET THE OBJECTIVE FUNCTION AND                                   *
c      STORE THE INTERPOLATED CALCULATED CONC.                         *
c     of  = value of the objective function                            *
c     spint = calculated conc., interpolated to measurements           *
c     yint(nopt+1,indmeas) = sequence of all spint's.                  *
c     sig(indmeas) is the standard deviation                           *
c     yof stores the objective function values                         *
c     indmeas = current next slot to be filled in yint                 *
c***********************************************************************
          do k=1,nrspmeas(t)
             j=idspmeas(t,k) ! next measured component
             call objf(of,spint,j,k,t)
             do i=1,nrxmeas(t,k)
              imeas=imeas+1 ! index for all measurements
              yint(jj,imeas) = spint(i)
            end do
          end do ! end measured species loop
        end do ! end time loop

       return
      end

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       subroutine olmcalcco(alpha,beta,yint,delp,sig)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'

        real*8 alpha(nopt,nopt),beta(nopt),delp(nopt)
       real*8 yint(nopt+1,maxxmeas*maxspmeas*ntopt)
       real*8 sig(maxxmeas*maxspmeas*ntopt)
        integer icount,j,t,k,i

c***********************************************************************
c      calculate the derivatives dydp, i.e. the effect of the parameter*
c      on the measurements, to create the matrix alpha                 *
c      dydp(total nr. measurements,nopt): change of conc. w. parameters*
c      calculate the derivatives dyofdp to make the vector beta, i.e.  *
c      the effect  of the parameters on the objective function         *
c      dyofdp(nopt,nopt): change of the obj. function w. parameters    *
c***********************************************************************
        do j=1,nopt
         icount = 0
         beta(j) = 0.d0
         do t=1,ntopt
           do k=1,nrspmeas(t)
              do i=1,nrxmeas(t,k)
               icount = icount+1
               beta(j)= beta(j) + (spmeas(t,k,i)-yint(1,icount))
     +             /sig(icount)/sig(icount)
     +             * (yint(j+1,icount)-yint(1,icount))/delp(j)
             end do ! end measured species loop
           end do   ! end component loop
         end do ! end time loop

         ! calculate alpha
         do k=1,nopt
            alpha(j,k) = 0.d0
           do i=1,icount
             alpha(j,k) = alpha(j,k) + (yint(j+1,i)-yint(1,i))/delp(j)*
     +        (yint(k+1,i)-yint(1,i))/delp(k) / sig(i)/sig(i)
            end do
          end do
        end do

       return
      end

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       subroutine olmstop(ofold,yof,yofopt,optpar,idegf,iter)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

       integer i,idegf,iter
      real*8 yof(nopt),yofopt,ofold,yid(nopt,nopt),optpar(nopt)
      real*8 rtol,gammp,prob

      external gammp

c***********************************************************************
c      stopping condition                                              *
c      1. is the difference between old and new objective functions    *
c          relative to their average magnitude small enough?           *
c      2. is chi^2 (=obj.function) smaller than the chi^2 for a given  *
c         confidence level and the relevant degrees of freedom?        *
c         (is the calculated chi^2 indicating that the concentrations  *
c          obtained with the current parameter differs from the data   *
c          only due to random noise - assuming normal distribution     *
c      if 1 or 2 is fulfilled, check if a better profile was found once*
c      idegf: degree of freedoms = nr. measured datapoints - nopt      *
c      prob: probability that the observed chi^2 will exceed the       *
c            value chi^2 (i.e. obj.f.) by chance even for a correct    *
c            model than current chi^2 -> accept as final solution      *
c            if prob < problim                                         *
c***********************************************************************
       rtol=2.d0*dabs(ofold-yof(1))/(dabs(ofold)+dabs(yof(1)))
      ofold = yof(1)
      prob = gammp(dfloat(idegf)/2.d0,ofold/2.d0)

!      if((rtol.lt.optftol).or.(prob.lt.problim).or.(iter.ge.itmax))then
       if((rtol.lt.optftol).or.(prob.lt.problim)) then
         if (ofold.gt.yofopt) then
          write(*,*) 'stopping, but better profile elsewhere!'
          write(*,*) 'ofold, yofopt: ', ofold, yofopt
          write(3,*) 'stopping, but better profile elsewhere!'
          write(3,*) 'ofold, yofopt: ', ofold, yofopt
          prob = gammp(dfloat(idegf)/2.d0,yofopt/2.d0)
          do i=1, nopt
            par(idpar(i)) = optpar(i)
          end do
        end if

c***********************************************************************
c      put the optimized parameters back into the complete list        *
c      assign them to their real names (transferback)                  *
c      and give final output                                           *
c      this includes computing the covariance matrix.                  *
c***********************************************************************
         call transferback()
        write(*,*) 'iterations: ', iter
         iter = itmax ! signals stop

        write (*,*) 'BEST PARAMETER VALUES, prob:',prob
        write (3,*) 'BEST PARAMETER VALUES, prob:', prob
        do i=1,nopt
          write (*,199) i, par(idpar(i))
          write (3,199) i, par(idpar(i))
        end do
        write(3,*) ' '
199     format(i5,2x,100(e14.7, 2x))

       end if ! end rtol

       return
      end

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       subroutine olmcov(alpha)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

       integer k,kk,indxx(nopt)
      real*8 yid(nopt,nopt),alpha(nopt,nopt)

        do k=1,nopt ! set up identity matrix
          do kk=1,nopt
            yid(k,kk) = 0.d0
          end do
          yid(k,k) = 1.d0
        end do

c     LUD choice
      if(ilud.eq.1)then
         call ludcmp(alpha,nopt,nopt,indxx,det) ! LUdecompose alpha
        do k=1,nopt
          call lubksb(alpha,nopt,nopt,indxx,yid(1,k))
         end do
      else
           call DGESV(nopt,nopt,alpha,nopt,indxx,yid,nopt,info)
           if(info.ne.0)then
                if(info.lt.0)write(*,*)'illegal value in LUD', info
                if(info.gt.0)write(*,*)'singular matrix in LUD'
           endif
      endif

        write(3,*) 'parameter covariances: sigma2 in diagonals'
        do k=1,nopt
          write(3,199) (alpha(k,kk), kk=1,nopt)
        end do
199     format(100(e14.7, 2x))

       return
      end

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
