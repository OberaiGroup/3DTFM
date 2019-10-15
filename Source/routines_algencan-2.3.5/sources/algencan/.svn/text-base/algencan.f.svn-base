C     ******************************************************************
C     ******************************************************************

      subroutine algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,
     & equatn,linear,coded,checkder,fu,cnormu,snorm,nlpsupn,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,iprint,m,n,ncomp
      double precision cnormu,epsfeas,epsopt,fu,nlpsupn,snorm

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"
#include "machconst.com"
#include "algconst.par"
#include "counters.com"
#include "outtyp.com"
#include "algparam.com"
#include "scaling.com"
#include "slacks.com"
#include "fixvar.com"

C     PARAMETERS
      integer n0,n1,n2
      parameter ( n0 =   500 )
      parameter ( n1 = 10000 )
      parameter ( n2 = 20000 )

C     LOCAL SCALARS
      logical innfail,lss,scl
      integer alinfo,geninfo,i,iter,j,maxit,nwcalls,nwtotit,outiter,
     +        solinfo,totiter,msqcalls,msqtotit
      double precision cnorm,cnormb,cnormub,f,fb,fub,ncsupn,nlpsupnb,
     +        rsupn
      real time

C     LOCAL ARRAYS
      character * 4 lsssub,sclsub
      double precision nl(nmax),c(mmax),cu(mmax),rho(mmax)
      real dum(2)

C     DATA STATEMENTS
      data dum/0.0,0.0/

C     EXTERNAL FUNCTIONS AND SUBROUTINES
      external fparam,checkd,auglag,gencan,lss,scl

C     ==================================================================
C     Start timing
C     ==================================================================

      !time = dtime(dum)
      call dtime(dum,time) 
C     ==================================================================
C     Initialization
C     ==================================================================

C     Set machine-dependent constants

      bignum    = 1.0d+99
      macheps   = 1.0d-16
      macheps12 = sqrt( macheps )
      macheps13 = macheps ** ( 1.0d0 / 3.0d0 )
      macheps23 = macheps ** ( 2.0d0 / 3.0d0 )

C     Set global counters

      fcnt    = 0
      efcnt   = 0
      efccnt  = 0
      egcnt   = 0
      egjccnt = 0
      ehcnt   = 0
      ehlcnt  = 0
      ehlpcnt = 0

      do j = 1,m
          eccnt(j)  = 0
          ejccnt(j) = 0
          ehccnt(j) = 0
      end do

C     ==================================================================
C     Set default values for algoritmic parameters
C     ==================================================================

C     Set user-provided subroutines indicators

      fcoded    = coded(1)
      gcoded    = coded(2)
      hcoded    = coded(3)
      ccoded    = coded(4)
      jaccoded  = coded(5)
      hccoded   = coded(6)
      fccoded   = coded(7)
      gjaccoded = coded(8)
      hlcoded   = coded(9)
      hlpcoded  = coded(10)

      innercall = .false.
      useustp   = .false.

C     Set indicator of whether the true Hessian of the Lagrangian can be
C     computed or not

      truehl = .false.
      if ( hlcoded .or. ( hcoded .and. ( hccoded .or. m .eq. 0 ) )) then
          truehl = .true.
      end if

C     Hessian-vector product strategy: HAPPRO, INCQUO or TRUEHP. TRUEHP
C     is the default option. If the proper subroutines were not coded by
C     the user, then HAPPRO is used instead. HAPPRO reduces to INCQUO in
C     the unconstrained and bound-constrained cases and switches between
C     INCQUO and the product by an approximation of the Hessian in the
C     constrained case.

      if ( truehl .or. hlpcoded ) then
          hptype = 'TRUEHP'
      else
          hptype = 'HAPPRO'
      end if

C     Inner solver method: TR (trust regions, i.e. BETRA), NW (Newtonian
C     system with the use of a direct solver, i.e. MA57AD), TN + TRUEHP
C     (Truncated Newton and true Hessian-vector product) or TN + HAPPRO
C     (Truncated Newton and product with Hessian approximation). If the
C     proper subroutines were not provided by the user then TN + HAPPRO
C     is used.

      if ( truehl .and. lss(lsssub) .and. n .le. n0 ) then
          innslvr = 'TR'
      else if ( truehl .and. lss(lsssub) .and. n .le. n1 ) then
          innslvr = 'NW'
      else
          innslvr = 'TN'
          if ( n .gt. n2 ) then
              hptype  = 'HAPPRO'
          end if
      end if

C     Scaling of linear systems

      if ( lsssub .eq. 'MA57' .or. scl(sclsub) ) then
          sclsys = .true.
      else
          sclsys = .false.
      end if

C     Ignore objective function (to only find a feasible point by
C     minimizing 1/2 of the squared infeasibility)
      ignoref = .false.

C     Acceleration step
      if ( truehl .and. lss(lsssub) ) then
          skipacc = .false.
      else
          skipacc = .true.
      end if

C     Slacks for inequality constraints
      slacks = .false.

C     Remove fixed variables (with identical lower and upper bounds)
      rmfixv = .true.

C     Scale objective function and constraints
      if ( m .gt. 0 .and. .not. ignoref ) then
          scale = .true.
      else
          scale = .false.
      end if

C     ==================================================================
C     Main output control (silent-mode?)
C     ==================================================================

      iprintctl(1) = .true. ! Banner
      iprintctl(2) = .true. ! Parameters and problem processing
      iprintctl(3) = .true. ! Warnings and errors messages
      iprintctl(4) = .true. ! Screen-mirror file algencan.out
      iprintctl(5) = .true. ! Solution file solution.txt
      iprintctl(6) = .true. ! Statistic files with table lines

      open(10,err=100,file='.silent',status='old')
      close(10)

      do i = 1,6
          iprintctl(i) = .false.
      end do

      iprint = 0

 100  continue

      if ( iprintctl(4) ) then
          open(unit=10,file='algencan.out',status='replace')
      else
          open(unit=10,                    status='scratch')
      end if

C     ==================================================================
C     Set solver arguments using the specification file
C     ==================================================================

      call fparam(epsfeas,epsopt,iprint,ncomp)

C     Outer and inner iterations output detail

      iprintout = iprint / 10
      iprintinn = mod( iprint, 10 )

C     Error tracker

      inform = 0

C     ==================================================================
C     Initialize problem data structures
C     ==================================================================

      call sinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

      nprint = min( n, ncomp )
      mprint = min( m, ncomp )

C     ==================================================================
C     Call the solver
C     ==================================================================

C     ALGENCAN for PNL problems

      if ( .not. ignoref .and. m .gt. 0 ) then
          call auglag(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,f,c,
     +    cnorm,snorm,nl,nlpsupn,fu,cu,cnormu,fub,cnormub,fb,cnormb,
     +    nlpsupnb,ncsupn,rsupn,outiter,totiter,nwcalls,nwtotit,
     +    msqcalls,msqtotit,innfail,alinfo,inform)

          solinfo = alinfo

C     GENCAN for box-constrained problems and feasibility problems

      else
          maxit = 999999999

C         Used in feasibility problems (ignoref=true). With lambda=0 and
C         rho=1, to minimize 1/2 of the squared infeasibility coincides
C         with minimizing the augmented Lagrangian.
          do j = 1,m
              lambda(j) = 0.0d0
              rho(j)    = 1.0d0
          end do

          call gencan(n,x,l,u,m,lambda,equatn,linear,rho,epsfeas,epsopt,
     +    maxit,iter,f,nl,nlpsupn,cnorm,cnormu,geninfo,inform)

          solinfo  = geninfo

          ncsupn   = 0.0d0
          rsupn    = 0.0d0

          outiter  = 0
          totiter  = iter
          nwcalls  = 0
          nwtotit  = 0
          msqcalls = 0
          msqtotit = 0
          innfail  = .false.

          if ( ignoref ) then
              f = 0.0d0
          end if
          fu = f

          fb       = f
          fub      = fu
          cnormb   = cnorm
          cnormub  = cnormu
          nlpsupnb = nlpsupn
      end if

      if ( inform .lt. 0 ) return

C     ==================================================================
C     End problem data structures
C     ==================================================================

      call sendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     Close output file
C     ==================================================================

      close(10)

C     ==================================================================
C     Stop timing
C     ==================================================================

      !time = dtime(dum)
      call dtime(dum,time) 
      time = dum(1)

      if ( iprintout .ge. 1 ) then
          write(*,9000) time
      end if

C     ==================================================================
C     Write statistics
C     ==================================================================

      if ( iprintctl(6) ) then
          open(20,file='algencan-tabline.out')
          write(20,9040) fu,cnormu,f,cnorm,nlpsupn,fub,cnormub,fb,
     +                   cnormb,nlpsupnb,ncsupn,rsupn,inform,solinfo,
     +                   innfail,n,m,outiter,totiter,fcnt,nwcalls,
     +                   nwtotit,msqcalls,msqtotit,time
          close(20)
      end if

C     ==================================================================
C     NON-EXECUTABLE STATEMENTS
C     ==================================================================

 9000 format(/,1X,'Total CPU time in seconds: ',F8.2)
 9040 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D7.1,1X,1P,D7.1,1X,I3,1X,I1,1X,L1,1X,I6,1X,I6,1X,I3,
     +       1X,I7,1X,I7,1X,I2,1X,I7,1X,I7,1X,I7,0P,F8.2)

      end
