c******************************************************
      subroutine driverForward()
c     forward solve with the finite element code
      implicit none
      integer imeas
      double precision g
c------------------------------------------------------
      imeas=1! and imeas=nmeas=1 for the forward solve
      g=-1234567.0d0! dummy value; used for subroutine result
      

c     update idnum for the given measurement case
      call init_idnum(imeas)

c     Create the matrices (csr sparse format)
      call inicsr
      call timer(4)

      call forwardSolveEnforceBC(imeas)

      call result(0,g)
      end


c******************************************************
      subroutine driverInverse(iopt,niter,n,bfgsM,noutput,prefixOut)
c     gradient based optimization algorithm
c     this does not view any of the variables in MAINMEM or IOUNIT
      implicit none
      integer iopt,niter,n,bfgsM,noutput
      character*80 prefixOut
      double precision x(n),l(n),u(n)
      integer nF
c------------------------------------------------------

c     fill the state variable (x), the lower and upper bounds for each variable (l and u)
      call init_optvars2(x,l,u)

c     initialize the state vector in a better way (not implemented - 2010/04/20)

c     Write an additional text file for checking convergence
      nF=1
      do while (prefixOut(nF:nF).ne." ")
        nF=nF+1
      enddo
      prefixOut(nF:(nF+3))='.cvg'
      open(16,file=prefixOut(1:(nF+3)),status='unknown')
      prefixOut(nF:(nF+3))='    '
      write(16,*)'#This file helps checking about the convergence of',
     $   ' the optimization algorithm'
      write(16,*)'#    by printing infos ',
     $   'related to the values of the objective function and'
      write(16,*)'#    its components.'
      write(16,*)'#- Every line describes an ',
     $   'iteration of the optimization algorithm,'
      write(16,*)'#- on every line there are x*2+2+4 real numbers',
     $   ' where x equals the number'
      write(16,*)'#    of measured fields used,'
      write(16,*)'#- the first 2 numbers on a line are the dataMatch',
     $   ' and the forceMatch terms'
      write(16,*)'#    for the first measurement used,'
      write(16,*)'#- if x>1, numbers at position 3 and 4 are the ',
     $   'dataMatch and the forceMatch'
      write(16,*)'#    terms for the second measurement'
      write(16,*)'#- etc...'
      write(16,*)'#- the following numbers are: the regularization ',
     $   'term for each measurement' 
      write(16,*)'#    and the objective function value,'
      write(16,*)'#- followed by the L2 norms of',
     $   ' parameter1, its gradient, parameter2 and its'
      write(16,*)'#    gradient.'
      write(16,*)'#note by JFD: gnuplot is convenient to print a column'
      write(16,*)'#             plot "*.out.cvg" u 3'
      write(16,*)'#             replot "*.out.cvg" u 6'

c     choose the optimization method depending on the iopt value
      if (iopt.eq.1) then
c       Steepest descent optimization (not implemented yet - 2010/04/20)
        write(*,*) 'the chosen optimization is not implemented ...',
     $       ' exiting!'
c        stop
        call driversd(niter,n,noutput,x,l,u,prefixOut)
      elseif (iopt.eq.2) then
c       L-BFGS-B optimization
        call driverbfgs(niter,n,bfgsM,x,l,u,prefixOut)
      elseif (iopt.eq.3) then
c       Newton optimization (not implemented yet 2010/04/20)
c         (very slow as the hessian is computed with finite differences...)
        write(*,*) 'the chosen optimization is not implemented ...',
     $       ' exiting!'
        stop
      elseif (iopt.eq.4) then
c       ASA_CG optimization
        call driverASACG(niter,n,x,l,u,prefixOut)
      elseif (iopt.eq.5) then
c       ASA_CG optimization
        call driverAlgencan(niter,n,x,l,u,prefixOut)
      else
        write(*,*) 'the chosen optimization does not exist yet ...',
     $       ' exiting!'
        stop
      endif

      close(16,status='keep')! close the cvg file
      
 
      return
      end


c******************************************************
c                SIMPLE DRIVER FOR L-BFGS-B (version 2.1)
c
c        L-BFGS-B is a code for solving large nonlinear optimization
c             problems with simple bounds on the variables.
c
c        The code can also be used for unconstrained problems and is
c        as efficient for these problems as the earlier limited memory
c                          code L-BFGS.
c
c        This is the simplest driver in the package. It uses all the
c                    default settings of the code.
c
c     References:
c        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c        memory algorithm for bound constrained optimization'',
c        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c        Subroutines for Large Scale Bound Constrained Optimization''
c        Tech. Report, NAM-11, EECS Department, Northwestern University,
c        1994.

c          (Postscript files of these papers are available via anonymous
c           ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

c                              *  *  *

c        NEOS, November 1994. (Latest revision June 1996.)
c        Optimization Technology Center.
c        Argonne National Laboratory and Northwestern University.
c        Written by
c                           Ciyou Zhu
c        in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

      subroutine driverbfgs(niter,n,m,x,l,u,prefixOut)
      implicit none
      integer          niter,n,m
      double precision x(n),l(n),u(n)
      character*80     prefixOut
c     A description of all these variables is given at the end of the driver
      character*60     task, csave
      logical          lsave(4)
      integer          iprint
      integer          nbd(n), iwa(3*n), isave(44)
      double precision f, factr, pgtol
      double precision g(n), dsave(29), wa(2*m*n+4*n+12*m*m+12*m)
c     Additional variables for this problem
      integer          i
c------------------------------------------------------

c     Output at every iteration
      iprint = 1
      
c     Specify the tolerances in the stopping criteria
      factr=1.0d-15! make factr really small as the code shold stop in lmain.f/gradfun
      pgtol=1.0d-20! make pgtol really small as the code shold stop in lmain.f/gradfun

      f = 1.0d8! initialization of the objective function (random large number)
 
c     Provide nbd which defines the bounds on the variables:
      do i=1,n
         nbd(i)=2! 2 means that there is both lower and upper bounds for each variable
      enddo

c     Start the iteration by initializing task
      task = 'START'

c     Beginning of the loop 

      call timerIteration(0) ! initialization

 111  continue
      
c     Call the L-BFGS-B code
      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave,prefixOut)

      if (task(1:2) .eq. 'FG') then
c        The minimization routine has returned to request the
c        function f and gradient g values at the current x
c        send the iteration number for file&screen printing
         call gradfun(x,f,g,isave(30))
  
c        Go back to the minimization routine
         goto 111
      endif

c      if (task(1:5).eq.'NEW_X'.and.isave(30).lt.niter)  goto 111
      if (task(1:5).eq.'NEW_X')  goto 111
c     The minimization routine has returned with a new iterate,
c         and we have opted to continue the iteration

c     End of the loop
c      print*,"for debug: task(1:5)=",task(1:5)
      ! if the gradient should be printed at the end of the simulation, something should be done here
      

c     If task is neither FG nor NEW_X we terminate execution

      return
      end

c     --------------------------------------------------------------
c     DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave
c
c     END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------


c******************************************************
      subroutine driversd(niter,n,noutput,x,l,u,prefixOut)
c     gradient based optimization algorithm
c     this does not view any of the variables in MAINMEM
      implicit none
      integer niter,n,noutput
      double precision x(n),l(n),u(n)
      character*80     prefixOut

      integer iter, nextOutp, i, cnt,it
      double precision f,f0, f1, f2, f3,fT, lineSearchDistance, xtmp
      double precision lambda0, lambda1, lambda2, lambdaT
      double precision g(n),g0(n), g1(n), g2(n), g3(n)
      double precision x0(n), x1(n), x2(n), x3(n)
      double precision aCoeff, bCoeff, g0InProd
c------------------------------------------------------
      f = 0.0d0
      g(:) = 0.0d0     
      f0 = 0.0d0
      g0(:) = 0.0d0
      f1 = 0.0d0
      g1(:) = 0.0d0
      f2 = 0.0d0
      g2(:) = 0.0d0
      f3 = 0.0d0
      g3(:) = 0.0d0
      nextOutp=0! output at iteration 0 (to check the loading of the material)

      do iter=0,niter
c        exit
        print *, 'itern=', iter
        call gradfun(x, f0, g0, iter)
c       output the results
        if (iter.eq.nextOutp) then
          call result(iter,g0)
          nextOutp = nextOutp+noutput
        endif
        lambda0 = 1
        do i=1,n
          x1(i) = x(i) - lambda0*g0(i)
          if (x1(i)<l(i)) x1(i)=l(i)
          if (x1(i)>u(i)) x1(i)=u(i)
        enddo
        call gradfun(x1, f1, g1, iter)
        g0InProd = 0
        do i=1,n
          g0InProd=g0InProd+g0(i)**2
        enddo
        cnt = 0
        lambdaT=lambda0
        fT=f1
        do  while(fT> (f0 -0.0001*lambdaT*g0InProd))
          if(cnt.eq.0) then
            print *, 'cnt=', cnt
            lambda1 = g0InProd/(2*(f1-f0+g0InProd))
            if(lambda1<0.1*lambda0) lambda1=0.1*lambda0
            if(lambda1>0.5*lambda0) lambda1=0.5*lambda0
            do i=1,n
              x2(i) = x(i) - lambda1*g0(i)
              if (x2(i)<l(i)) x2(i)=l(i)
              if (x2(i)>u(i)) x2(i)=u(i)
            enddo
            call gradfun(x2, f2, g2, iter)
            lambdaT=lambda1
            fT=f2
          else
            print *, 'cnt=', cnt
            aCoeff=1/(lambda1-lambda0)*
     &       ((f2-f0+g0InProd*lambda1)/lambda1**2
     &       -(f1-f0+g0InProd*lambda0)/lambda0**2)
            bCoeff=1/(lambda1-lambda0)*
     &       (-(f2-f0+g0InProd*lambda1)*lambda0/lambda1**2
     &       +(f1-f0+g0InProd*lambda0)*lambda1/lambda0**2)
            lambda2=(-bCoeff+sqrt(bCoeff**2+3*aCoeff*g0InProd))
     &        /(3*aCoeff)
            if(lambda2<0.1*lambda1) lambda2=0.1*lambda1
            if(lambda2>0.5*lambda1) lambda2=0.5*lambda1
            do i=1,n
              x3(i) = x(i) - lambda2*g0(i)
              if (x3(i)<l(i)) x3(i)=l(i)
              if (x3(i)>u(i)) x3(i)=u(i)
            enddo
            call gradfun(x3,f3,g3,iter)
            lambda0 = lambda1
            lambda1 = lambda2
            lambdaT=lambda2
            f1=f2
            f2=f3
            fT=f3
          endif
          cnt = cnt+1
          print *, 'lambdaT=', lambdaT
          print *, 'fT=', fT
          print *, 'rhs=', (f0 -0.0001*lambdaT*g0InProd)
          print *, 'f0=', f0
          print *, 'g0InProd=', g0InProd
        enddo
        do i=1,n
          x(i) = x(i) - lambdaT*g0(i)
          if (x(i)<l(i)) x(i)=l(i)
          if (x(i)>u(i)) x(i)=u(i)
        enddo
      enddo
      
      do it=0,niter
        exit
        print *, 'before',' f=',f
        call gradfun(x,f,g,it)
        print *, 'after', ' f=',f
        lineSearchDistance=1.0d-2
        do i=1,n
          x(i)=x(i)-lineSearchDistance*g(i)
          if (x(i)<l(i)) x(i)=l(i)
          if (x(i)>u(i)) x(i)=u(i)
        enddo
      enddo


      return
      end


c******************************************************
      subroutine driverASACG(niter,n,x,l,u,prefixOut)
c     gradient based optimization algorithm
c     this does not view any of the variables in MAINMEM
      implicit none
      integer niter,n
      double precision x(n),l(n),u(n)
      character*80     prefixOut

c     call a c function
      call driver_asa_cg(n,x,l,u,1000000.0d0)
c     note that the max number of iterations is only approximate as ASA_CG counts iterations in a weird way.

      return
      end


c******************************************************
      subroutine driverAlgencan(niter,n,x,l,u,prefixOut)
c     gradient based optimization algorithm
c     this does not view any of the variables in MAINMEM
c     modified from algencanma.f (the main program for algencan)
      implicit none

C     ARGUMENTS
      integer niter,n
      double precision x(n),l(n),u(n)
      character*80     prefixOut

C     LOCAL VARIABLES
      logical checkder
      integer inform,iprint,m,ncomp
      double precision cnorm,epsfeas,epsopt,f,nlpsupn,snorm

C     LOCAL ARRAYS
      logical coded(10),equatn(0),linear(0)! 0 refers to the value of m
      double precision lambda(0)! 0 refers to teh value of m

C     EXTERNAL SUBROUTINES
c      external algencan!,param

C     SET UP PROBLEM DATA

c      call param(epsfeas,epsopt,iprint,ncomp)! set some parameters
      ! make the tolerances really small to let lmain.f:gradfun define the stopping criterion
      epsfeas  = 1.0d-14! feasibility stopping criterion: threshold for the sup-norm of the vector of infeasibility
      epsopt   = 1.0d-14! optimality (KKT) stopping criterion
      iprint   = 10! first digit is for outer-loop printing and second is for the inner-loop; the higher the digits, the more printing on the screen (default=10 that is 1 and 0)
      ncomp    = 6! number of first entries of x, l or u that will be printed on the screen (default=6)
      m=0! number of additional constraints
      f=0.0d0
      cnorm=0.0d0! ask EG
      snorm=0.0d0! ask EG
      nlpsupn=0.0d0! ask EG
      inform=0! flag for normal or abnormal termination


      call inip(n,x,l,u,m,lambda,equatn,linear,coded,checkder)

c      print*,"epsfeas=",epsfeas,", epsopt=",epsopt
c      print*,"iprint=",iprint,", ncomp=",ncomp
c      print*,"n=",n,", m=",m
c      print*,"x=",x
c      print*,"l=",l
c      print*,"u=",u
c      print*,"lambda=",lambda
c      print*,"equatn=",equatn
c      print*,"linear=",linear
c      print*,"coded=",coded
c      print*,"checkder",checkder
c      print*,"f=",f

      call algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,
     &  equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)

c      call endp(n,x,l,u,m,lambda,equatn,linear)! for post-processing if needed

      end


