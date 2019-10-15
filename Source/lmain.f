c*******************************************************
      program FEACOREADJOINT
c     written by Assad Oberai
c     modified by Sevan Goenezen, Jean-Francois Dord
c     
c     Notes: 
c     (1) "The finite element method," by TJRH
c     (2) "File formats for vtk version 4.2" 
c         (http://public.kitware.com/VTK/pdf/file-formats.pdf)
c     
c     Acknowledgements:
c     (1) Peter M Pinsky, Stanford University
c     (2) Manish Malhotra
c     (3) Carlos Rivas, Boston University
c******************************************************
      USE MAINMEM
      USE IOUNIT
      implicit none
      integer nargin,nF
      character*80 argin1
      logical ex
      double precision g
c------------------------------------------------------
cC     The following line can be uncommented to test openmp 
c      INTEGER NTHREADS, TID! variable for test
cc      INTEGER OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM! variables for test
cC     Fork a team of threads giving them their own copies of variables
c!$OMP PARALLEL PRIVATE(NTHREADS, TID)
cC     Obtain thread number
c      TID = OMP_GET_THREAD_NUM()
c      PRINT *, 'Hello World from thread = ', TID
cC     Only master thread does this
c      IF (TID .EQ. 0) THEN
c        NTHREADS = OMP_GET_NUM_THREADS()
c        PRINT *, 'Number of threads = ', NTHREADS
c      END IF
cC     All threads join master thread and disband
c!$OMP END PARALLEL
c      stop

c     BEGINNING OF NLACE
      call timerMain(0) ! initialize timing
  
c     GET THE NAME OF AN INPUT FILE GIVEN ON THE COMMAND LINE
      argin1(1:40)='                                        '
      argin1(41:80)='                                        '
c      nargin=iargc() ! number of input arguments - f77 syntax
      nargin=command_argument_count() ! number of input arguments - f90 syntax
      if (nargin.eq.0) then
        argin1(1:7)="data.in" ! use default name
        nF=7
      elseif (nargin.eq.1) then
        call getarg(1,argin1) ! get input argument number 1
        nF=1 ! number of non-empty characters
        do while (argin1(nF:nF).ne." ")
          nF=nF+1
        enddo
        nF=nF-1
      else
        write(iwrit,*)
     +   'Too many input arguments on the command line ... exiting'
        stop
      endif

      if (argin1(1:4).eq."-gui") then! want to use the GUI of imageJ
        print*,"starting imageJ GUI for preprocessing"
        print*,"load an image with the menu"
        print*,"choose the frames, the BCs, the parameters, ..."
        print*,"write a text input file with all the info for NLACE"
        ! find a file name (nF+argin1)
        ! question could this text input file be reload for modifications.... great!
        print*,"exit the preprocessing gui (gray some button, slider,",
     &         " etc)"
        print*,"start running NLACE"
        ! example: just open and plot an image
        !"~/ImageJ/ImageJ.png"
        stop
      endif!if (argin1(1:4).eq."-gui")

c     READ DATA, INTIALIZE AND FILL IN SOME ARRAYS
      call data(nF,argin1)

c     if no output name has been specified, use the input file's name
      if (prefixOut(1:1).eq." ") then
        prefixOut(1:80)=argin1
      endif
c     check that the output file names chosen does not exist already; if so, modify output name
      nF=1
      do while (prefixOut(nF:nF).ne." ")
        nF=nF+1
      enddo
      nF=nF-1
      prefixOut((nF+1):(nF+4))='.ite' ! check with file .ite - for inverse solve
      inquire(file=prefixOut(1:(nF+4)),exist=ex)
      do while (ex)
        prefixOut((nF+1):(nF+5))='2.ite'
        nF=nF+1
        inquire(file=prefixOut(1:(nF+4)),exist=ex)
      enddo
      prefixOut((nF+1):(nF+4))='    ' ! remove .ite
      prefixOut((nF+1):(nF+6))='i0.vtk' ! check with file i0.vtk - for forward solve
      inquire(file=prefixOut(1:(nF+6)),exist=ex)
      do while (ex)
        prefixOut((nF+1):(nF+7))='2i0.vtk'
        nF=nF+1
        inquire(file=prefixOut(1:(nF+6)),exist=ex)
      enddo
      prefixOut((nF+1):(nF+6))='      ' ! remove i0.vtk

      if ((buildSymmetricMatrices).and.(solveMethod.gt.1)) then
        ! build full matrices even if symmetric. future work: implement the CG solver from MKL
        print*,"WARNING: fgmres solvers are meant for non-symmetric ",
     &         "matrices. Using with symmetric matrices may lead to a",
     &         " seg fault... what about implementing MKL-CG? ;-)"
      endif


      write(iwrit,'(/,A)') 
     +   ' ## Source/lmain.f: Input file read & memory filled'
      call timerMain(1)

c     THE MAIN PROGRAM, GRADIENT BASED MINIMIZATION or FORWARD PROBLEM
      if (problemType.eq.1) then
! 999    continue ! for inverse statics
        write(iwrit,*)" "
        write(iwrit,*)" ++ Running the inverse code ",
     +                "(optimization driver)"
        if (iopt.eq.2) then
          write(iwrit,*)" ++ L-BFGS-B selected"
        elseif (iopt.eq.4) then
          write(iwrit,*)" ++ ASA_CG selected"
        else
        endif
        write(iwrit,*)" "
        write(iwrit,*)" ++ iteration X0: loading the material with ",
     +          " the initial guess of the state vector"
        call driverInverse(iopt,niter,noptvar,bfgsM,noutput,prefixOut)
c        if (mod(niter,5).eq.0) goto 999 ! for inverse statics
c       the optimization algorithm has terminated before the NLACE criteria said so
        write(iwrit,*) "the opt algo has terminated before the NLACE",
     +      " criteria said so ..."
        g=-1234567.0d0! default intialization
        call result(currentOptIt+1,g)
      elseif (problemType.eq.2) then
        write(iwrit,*)"  ++ Running the forward code "
        call driverForward()
        write(iwrit,*)"## Source/lmain.f: Forward solve terminated"
        call timerMain(2)
      endif

c     FREE UP MEMORY
      call MEMORY_DEALLOC

      stop
      end


c******************************************************
      subroutine init_optvars(n,x)
c     initialize optimization variables for sd
      USE MAINMEM
      implicit none
      integer n
      double precision x(*)
      integer ipoin,iset,itemp
c------------------------------------------------------
      n = noptvar
      do ipoin = 1,npoin
        do iset = 1,nset_nodal
          itemp = ipoin_nod_grad(iset,ipoin)
          if (itemp.ne.0) then
            x(itemp) = adata_nodal(iset,ipoin)
          endif
        enddo
      enddo
      return
      end
      

c******************************************************
      subroutine init_optvars2(x,l,u)
c     initialize optimization variables for l-bfgs-b
      USE IOUNIT
      USE MAINMEM
      implicit none
      double precision lo,up
      double precision x(noptvar),l(noptvar),u(noptvar)
      integer ipoin,iset,itemp
c------------------------------------------------------
      do ipoin = 1,npoin
        do iset = 1,nset_nodal
          itemp = ipoin_nod_grad(iset,ipoin)
          lo=lowerb(iset)
          up=upperb(iset)
          if (itemp.ne.0) then
            if (iset.eq.1) then
              x(itemp) = adata_nodal(iset,ipoin)/scaleOptVar
              l(itemp) = lo/scaleOptVar
              u(itemp) = up/scaleOptVar
            else
              x(itemp) = adata_nodal(iset,ipoin)
              l(itemp) = lo
              u(itemp) = up
            endif
          endif
        enddo
      enddo
      return
      end
      

c******************************************************
      subroutine gradfun(x,f,g,itNum)
c     provide the gradient g and objective function value f at the state x
      USE IOUNIT
      USE MAINMEM
      implicit none
      double precision f,x(*),g(*)
      integer itNum! iteration number
      integer itNewton,imeas,lstep
      integer ipoin,iset,itemp! for scaling the gradient at the end of the function
      double precision stepL,accumStepL, ptmp
      double precision g1, g2, minX, maxX, minY, maxY
      double precision g2zones(optPnZones1,optPnZones2)
      integer nPtsZone(optPnZones1,optPnZones2)
      integer izone, jzone
      double precision xdiffn, gdiffn, gn! norms of vectors
      double precision maxfdiff, maxfdifft! for the stopping criterion
      logical stoppingHere
      character*40 stoppingMessage
      integer nReduceTimeStepAmplitude
      integer ieq
c------------------------------------------------------
c     some optimization software ask for f and g separately; to avoid recomputations
c      the previous results, f and g are stored in fmem and gmem; x is tored in xmem.
c      The first step is now to compare x and xmem; if identical or close enough, do 
c      not recompute but return the stored values fmem and gmem
      xdiffn=(x(1)-xmem(1))**2
      do ipoin=2,noptvar! ipoin is re-used
        xdiffn=xdiffn+(x(ipoin)-xmem(ipoin))**2
      enddo
      xdiffn=sqrt(xdiffn)
      if (xdiffn.lt.1.0e-15) then
        f=fmem(1)
        g(1:noptvar)=gmem
        write(iwrit,*) "gradfun: using the previous results to",
     &                 " avoid recomputing, xdiffn=",xdiffn
        return
      else
c        write(iwrit,*) "xdiffn=",xdiffn
c       determine the current iteration number (gradfun can be called several times per iteration)
        if (itNum.gt.lastOptIt) then
          currentOptIt=currentOptIt+1
          lastOptIt=itNum
          newIt=.true.
        else
          newIt=.false.
        endif

        if (countFGcalls.gt.0) then
          write(iwrit,*)"  "
          write(iwrit,*)" ++ iteration ",currentOptIt
        endif

        write(iwrit,'(A,A,I2,A)')'    -- Updating objective function ',
     $       'and gradient, ',nmeas,' measurements are used'
        call timer(0)

c       update adata_nodal and delta_adata_nodal using x
        call update_datn(x) 
  
c       plot image of the optimized fields every time gradfun is called
        if (monitorRun) then
          ! try to show the last 3-4 iterations ...
          print*,"if no imageJ window is open, open a imageJ window"
          print*,"split the image in nset_nodal*4 subplot"
          print*,"update the image by shifting old reconstructions to ",
     &           "the left"
          print*,"add the new reconstructions on the right column"
        endif

c       countFGcalls determines if step by step loading is used or the continuation method
        countFGcalls = countFGcalls + 1 
        itNewton=0
        convRes = 0.0d0! initialize to avoid a memory error
        convSol = 0.0d0! same
     

c JFD paralleliser ici avec mpi...

        do imeas = 1,nmeas ! for all measured fields
c         update idnum for the given measurement case
          call init_idnum(imeas)
          ifield = imeas
c         Create the matrices (csr sparse format)
          call inicsr
c          if (buildSymmetricMatrices) call reduceSymm JFD to remove if other approach works
          call timer(4)

c         Solve the primal problem (impose the BCs or increment the material properties)
          if (countFGcalls.eq.1) then ! first time through: finding balance with imposed bcs
            call forwardSolveEnforceBC(imeas)
          else ! every other time, a new balance is computed with updated material properties (material continuation)
            bc_step = 1.d0
            primal(:,:)=total_primal(imeas,:,:)
            lstep = 0
            accumStepL=0.0d0! accumulated step length (0<accumStepL<1)
            stepL=1.0d0/dble(ncontin)! step length
            nReduceTimeStepAmplitude = 0! number of times the time step amplitude can be reduced
            do while (accumStepL.lt.0.999999d0)
              
              lstep = lstep + 1
              if (timingpr) call timerProfiling(0)
              accumStepL=accumStepL+stepL
              if (accumStepL.gt.1.0d0) accumStepL=1.0d0! don't go beyond the limit
              call ContinIncrement(accumStepL)! modifies adata_nodal
              if (timingpr) call timerProfiling(7)
              convRes = 5*tol
              itNewton=0
              do while (convRes.gt.tol)! Newton iterations
                itNewton=itNewton+1
                negJac = .true.
                call stiff_residl(imeas) ! build rhs and tangent stiffness
                if (negJac) then
                  if (timingpr) call timerProfiling(2)
c                  call CSR_MATLAB (atang,jdiag,jcsr,aload,neqns,1)            
                  call residual(1) ! computes the convergence
                  if (timingpr) call timerProfiling(3)
                  call solve(solveMethod,itNewton) ! solve primal
                  if (timingpr) call timerProfiling(4)
                  call residual(2)
                  if (timingpr) call timerProfiling(5)
                  call update_primal ! update primal 
                  if (timingpr) call timerProfiling(6)
                  write(iwrit,'(A,I3,A,I3,A,1p,E10.4,A,1p,E10.4)')
     +             '       Loading step (mat):',lstep,', Newton it.:',
     +             itNewton,', residual(L2) = ',convRes,
     +             ', increment(L2) = ',convSol
                endif
                if ((itNewton.gt.50).or.(isnan(convRes)).or.
     +           (negJac.eqv. .false.)) then
                 if (nReduceTimeStepAmplitude.lt.20) then! this has not been tested yet - JFD 2011-04-12
                  write(iwrit,'(A,A)')'      Warning: Reducing the ',
     +               'step length to allow convergence'
                  nReduceTimeStepAmplitude = nReduceTimeStepAmplitude+1
                  itNewton=0
                  convRes=5*tol
c                 modify the step length and restart the computation
                  stepL=stepL/2.0d0
                  accumStepL=accumStepL-stepL
                  call ContinIncrement(accumStepL)! modifies adata_nodal
                  primal(:,:)=total_primal(imeas,:,:)
                 else
                  write(iwrit,'(A,A)') '      Warning: The specified ',
     +             'tolerance can not be reached. Going on as is...'
                  exit
                 endif
                endif
              enddo!do while (convRes.gt.tol)
             total_primal(imeas,1:npoin,1:mdofn)=primal(1:npoin,1:mdofn)! save the field obtained
            enddo! lstep
          endif
          call timer(1)

c         construct the rhs for the dual/adjoint problem
          call residl_dual(imeas)

c         transpose global stiffness matrix in pardiso format
          if (.not.buildSymmetricMatrices) call transpose_stiff
          call solve(solveMethod,10) ! solve dual - JFD: second argument is integer>5
          call update_dual(imeas) ! update dual
          call timer(2)
        
c         computing the gradient
          call calcgrad(imeas,f,g)
          call timer(3)
  
c         deallocate the memory allocated in inicsr
          call MEMORY_DEALLOCBCD

        enddo! imeas
      
        if (optPrecondition) then
          if (ndime.ne.2) then
            write(iwrit,*) "preconditioning not allowed for ndime.ne.2",
     &              " ... exiting"
            stop
          endif
c         scale the gradient and get the l2 norm of its components
          g1=0.0d0
          g2=0.0d0
          ! get the min and max x and y coordinate
          minX=coord(1,1)
          maxX=coord(1,1)
          minY=coord(2,1)
          maxY=coord(2,1)
          do ipoin=1,npoin
            if (coord(1,ipoin).lt.minX) minX=coord(1,ipoin)
            if (coord(1,ipoin).gt.maxX) maxX=coord(1,ipoin)
            if (coord(2,ipoin).lt.minY) minY=coord(2,ipoin)
            if (coord(2,ipoin).gt.maxY) maxY=coord(2,ipoin)
          enddo
c          write(iwrit,*)"minX=",minX,", maxX=",maxX,
c         &              ", minY=",minY,", maxY=",maxY
          g2zones(:,:)=0.0d0
          nPtsZone(:,:)=0
          do ipoin = 1,npoin
            ! get the x-y-z coordinates to compute L2-norms of gradient on sub-images
            izone=1+int(dble(optPnZones1)*(coord(1,ipoin)-minX)/
     &                (1.000001*(maxX-minX)))
            jzone=1+int(dble(optPnZones2)*(coord(2,ipoin)-minY)/
     &                (1.000001*(maxY-minY)))
       
c            write(iwrit,*)" node coord: (",coord(1,ipoin),", ",
c     &         coord(2,ipoin),"), zone=(",izone,", ",jzone,")"
            nPtsZone(izone,jzone)=nPtsZone(izone,jzone)+1;
            do iset = 1,nset_nodal
              itemp = ipoin_nod_grad(iset,ipoin)
              if (itemp.ne.0) then
                if (iset.eq.1) then
                  g(itemp)=g(itemp)*scaleOptVar
                  g1=g1+g(itemp)*g(itemp)
                else
                  g2=g2+g(itemp)*g(itemp)
                  g2zones(izone,jzone)=g2zones(izone,jzone)+
     &                   g(itemp)*g(itemp)
                endif
              endif
            enddo
          enddo
          g1=sqrt(g1)
          g2=sqrt(g2)
          g2zones=sqrt(g2zones)
          write(iwrit,*)"   -- For scaling: g1=",g1," and g2=",g2
          do ipoin=1,optPnZones2! print so that it be visually nice
            write(iwrit,*) g2zones(:,optPnZones2+1-ipoin)
          enddo
          do ipoin=1,optPnZones2! print so that it be visually nice
            write(iwrit,*) nPtsZone(:,optPnZones2+1-ipoin)
          enddo
          do ipoin=1,optPnZones2! print so that it be visually nice
            write(iwrit,*) g2zones(:,optPnZones2+1-ipoin)/
     &             dble(nPtsZone(:,optPnZones2+1-ipoin))
          enddo
        endif
    
c       stopping criteria for NLACE
        stoppingHere=.false.
        ! stopping criterion based on the max number of iterations specified in the input file
        if (niter.eq.currentOptIt) then
          stoppingHere=.true.
          stoppingMessage="maximum number of iterations is reached"
        endif
        ! stopping criterion based on the relative change in the objective function value
        maxfdiff=abs((f-fmem(1))/fmem(1))
        do ipoin=2,nfStored! ipoin is re-used
          maxfdifft=abs((f-fmem(ipoin))/fmem(ipoin))
          if (maxfdifft.gt.maxfdiff) maxfdiff=maxfdifft
        enddo
        if (maxfdiff.lt.1.0e-8) then
          stoppingHere=.true.
          stoppingMessage="objective function value has converged"
        endif
        ! stopping criterion based on the relative change in the gradient norm
        gdiffn=(g(1)-gmem(1))**2
        gn=g(1)**2
        do ipoin=2,noptvar! ipoin is re-used
          gdiffn=gdiffn+(g(ipoin)-gmem(ipoin))**2
          gn=gn+g(ipoin)**2
        enddo
        gdiffn=sqrt(gdiffn)
        gn=sqrt(gn)
        if (abs((gdiffn-gmemn)/gmemn).lt.1.0e-8) then
          stoppingHere=.true.
          stoppingMessage="gradient L2 norm has converged"
        endif

        write(iwrit,'(A,A,1p,E10.4,A,1p,E10.4)')"   -- Objective ",
     &     "function and gradient created, f=",f," and ||g||_2=",gn

c       Output of results 
        if ((currentOptIt.eq.nextOutput).or.stoppingHere) then
          if (stoppingHere) then
            call postprocess
          endif
          call result(currentOptIt,g)
          nextOutput = nextOutput+noutput
        endif
        
        call timerIteration(1)
  
        if (stoppingHere) then
          write(iwrit,*) " "
          write(iwrit,*) "    ",stoppingMessage,": stopping here"
          write(iwrit,'(A,1p,E10.4,A,1p,E10.4)')"     fmem(1)=",fmem(1),
     &               ", maxfdiff=",maxfdiff
          write(iwrit,'(A,1p,E10.4,A,1p,E10.4)')"     gmemn=",gmemn,
     &               ", gdiffn=",gdiffn
          write(iwrit,*) "    number of function calls = ",countFGcalls
          write(iwrit,*) " "
c         make the last screen print
          write(iwrit,*)"## Source/lmain.f: Optimization terminated",
     &                  " normally"
          close(16,status='keep')! close the cvg file
          call timerMain(2)
          if (iopt.eq.5) then
            close(10)! close algencan.out
            close(20)! close some file
            close(50)! close some other file
          endif
          stop
        endif 

        if (newIt) then
          
c          if((mod(currentOptIt+1,5).lt.0.01).and.newIt) then ! update coord, for inverse statics
c           coord(:,:) = coord_copy(:,:) - total_primal(1,:,:)
c           print*,'Update coordinates for inverse statics!'
c           return
c          endif

c         update the memroy values with the function's arguments for this computation
          xmem=x(1:noptvar)
          do ipoin=2,nfStored
            fmem(nfStored-ipoin+2)=fmem(nfStored-ipoin+1)
          enddo
          fmem(1)=f
          gmem=g(1:noptvar)
          gmemn=gn
        endif

        return
      endif!if (xdiff.lt.1.0e-13)
      end


c******************************************************
      subroutine forwardSolveEnforceBC(imeas)
c     loads the BCs incrementally to get the displacements
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer itNewton,imeas,lstep
      double precision stepL,accumStepL
      integer nReduceTimeStepAmplitude! number of times the time step amplitude can be reduced
c------------------------------------------------------

c     Solve the primal problem (impose the BCs or increment the material properties)
      lstep = 0
      accumStepL=0.0d0! accumulated step lengths (0<accumStepL<1)
      stepL=1.0d0/dble(lsteps)! step length
      nReduceTimeStepAmplitude = 0
      do while (accumStepL.lt.0.999999d0)
        lstep = lstep + 1
        if (timingpr) call timerProfiling(0)
        accumStepL=accumStepL+stepL
        if (accumStepL.gt.1.0d0) accumStepL=1.0d0! don't go beyond the limit
        bc_step = accumStepL
        call LoadIncrement(imeas,accumStepL)
        call ContinIncrementTrac(accumStepL)! modifies adata_nodal
        if (timingpr) call timerProfiling(1)
        convRes=5*tol
        itNewton=0
        do while (convRes.gt.tol)
          itNewton=itNewton+1
          call stiff_residl(imeas) ! build rhs and tangent stiffness
          if (timingpr) call timerProfiling(2)
c          call CSR_MATLAB (atang,jdiag,jcsr,aload,neqns,1)            
          call residual(1) ! computes the convergence
          if (timingpr) call timerProfiling(3)
          call solve(solveMethod,itNewton) ! solve primal
          if (timingpr) call timerProfiling(4)
          call residual(2)
          if (timingpr) call timerProfiling(5)
          call update_primal ! update primal 
          if (timingpr) call timerProfiling(6)
          write(iwrit,'(A,I3,A,I3,A,1p,E10.4,A,1p,E10.4)')
     +        '       Loading step (bc)=',lstep,', Newton it.=',
     +           itNewton,', residual(L2)= ',convRes,
     +           ', increment(L2)= ',convSol
c          exit
          if ((itNewton.gt.50).or.(isnan(convRes))) then! reduce the step length to allow convergence
           if (nReduceTimeStepAmplitude .lt. 10) then! this has not been tested yet... JFD 2011-04-12
            write(iwrit,'(A,A)')'      Warning:reducing the step',
     +           ' length to allow convergence'
            nReduceTimeStepAmplitude = nReduceTimeStepAmplitude+1
            itNewton=0
            convRes=5*tol
c           modify the step length and restart the computation
            primal(:,:)=total_primal(imeas,:,:)   ! original 
            stepL=stepL/2.0d0
            accumStepL=accumStepL-stepL
            call LoadIncrement(imeas,accumStepL)
            call ContinIncrementTrac(accumStepL)! modifies adata_nodal
           else
            write(iwrit,'(A,A)') '      Warning: The specified ',
     +          'tolerance can not be reached. Going on as is...'
            exit
           endif
          endif
c          total_primal(imeas,1:npoin,1:mdofn)=primal(1:npoin,1:mdofn)! save the field obtained
c          call result(itNewton,-123456.0d0)
        enddo!do while (convRes.gt.tol)
        total_primal(imeas,1:npoin,1:mdofn)=primal(1:npoin,1:mdofn)! save the field obtained
      enddo! lstep
      end


c*********************************************************
      subroutine init_idnum(imeas)
c     initializes idnum, primal and dual
c     idnum(npoin,mdofn): global dof number for each point
      USE MAINMEM
      implicit none
      integer imeas
      integer ipoin, nodfx,idofn
c---------------------------------------------------------
c     initialize to 0
      idnum(:,:)  = 0
      primal(:,:) = 0.0d0
      dual(:,:)   = 0.0d0

c     if bc is specified, update to the appropriate value

      do ipoin = 1,mpoinbc
        nodfx = ipoin_bc(imeas,ipoin)
        do idofn = 1,ldofn(nodfx)
          if (icode_bc(imeas,ipoin,idofn).ne.0) then
            idnum(nodfx,idofn) = icode_bc(imeas,ipoin,idofn)
            primal(nodfx,idofn) = value_bc(imeas,ipoin,idofn)
          endif
        enddo
      enddo
      
      return
      end


c******************************************************
      subroutine LoadIncrement(imeas,factor2)
c     increase boundary conditions incremental
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer imeas
      double precision factor2
      integer ipoin,nodfx,idofn
c--------------------------------------------------------
      do ipoin = 1,mpoinbc
        nodfx = ipoin_bc(imeas,ipoin)
        do idofn = 1,ldofn(nodfx)
          if (icode_bc(imeas,ipoin,idofn).ne.0) then
            primal(nodfx,idofn) = value_bc(imeas,ipoin,idofn)
     $                                      *factor2
          endif
        enddo
      enddo
      
      return
      end


c*********************************************************
      subroutine ContinIncrement(factor2)
c     the material parameters are increased at each step 
c     to reach the current material update step by step
      USE IOUNIT
      USE MAINMEM
      implicit none
      double precision factor2
      integer ipoin, itemp, iset
c--------------------------------------------------------
c     adata_nodal2 store the final state at the previous opt iteration
c     delta_adata_nodal is the total step to make at the current opt iteration
      adata_nodal(:,:) = delta_adata_nodal(:,:)*factor2
     $                   + adata_nodal2(:,:)
      do ipoin = 1,npoin
        do iset = 1,nset_nodal
          itemp = ipoin_nod_grad(iset,ipoin)
          if (itemp.ne.0) then
          endif
        enddo
      enddo

      return
      end


c*********************************************************
      subroutine ContinIncrementTrac(factor2)
c     the material parameters are increased at each step 
c     to reach the current material update step by step
      USE IOUNIT
      USE MAINMEM
      implicit none
      double precision factor2
      integer ipoin, itemp, iset
c--------------------------------------------------------
c     adata_nodal2 store the final state at the previous opt iteration
c     delta_adata_nodal is the total step to make at the current opt iteration
      do ipoin = 1,npoin
       do iset = 3,nset_nodal
        adata_nodal(iset,ipoin) = adata_nodal2(iset,ipoin)*factor2
       enddo
      enddo !for elem308 only
      do ipoin = 1,npoin
        do iset = 1,nset_nodal
          itemp = ipoin_nod_grad(iset,ipoin)
          if (itemp.ne.0) then
          endif
        enddo
      enddo

      return
      end


c*********************************************************
      subroutine update_datn(x) 
c     update the nodal data stored in adata_nodal using x
c       (the vector of optimized variables) and fill
c       delta_adata_nodal
      USE IOUNIT
      USE MAINMEM
      implicit none
      double precision x(*)
      double precision factor1
      integer ipoin, itemp, iset
c---------------------------------------------------------
c    for check
c      do ipoin = 1,npoin
c        do iset = 1,nset_nodal
c          Print*,"initial ",adata_nodal(iset,ipoin)
c          Print*,"current ",x(ipoin_nod_grad(iset,ipoin))
c        enddo
c      enddo
c      stop
c     end check
      do ipoin = 1,npoin
        do iset = 1,nset_nodal
          itemp = ipoin_nod_grad(iset,ipoin)
          if (itemp.ne.0) then
            adata_nodal2(iset,ipoin) = adata_nodal(iset,ipoin)! keep a copy of the previous state
            if (iset.eq.1) then! state where we want to go
              adata_nodal(iset,ipoin) = x(itemp)*scaleOptVar
            else
              adata_nodal(iset,ipoin) = x(itemp)
            endif
          endif
        enddo
      enddo
      
      delta_adata_nodal(:,:) = adata_nodal(:,:)-adata_nodal2(:,:)! step to make at the current opt iteration

      return
      end


c*********************************************************
      subroutine stiff_residl(imeas)
c     Compute the tangent stiffness matrix and the right hand 
c      side for the primal problem
      USE MAINMEM
      implicit none
      integer imeas
      integer ielem
      integer nnode, nevab
      double precision pelem(mprops)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,mnode)
      double precision elemdata_nodal(nset_nodal,mnode)
      double precision uelem(mdofn,mnode)
      double precision uelem_meas(mdofn,mnode) 
      integer ldnum(mevab)
c      integer TID,OMP_GET_THREAD_NUM! for debugging openmp
c---------------------------------------------------------
c     Initialize global variables
      aload(1:neqns) = 0.0d0! rhs/residual
      atang(1:nnz)   = 0.0d0! (tangent) stiffness
c     Initialize local variables
      nnode = 0! number of nodes in the element
      nevab = 0
      pelem(:)=0.0d0! material properties of the element
      estif(:,:) = 0.0d0! elemental (tangent) stiffness matrix
      eforc(:) = 0.0d0! elemental rhs/residual
      elnods(:) = 0! nodes numbers of the element
      xelem(:,:) = 0.0d0! coordinates of the nodes in the element
      elemdata_nodal(:,:) = 0.0d0! value of the parameters at the nodes of the element
      uelem(:,:) = 0.0d0! value of the displacements at the nodes of the element
      uelem_meas(:,:) = 0.0d0! value of the measured disp. at the nodes of the element
      ldnum(:) = 0

!$OMP PARALLEL DO
!$OMP& PRIVATE(ielem,pelem,nnode,nevab,elnods,xelem,elemdata_nodal,
!$OMP&   uelem,uelem_meas,elemvec_ndofn,estif,eforc,ldnum)
c!$OMP&   TID)
!$OMP& SHARED(nelem,imeas,lrefn,lnods,nset_nodal,adata_nodal,ldofn,
!$OMP&   idnum,primal,dual,meas,ndime,coord,mprops,props,mevab,jdiag,
!$OMP&   jcsr,buildSymmetricMatrices,bc_step)
!$OMP% DEFAULT(none)
!$OMP& SCHEDULE(STATIC,20)
!$OMP& REDUCTION(+:atang,aload)

c     element-level computations
      do ielem = 1,nelem
c       The following lines test openmp, uncomment for debug (don't forget the line above!)
c        TID = OMP_GET_THREAD_NUM()
c        arint*,'Thread ',TID,', imeas = ',imeas
c        print*,'Thread ',TID,', ielem = ',ielem,'/ ',nelem

c       Copy the information pertaining to the element into local variables
        call eloca (imeas,ielem,pelem,nnode,nevab,elnods,xelem,
     $          elemdata_nodal,uelem,uelem_meas,ldnum)

c       contribution from sources in the volume
        call elmlib (ielem,5,pelem,nnode,estif,eforc,elnods,xelem,
     $         elemdata_nodal,uelem,uelem_meas)
        call addv (eforc,aload,nevab,ldnum)

c       contributions from traction on boundary
c        call elmlib(ielem,5,pelem,nnode,estif,eforc,elnods,xelem,
c     $         elemdata_nodal,uelem,uelem_meas)
c        call addv (eforc,aload,nevab,ldnum)

c       residual contribution and its Hessian (tangent stiffness)
        call dpart (ielem,nnode)! JFD this step should already be ok ...
        call elmlib (ielem,3,pelem,nnode,estif,eforc,elnods,xelem,
     $         elemdata_nodal,uelem,uelem_meas)

        call subv  (eforc,aload,nevab,ldnum)
  
        if (buildSymmetricMatrices) then
c          call csr_addm (estif,atang,mevab,nevab,
c     $                      ldnum,jdiag,jcsr)! lower triangular
          call csr_addm2 (estif,atang,mevab,nevab,
     $                      ldnum,jdiag,jcsr)! upper triangular
        else
          call csr_addum (estif,atang,mevab,nevab,
     $                      ldnum,jdiag,jcsr)
        endif
      enddo
!$OMP END PARALLEL DO
        
c      Print*, "atang  ", atang
c      Print*, "aload  ", aload
c      stop
     
      call bfadd

      call weakadd(imeas)

      return
      end


c****************************************************************
      subroutine residual (rtask)
c     checks the convergence of the Newton-Iteration
c     taking L2-norm of the right-handside aload
      USE MAINMEM
      implicit none
      integer k, rtask
c----------------------------------------------------------------

      if (rtask.eq.1) then
        convRes = 0.0d0
        do k = 1, neqns
          convRes = convRes + aload(k)**2.0d0
        end do
        convRes=sqrt(convRes)

      else if (rtask.eq.2) then
        convSol = 0.0d0
        do k = 1, neqns
          convSol = convSol + aload(k)**2.0d0
        end do
        convSol=sqrt(convSol)
      endif

      return
      end


c*********************************************************
      subroutine residl_dual(imeas)
c     Compute the right hand side for the dual problem
      USE MAINMEM
      implicit none
      integer imeas
      integer ielem,i
      integer expon
c     Localized variables
      integer nnode,nevab
      double precision pelem(mprops)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,mnode)
      double precision elemdata_nodal(nset_nodal,mnode)
      double precision uelem(mdofn,mnode)
      double precision uelem_meas(mdofn,mnode)
      integer ldnum(mevab)
c---------------------------------------------------------
      aload(1:neqns) = 0.0d0
      expon = 2
c     Initialize local variables
      nnode = 0
      nevab = 0
      pelem(:) = 0.0d0
      estif(:,:) = 0.0d0
      eforc(:) = 0.0d0
      elnods(:) = 0
      xelem(:,:) = 0.0d0
      elemdata_nodal(:,:) = 0.0d0
      uelem(:,:) = 0.0d0
      uelem_meas(:,:) = 0.0d0
      forceM = 0.0d0
      ldnum(:) = 0
c JFD paralleliser ici avec openmp
      
c     Do-loop computes forceM needed for the following do-loop 
      do ielem = 1,nelem
        if (lrefn(ielem,1).eq.306) then
            call eloca_dual (imeas,ielem,expon,pelem,nnode,nevab,elnods,
     $            xelem,elemdata_nodal,uelem,ldnum)
            call elmlib (ielem,4,pelem,nnode,estif,eforc,elnods,xelem,
     $            elemdata_nodal,uelem,uelem_meas)           
        endif        
      enddo
      
c     element-level computations
      do ielem = 1,nelem
        call eloca_dual (imeas,ielem,expon,pelem,nnode,nevab,elnods,
     $         xelem,elemdata_nodal,uelem,ldnum)
         
c       Build the elemental RHS/residual for the dual problem
        call elmlib (ielem,6,pelem,nnode,estif,eforc,elnods,xelem,
     $         elemdata_nodal,uelem,uelem_meas)

        call addv  (eforc,aload,nevab,ldnum)
      enddo

      return
      end


c*********************************************************
      subroutine calcgrad(imeas,f,g)
c     Compute the gradient given the primal and the dual
      USE MAINMEM
      USE IOUNIT
      implicit none
      integer imeas
      double precision f, g(*)
      integer ielem,inode,ipoin,iset,itemp,expon
      integer nnode, nevab
      double precision pelem(mprops)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,mnode)
      double precision elemdata_nodal(nset_nodal,mnode)
      double precision uelem(mdofn,mnode)
      double precision uelem_meas(mdofn,mnode)
      integer ldnum(mevab)
      double precision suml21,suml22
c---------------------------------------------------------
      if (imeas.eq.1) then
        g(1:noptvar)   = 0.0d0
c        grad_nodal(:,:)  = 0.0d0! not used anymore
        dataMatch      = 0.0d0
        forceMatch     = 0.0d0
        regularization = 0.0d0
        dataMatchMem   = 0.0d0
        forceMatchMem  = 0.0d0
        regularizationMem= 0.0d0
        Gxi            = 0.0d0
      endif
      suml2grad1     = 0.0d0! JFD this is not optimal... should be done in this routine not in the elems
      suml2grad2     = 0.0d0
c     Initialize local variables
      nnode = 0
      nevab = 0
      pelem(:) = 0.0d0
      estif(:,:) = 0.0d0
      eforc(:) = 0.0d0
      elnods(:) = 0
      xelem(:,:) = 0.0d0
      elemdata_nodal(:,:) = 0.0d0
      uelem(:,:) = 0.0d0
      uelem_meas(:,:) = 0.0d0
      ldnum(:) = 0
      
      expon = 1
      
      
c JFD paralleliser ici avec openmp

c     element-level computations
      do ielem = 1,nelem
        call eloca_dual (imeas,ielem,expon,pelem,nnode,nevab,elnods,
     $         xelem,elemdata_nodal,uelem,ldnum)

        call elmlib (ielem,7,pelem,nnode,estif,eforc,elnods,xelem,
     $         elemdata_nodal,uelem,uelem_meas)


        do inode = 1,nnode
          ipoin = lnods(ielem,inode)! = elnods(inode)
          do iset = 1,nset_nodal
            itemp = ipoin_nod_grad(iset,ipoin)
            if (itemp.ne.0) then
c              grad_nodal(iset,ipoin) = grad_nodal(iset,ipoin)+
c     $              egrad(iset,inode)! not used anymore
              g(itemp) = g(itemp) + egrad(iset,inode)
            endif
          enddo! iset
        enddo! inode
      enddo! ielem

      forceMatch = forceMatch+0.5d0*Gxi*(forceM-Gforce)**2.0d0
      
      f = dataMatch+regularization+forceMatch

      suml21=0.0d0
      if (1.lt.nset_nodal) then
        suml22=0.0d0
      endif
      do inode=1,npoin
        suml21 = suml21 + adata_nodal(1,inode)*adata_nodal(1,inode)
        if (1.lt.nset_nodal) then
          suml22 = suml22 + adata_nodal(2,inode)*adata_nodal(2,inode)
        endif
      enddo
      suml21 = sqrt(suml21)
      suml2grad1 = sqrt(suml2grad1)
      if (1.lt.nset_nodal) then
        suml22 = sqrt(suml22)
        suml2grad2 = sqrt(suml2grad2)
      endif

c      if (((iopt.eq.2).and.(newIt.or.countFGcalls.le.2)).or.
c     $     (iopt.ne.2)) then
      if (newIt.or.countFGcalls.le.2) then
        if (imeas.eq.nmeas) then
          if (nset_nodal.eq.1) then
            write(16,'(1p,E13.5,1p,E13.5,1p,E13.5,1p,E19.11,1p,E13.5,
     $        1p,E13.5)') dataMatch-dataMatchMem,
     $        forceMatch-forceMatchMem,regularization-regularizationMem,
     $        f,suml21,suml2grad1
          else!if (nset_nodal.eq.2) then
            write(16,'(1p,E13.5,1p,E13.5,1p,E13.5,1p,E19.11,1p,E13.5,
     $        1p,E13.5,1p,E13.5,1p,E13.5)') dataMatch-dataMatchMem,
     $        forceMatch-forceMatchMem,regularization-regularizationMem,
     $        f,suml21,suml2grad1,suml22,suml2grad2
          endif
        else
          write(16,'(1p,E13.5,1p,E13.5)',advance="no")
     $      dataMatch-dataMatchMem,forceMatch-forceMatchMem
c       store the values of the relevant quantities for the next set of measures
          dataMatchMem=dataMatch
          forceMatchMem=forceMatch
           regularizationMem=regularization
        endif
      endif
      
      return
      end
     

c******************************************************
      subroutine eloca_dual (imeas,ielem,expon,pelem,nnode,nevab,elno,
     $   xel,elemd,uel,ldn)
c     localize element arrays
      USE MAINMEM
      implicit none 
      integer ielem,imeas
      integer ievab,inode,nodei,idofn,idime,iprop, iset
      integer matno
      integer ii,jdofn
      integer expon
      integer nnode,nevab
      double precision pelem(*)
      integer elno(*)! = elnods 
      double precision xel(ndime,*)! = xelem 
      double precision elemd(nset_nodal,*)! = elemdata_nodal
      double precision uel(mdofn,*)! = uelem
      double precision TT_nodal(mdofn,mdofn)
      double precision udiff_copy(mdofn)
      integer ldn(*)! = ldnum
c-------------------------------------------------------
      ievab = 0
      nevab = 0
      nnode = lrefn(ielem,2)
      uelem_diff(:,:) = 0.0d0
     
      do inode = 1,nnode
        nodei = lnods(ielem,inode)
        elno(inode) = lnods(ielem,inode)
        do iset = 1,nset_nodal
          elemd(iset,inode) = adata_nodal(iset,nodei)
          if (ipoin_nod_grad(iset,nodei).ne.0) 
     $         ielem_nod_grad(iset,inode) = 1
        enddo
         
c       localize equation numbers and dof values
        elemvec_ndofn(inode) = ldofn(nodei)
         
c     construct the nodal w_n*T-matrix
        ii = 0
        do idofn = 1,elemvec_ndofn(inode)
           TT_nodal(idofn,idofn) = fac_meas(imeas,nodei,idofn)
           do jdofn = idofn+1,elemvec_ndofn(inode)
              ii = ii+1 
              TT_nodal(idofn,jdofn) = fac_meas_nd(imeas,nodei,ii)
              TT_nodal(jdofn,idofn) = TT_nodal(idofn,jdofn)
           enddo
        enddo
c     end construct

        do idofn = 1,ldofn(nodei)
          ievab = ievab + 1
          ldn(ievab) = idnum(nodei,idofn)
          uel(idofn,inode) = primal(nodei,idofn)
          uelem_dual(idofn,inode) = dual(nodei,idofn)
          do jdofn = 1,ldofn(nodei)
             uelem_diff(idofn,inode) = uelem_diff(idofn,inode)+ 
     $            TT_nodal(idofn,jdofn)*
     $            (primal(nodei,jdofn)-
     $             meas(imeas,nodei,jdofn))
          enddo
        enddo
        if (expon.eq.2) then
          udiff_copy(1:mdofn) = uelem_diff(1:mdofn,inode)
          uelem_diff(1:mdofn,inode) = 0.0d0
          do idofn=1,ldofn(nodei)
            do jdofn = 1,ldofn(nodei)
               uelem_diff(idofn,inode) = uelem_diff(idofn,inode)+ 
     $              TT_nodal(jdofn,idofn)*udiff_copy(jdofn)
             enddo
          enddo
        endif
         
c       localize nodal coordinates
        do idime = 1,ndime
          xel(idime,inode) = coord(idime,nodei)
        enddo
      enddo

      nevab = ievab
c     localize element properties
      matno = lrefn(ielem,3)
      do iprop = 1,mprops! 1,nprop !(with nprop=lprop(matno))
        pelem(iprop) = props(matno,iprop)
      enddo
      
      return
      end
      
      
c******************************************************
      subroutine eloca (imeas,ielem,pelem,nnode,nevab,elno,xel,
     $      elemd,uel,uel_meas,ldn)
c     localize element arrays
      USE MAINMEM
      implicit none 
      integer imeas,ielem
      integer nnode,nevab
      double precision pelem(*)
      integer elno(*)! = elnods
      double precision xel(ndime,*)! = xelem
      double precision elemd(nset_nodal,*)! = elemdata_nodal
      double precision uel(mdofn,*)! = uelem
      double precision uel_meas(mdofn,*)! = uelem_meas
      integer ldn(*)! = ldnum
      integer ievab,inode,nodei,idofn,idime,iprop, iset
      integer matno
c-------------------------------------------------------
      ievab = 0
      nevab = 0
      nnode = lrefn(ielem,2)
      
      do inode = 1,nnode
        nodei = lnods(ielem,inode)
        elno(inode) = lnods(ielem,inode)! nodes numbers in the element
        do iset = 1,nset_nodal
          elemd(iset,inode) = adata_nodal(iset,nodei)
        enddo
         
c       localize equation numbers and dof values
c        elemvec_ndofn(inode) = ldofn(nodei)! incompatible with openmp, make it local
        do idofn = 1,ldofn(nodei)
          ievab = ievab + 1
          ldn(ievab) = idnum(nodei,idofn)
          uel(idofn,inode) = primal(nodei,idofn)
          uel_meas(idofn,inode) = meas(imeas,nodei,idofn)
        enddo
c       localize nodal coordinates
        do idime = 1,ndime
          xel(idime,inode) = coord(idime,nodei)
        enddo
      enddo
      
      nevab = ievab
         
c     localize element properties
      matno = lrefn(ielem,3)
      do iprop = 1,mprops! 1,nprop !(with nprop=lprop(matno))
        pelem(iprop) = props(matno,iprop)
      enddo
      
      return
      end
      
      
c********************************************************
      subroutine update_primal
c     update solution vector
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ipoin,idofn,k
c--------------------------------------------------------
      do ipoin = 1,npoin
        do idofn = 1,ldofn(ipoin)
          k = idnum(ipoin,idofn)
          if (k.ne.0) then
            primal(ipoin,idofn) = primal(ipoin,idofn) + aload(k)
          endif
        enddo
      enddo

      do idofn = 1,4
c        write(iwrit, *) "idofn=", idofn
        do ipoin = 1,npoin
c            write(iwrit, *) "primal(ipoin,idofn)=",
c     &                       primal(ipoin,idofn),
c     &                      " ipoin=", ipoin, " idofn=", idofn
c            write(iwrit, *) primal(ipoin,idofn)
        enddo
      enddo
      return
      end


c********************************************************
      subroutine transpose_stiff
c     transpose the stiffness matrix in pardiso format
      USE MAINMEM
      implicit none
      integer i, j, k, l
c--------------------------------------------------------
c      k = 0
c      jdiag_temp(:) = 0 
c      atang_temp(:) = 0.0d0
c      jcsr_temp(:)  = 0
c      jdiag_temp(1) = 1
c      do j = 1,neqns
c        do i = 1,nnz
c          if (jcsr(i).eq.j) then
c            k = k+1
c            atang_temp(k) = atang(i)
c            l=0
c            do while (jdiag(l+1).le.i)
c                    l = l+1
c            enddo
c            jcsr_temp(k)  = l
c          endif
c        enddo
c        jdiag_temp(j+1) = k + 1
c      enddo
c      atang(1:nnz) = atang_temp(1:nnz)
c      jcsr(1:nnz) = jcsr_temp(1:nnz)
c      jdiag(1:(neqns+1)) = jdiag_temp(1:(neqns+1))


c      Alternative way, more efficient and faster
       jdiag_temp(:) = 0
       atang_temp(:) = 0.0d0
       jcsr_temp(:)  = 0

cc       do j = 1,nnz
cc         jdiag_temp(jcsr(j)) = jdiag_temp(jcsr(j))+1
cc       end do
cc       jdiag_temp(2:(neqns+1)) = jdiag_temp(1:neqns)
cc      jdiag_temp(1) = 1

cc       do j = 1,neqns
cc         jdiag_temp(j+1) = jdiag_temp(j)+jdiag_temp(j+1)
cc       end do

cc       jdiag(1:(neqns+1)) = jdiag_temp(1:(neqns+1))
cc       jdiag_temp(:) = 0

       do j = 1,nnz
         k = jdiag(jcsr(j))+jdiag_temp(jcsr(j))
         jdiag_temp(jcsr(j)) = jdiag_temp(jcsr(j))+1
cc         jcsr_temp(k)  = jcsr(j)
         atang_temp(k) = atang(j)
       end do
cc       jcsr(1:nnz)  = jcsr_temp(1:nnz)
       atang(1:nnz) = atang_temp(1:nnz)

      return
      end


c********************************************************
      subroutine update_dual(imeas)
c     update the dual vector
      USE MAINMEM
      implicit none
      integer imeas
      integer ipoin,idofn,k
c--------------------------------------------------------
      do ipoin = 1,npoin
        do idofn = 1,ldofn(ipoin)
          k = idnum(ipoin,idofn)
          if (k.ne.0) then
            dual(ipoin,idofn) = aload(k)
          endif
        enddo
      enddo
c     update the total dual
      total_dual(imeas,1:npoin,1:mdofn) = dual(1:npoin,1:mdofn)
      
      return
      end

c********************************************************
      subroutine bfadd 
c     add the body forces into the rhs 
      USE MAINMEM
      implicit none
      integer ipoin,idofn,nevab,ldnum
      double precision bforc
c--------------------------------------------------------
      nevab = 1
      bforc = 0.0d0 
      do ipoin = 1,npoin
        do idofn = 1,2
          ldnum= idnum(ipoin,idofn) 
          bforc = bodyforce(ipoin,idofn)
          call addv (bforc,aload,nevab,ldnum)  
        enddo
      enddo
      
      return
      end

c********************************************************
      subroutine weakadd(imeas)
c     add the weak springs 
      USE MAINMEM
      implicit none
      integer ipoin,idofn,k,nevab,ldnum
      integer imeas,nodfx
      double precision fpart,fpart2
c--------------------------------------------------------
      nevab = 1 
      fpart = 0.0d0
      fpart2 = 0.0d0
      ldnum = 0 
      do ipoin = 1,mpoinbc
        nodfx = ipoin_ws(imeas,ipoin)
        do idofn = 1,ldofn(nodfx)
          if (icode_ws(imeas,ipoin,idofn).ne.0) then
          ldnum = idnum(nodfx,idofn)
          fpart = value_ws(imeas,ipoin,idofn) * delta_ws 
          fpart2 = primal(nodfx,idofn) * delta_ws
c         residual contribution
          call addv (fpart,aload,nevab,ldnum)  
          call subv (fpart2,aload,nevab,ldnum) 
c         stiffness contribution
          call csr_addm2 (delta_ws,atang,mevab,nevab,
     $                      ldnum,jdiag,jcsr)
          endif
        enddo
      enddo
      return
      end
c*********************************************************
      subroutine dpart (ielem,nnode)
c     partition and zero the element vector 
      USE MAINMEM
      implicit none
      integer ielem
      integer nnode
      integer ievab,idofn,k,inode,nodei
c----------------------------------------------------------
      ievab = 0
      do inode = 1,nnode
        nodei = lnods(ielem,inode)
        elemvec_ndofn(inode) = ldofn(nodei)
c        do idofn = 1,ldofn(nodei)! JFD: do not know why it was commented ... ldnum is no longer allocatable
c          ievab = ievab + 1
c          k = ldnum(ievab)
c          if (k.ne.0) then
c            uelem(idofn,inode) = 0.d0
c          endif
c        enddo
      enddo
      return
      end
c***********************************************************
      subroutine postprocess 
c     Post process data with element 307
      USE MAINMEM
      implicit none
      integer imeas,ielem, count
      integer nnode,nevab,ldn(mevab)
      double precision pelem(mprops),estif(mevab,mevab)
      integer elno(mnode)
      double precision eforc(mevab),xel(ndime,mnode)
      double precision elemdata_nodal(nset_nodal,mnode)
      double precision uel(mdofn,mnode)
      double precision uel_meas(mdofn,mnode)
      
      pforce  = 0.0d0 
      moment = 0.0d0
      count = 1
      imeas = 1 !!!!might need to change later!!!
      do ielem=1, nelem
        if (lrefn(ielem,1).eq.307) then
          count=0
          call eloca (imeas,ielem,pelem,nnode,nevab,elno,xel,
     $          elemdata_nodal,uel,uel_meas,ldn)
          call elmlib (ielem,4,pelem,nnode,estif,eforc,elno,xel,
     $            elemdata_nodal,uel,uel_meas)           
        endif
      enddo

      if (count.eq.0) then
        MuFactor=abs(ForceTotal/pforce)
        adata_nodal(2,:)=MuFactor*adata_nodal(2,:)
      endif

      return
      end
c***********************************************************
      subroutine termin (ifail)
c     termination of execution due to failure modes
      USE MAINMEM
      USE IOUNIT
      implicit none
      integer ifail
c-----------------------------------------------------------
      
      write(iwrit,800)
 800  format(1x/17x,30('*')/17x,'*',28(' '),'*'/
     +     17x,'*    execution terminated    *'/
     +     17x,'*',28(' '),'*'/17x,30('*'))
      
      write(iwrit,*)ifail

      stop
      end

