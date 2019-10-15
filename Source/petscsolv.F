c**********************************************************
      subroutine petscsolve
c     solve the sytem of equations
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"

c Petsc 2.2.0 Sles Gone
c#include "include/finclude/petscsles.h"

      Vec              x,b,u
      Mat              A 
      KSP              ksp
      PC               pc

c     Petsc 2.2.0 Sles Gone
c      SLES             sles

      double precision norm,tol,normb
      integer          ierr,i,n,col(200),its,flg,size,rank
      integer          nval,j,im1
      integer          ipoin,idofn,ii
      integer          irow(NEQNS)
      PetscScalar           none,one,valone,value(200)
      PetscScalar           xx_v(1)
      PetscOffset      xx_i
      PetscScalar        omega

c     
c----------------------------------------------------------
      
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      
      call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
      if (size .ne. 1) then
         call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
         if (rank .eq. 0) then
            write(6,*) 'This is a uniprocessor example only!'
         endif
         SETERRQ(1,' ',ierr)
      endif
      none = -1.0
      one  = 1.0
      n    = NEQNS
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr)
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Compute the matrix and right-hand-side vector that define
!         the linear system, Ax = b.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create matrix.  When using MatCreate(), the matrix format can
!  be specified at runtime.

c      call MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,A,   &
c     &               ierr)
c      call MatSetFromOptions(A,ierr)

      call MatCreateSeqAIJ(PETSC_COMM_WORLD,n,n,100,PETSC_NULL_INTEGER,
     $     A,ierr)


      do i = 1,n
         im1 = i-1
         nval = 0
         do j = jdiag(i),jdiag(i+1)-1
            nval = nval+1
            col(nval) = jcsr(j)-1
            value(nval) = atang(j)
c            write(*,*) i
c            write(10,*)i,col(nval)+1,value(nval)
         enddo
         call MatSetValues(A,1,im1,nval,col,value,INSERT_VALUES,ierr)
      enddo

            

!  Assemble matrix.
!   - Note that MatSetValues() uses 0-based row and column numbers
!     in Fortran as well as in C (as set here in the array "col").
c
c      value(1) = -1.0
c      value(2) = 2.0
c      value(3) = -1.0
c      do 50 i=1,n-2
c         col(1) = i-1
c         col(2) = i
c         col(3) = i+1
c         call MatSetValues(A,1,i,3,col,value,INSERT_VALUES,ierr)
c  50  continue
c      i = n - 1
c      col(1) = n - 2
c      col(2) = n - 1
c      call MatSetValues(A,1,i,2,col,value,INSERT_VALUES,ierr)
c      i = 0
c      col(1) = 0
c      col(2) = 1
c      value(1) = 2.0
c      value(2) = -1.0
c      call MatSetValues(A,1,i,2,col,value,INSERT_VALUES,ierr)
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      
!  Create vectors.  Note that we form 1 vector from scratch and
!  then duplicate as needed.

      call VecCreate(PETSC_COMM_WORLD,x,ierr)
      call VecSetSizes(x,PETSC_DECIDE,n,ierr)
      call VecSetFromOptions(x,ierr)
      call VecDuplicate(x,b,ierr)
      call VecDuplicate(x,u,ierr)

!  Set exact solution; then compute right-hand-side vector.

c      call VecSet(one,u,ierr)
c      call MatMult(A,u,b,ierr)

      do i = 1,n
         im1 = i-1
         valone = aload(i)
c         write(20,*)valone
         call VecSetValue(b,im1,valone,INSERT_VALUES,ierr)
      enddo
c     call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)

c     set the inital guess to the previous solution

c      do ipoin = 1,npoin
c         do idofn = 1,ldofn(ipoin)
c            ii = idnum(ipoin,idofn)
c            if (ii.ne.0) then
c               im1 = ii-1
c               valone = primal(ipoin,idofn)
c               call VecSetValue(x,im1,valone,INSERT_VALUES,ierr)
c            endif
c         enddo
c      enddo

c     end set value
               
               

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!          Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create linear solver context
      
      write(*,*)'solving now'

c     SLES GONE PETSC 2.2.0
c     call SLESCreate(PETSC_COMM_WORLD,sles,ierr)

      call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.

c     SLES GONE PETSC 2.2.0
c     call SLESSetOperators(sles,A,A,DIFFERENT_NONZERO_PATTERN,ierr)

      call KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)
      
!  Set linear solver defaults for this problem (optional).
!   - By extracting the KSP and PC contexts from the SLES context,
!     we can then directly directly call any KSP and PC routines
!     to set various options.
!   - The following four statements are optional; all of these
!     parameters could alternatively be specified at runtime via
!     SLESSetFromOptions();

c     SLES GONE!
c      call SLESGetKSP(sles,ksp,ierr)
c      call SLESGetPC(sles,pc,ierr)

      call KSPGetPC(ksp,pc,ierr)
      
c      call PCSetType(pc,PCJACOBI,ierr)
c      call KSPSetType(ksp,KSPBCGS,ierr)
      call KSPSetType(ksp,KSPCG,ierr)

c      omega = 1.7
c      call PCSetType(pc,PCEISENSTAT,ierr)
c      call PCEisenstatSetOmega(pc,omega,ierr)
      call PCSetType(pc,PCILU,ierr)
      call PCILUSetLevels(pc,1,ierr)
c      call PCSetType(pc,PCLU,ierr)
      tol = 1.d-6
      call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION,     &
     &     PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
      
!  Set runtime options, e.g.,
!      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  SLESSetFromOptions() is called _after_ any other customization
!  routines.

c     SLES GONE
c      call SLESSetFromOptions(sles,ierr)

      call KSPSetFromOptions(ksp,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c     SLES GONE
c      call SLESSolve(sles,b,x,its,ierr)
      
c      call KSPSetRhs(ksp,b,ierr)
c      call KSPSetSolution(ksp,x,ierr)
      call KSPSolve(ksp,b,x,ierr)


!  View solver info; we could instead use the option -sles_view

c      call SLESView(sles,PETSC_VIEWER_STDOUT_WORLD,ierr)
      
c      call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Check the error

c      call VecAXPY(none,u,x,ierr)
      write(*,*)'checking solve .....'
      call MatMult(A,x,u,ierr)
      call VecNorm(b,NORM_2,normb,ierr)
      call VecAXPY(none,u,b,ierr)
      call VecNorm(b,NORM_2,norm,ierr)
      call KSPGetIterationNumber(ksp,its,ierr)

      norm = norm/normb
      if (norm .gt. 1.e-12) then
        write(6,100) norm,its
      else
        write(6,200) its
      endif
 100  format(' Relative Norm of error = ',e10.4,',  Iterations = ',i5)
 200  format(' Relative Norm of error < 1.e-12,Iterations = ',i5)


      call VecGetArray(x,xx_v,xx_i,ierr)
      do i = 1,NEQNS
         ALOAD(i) = xx_v(i+xx_i)
c         write(*,*)'ALOAD==',ALOAD(i)
      enddo
      call VecRestoreArray(x,xx_v,xx_i,ierr)
 
!  Free work space.  All PETSc objects should be destroyed when they
!  are no longer needed.

      call VecDestroy(x,ierr)
      call VecDestroy(u,ierr)
      call VecDestroy(b,ierr)
      call MatDestroy(A,ierr)
C      call SLESDestroy(sles,ierr)
      call KSPDestroy(ksp,ierr)
c      write(*,*)'here-------------'
c      call PetscFinalize(ierr)

      end
