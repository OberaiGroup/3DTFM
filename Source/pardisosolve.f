********************************************************************************
*                              INTEL CONFIDENTIAL
*   Copyright(C) 2004-2010 Intel Corporation. All Rights Reserved.
*   The source code contained  or  described herein and all documents related to
*   the source code ("Material") are owned by Intel Corporation or its suppliers
*   or licensors.  Title to the  Material remains with  Intel Corporation or its
*   suppliers and licensors. The Material contains trade secrets and proprietary
*   and  confidential  information of  Intel or its suppliers and licensors. The
*   Material  is  protected  by  worldwide  copyright  and trade secret laws and
*   treaty  provisions. No part of the Material may be used, copied, reproduced,
*   modified, published, uploaded, posted, transmitted, distributed or disclosed
*   in any way without Intel's prior express written permission.
*   No license  under any  patent, copyright, trade secret or other intellectual
*   property right is granted to or conferred upon you by disclosure or delivery
*   of the Materials,  either expressly, by implication, inducement, estoppel or
*   otherwise.  Any  license  under  such  intellectual property  rights must be
*   express and approved by Intel in writing.
*
********************************************************************************
*   Content : MKL PARDISO Fortran-77 example
*
********************************************************************************
C----------------------------------------------------------------------
C Example program to show the use of the "PARDISO" routine
C for symmetric linear systems
C---------------------------------------------------------------------
C This program can be downloaded from the following site:
C www.pardiso-project.org
C
C (C) Olaf Schenk, Department of Computer Science,
C University of Basel, Switzerland.
C Email: olaf.schenk@unibas.ch
C
C---------------------------------------------------------------------

c      PROGRAM pardiso_sym
c      IMPLICIT NONE
c      include 'mkl_pardiso.f77'

      subroutine pardisosolve 
      USE IOUNIT! access the file ID for screen-printing
      USE MAINMEM! access atang, jdiag, jcsr, neqns, aload, tvect
      IMPLICIT NONE
C.. Internal solver memory pointer for 64-bit architectures
C.. INTEGER*8 pt(64)
C.. Internal solver memory pointer for 32-bit architectures
C.. INTEGER*4 pt(64)
C.. This is OK in both cases
      INTEGER*8 pt(64)
C.. All other variables
      INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
c      INTEGER n -> neqns
      INTEGER iparm(64)
c      INTEGER ia(9) -> jdiag
c      INTEGER ja(18) -> jcsr
c      REAL*8 a(18) -> atang
c      REAL*8 b(8) -> aload
c      REAL*8 x(8) -> tvect
      INTEGER i, idum
      REAL*8 waltime1, waltime2, ddum
C.. Fill all arrays containing matrix data.
c      DATA n /8/
      DATA nrhs /1/, maxfct /1/, mnum /1/
c      DATA ia /1,5,8,10,12,15,17,18,19/
c      DATA ja
c     1 /1,  3,    6,7,
c     2    2,3,  5,
c     3      3,        8,
c     4        4,    7,
c     5          5,6,7,
c     6            6,  8,
c     7              7,
c     8                8/
c      DATA a
c     1 /7.d0,     1.d0,          2.d0,7.d0,
c     2       -4.d0,8.d0,     2.d0,
c     3            1.d0,                    5.d0,
c     4                 7.d0,     9.d0,
c     5                      5.d0,1.d0,5.d0,
c     6                           -1.d0,     5.d0,
c     7                                11.d0,
c     8                                     5.d0/

c      integer omp_get_max_threads
c      external omp_get_max_threads

C..
C.. Set up PARDISO control parameter
C..
      if (timingpr) call timerPardiso(0)
      do i = 1, 64
         iparm(i) = 0
      end do
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
c      iparm(3) = omp_get_max_threads()!1 ! numbers of processors
      iparm(3) = nCores ! numbers of processors
      iparm(4) = 0 ! no iterative-direct algorithm
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 9 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) = 0 ! Output: Numbers of CG Iterations
      error = 0 ! initialize error flag
      msglvl = 0! JFD 1 ! print statistical information
      if (buildSymmetricMatrices) then
        mtype = -2 ! symmetric, indefinite
      else
        mtype = 11 ! unsymmetric real
      endif
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
      do i = 1, 64
         pt(i) = 0
      end do
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
      phase = 11 ! only reordering and symbolic factorization
      CALL pardiso (pt,maxfct,mnum,mtype,phase,neqns,atang,jdiag,jcsr,
     1 idum, nrhs,iparm,msglvl,ddum,ddum,error)
c      WRITE(*,*) 'Reordering completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         STOP 1
      END IF
      if (timingpr) call timerPardiso(1)
c      WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
c      WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
C.. Factorization.
      phase = 22 ! only factorization
      CALL pardiso (pt,maxfct,mnum,mtype,phase,neqns,atang,jdiag,jcsr,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
c      WRITE(*,*) 'Factorization completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         STOP 1
      ENDIF
      if (timingpr) call timerPardiso(2)
C.. Back substitution and iterative refinement
      iparm(8) = 2 ! max numbers of iterative refinement steps
      phase = 33 ! only factorization
c      do i = 1, neqns
c         aload(i) = 1.d0
c      end do
      tvect(:)=0.0d0
      CALL pardiso (pt,maxfct,mnum,mtype,phase,neqns,atang,jdiag,jcsr,
     1 idum, nrhs, iparm, msglvl, aload, tvect, error)
c      WRITE(*,*) 'Solve completed ... '
c      WRITE(*,*) 'The solution of the system is '
      DO i = 1, neqns
         aload(i) = tvect(i)
c         WRITE(*,*) ' tvect(',i,') = ', tvect(i)
      END DO
      if (timingpr) call timerPardiso(3)
C.. Termination and release of memory
      phase = -1 ! release internal memory
      CALL pardiso (pt,maxfct,mnum,mtype,phase,neqns,ddum,idum,idum,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
      if (timingpr) call timerPardiso(4)
      END
