c**********************************************************
      subroutine solve(solveMethod,itN)
c     Select the solver for the sytem of equations
      implicit none 
      integer solveMethod,itN
c----------------------------------------------------------
c     NOTE: the matrix is in CSR format in variables jdiag, jcsr and atang
c           the rhs is in aload
c           the solution should and will be put in aload after solving
c           tvect can be used a buffer or a copy of aload

      if (solveMethod.eq.0) then
c       Use PETSC solver (not test since 2009)
c        call petscsolve
      elseif (solveMethod.eq.1) then
c       Use the solver pardiso (U. of Basel or Intel MKL)
        call pardisosolve
      elseif (solveMethod.eq.2) then
c       Use the RCI F-GMRES solver + ILU0 preconditioning from Intel MKL
c        call fgmresilu0solve
      elseif (solveMethod.eq.3) then
c       Use the RCI F-GMRES solver + ILU0 preconditioning from Intel MKL
c        call fgmresilutsolve
      endif

      return
      end

