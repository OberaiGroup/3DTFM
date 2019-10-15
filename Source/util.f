c**********************************************************
      subroutine addv (e,a,nevab,ldn)
c     assemble element vector { e } into global vector { a }
      implicit none
      integer nevab,ldn(*)
      double precision e(*),a(*)
      integer ievab,k
c------------------------------------------------------------
      do ievab = 1,nevab
        k = ldn(ievab)
        if (k.gt.0) then
          a(k) = a(k) + e(ievab)
        endif
      enddo

      return
      end


c**********************************************************
      subroutine subv (e,a,nevab,ldn)
c     assemble element vector { -e } into global vector { a }
      implicit none
      integer nevab,ldn(*)
      double precision e(*),a(*)
      integer ievab,k
c-----------------------------------------------------------
      do ievab = 1,nevab
        k = ldn(ievab)
        if (k.gt.0) then
          a(k) = a(k) - e(ievab)
        endif
      enddo

      return
      end
      

c**************************************************************
      subroutine csr_addm (estif,astif,mevab,nevab,ldn,icsr,jcsr)
c     assemble symmetric, real element matrix into global 
c                array stored in csr form
c     makes a lower triangular matrix
      implicit none
      integer mevab,nevab
      integer ldn(*),icsr(*),jcsr(*)
      double precision estif(mevab,*),astif(*)
      integer ievab,i,jst,jen,j,jpos,jevab
c---------------------------------------------------------------
      do ievab = 1,nevab
        i = ldn(ievab)
        if (i.ne.0) then
          jst = icsr(i)
          jen = icsr(i+1)-1
          do jevab = 1,nevab
            j = ldn(jevab)
            if ((i.ge.j).and.(j.ne.0)) then

              jpos = jst
              do while ( (jcsr(jpos).ne.j).and.(jpos.le.jen))
                jpos = jpos+1
              enddo
              if ((jcsr(jpos).ne.j).and.(jpos.eq.jen)) then
                write(*,*)' Error in c_csr_addm'
              endif
              astif(jpos) = astif(jpos)+estif(ievab,jevab)

            endif
          enddo
        endif
      enddo

      return
      end


c**************************************************************
      subroutine csr_addm2 (estif,astif,mevab,nevab,ldn,icsr,jcsr)
c     assemble symmetric, real element matrix into global 
c                array stored in csr form
c     makes a upper triangular matrix
      implicit none
      integer mevab,nevab
      integer ldn(*),icsr(*),jcsr(*)
      double precision estif(mevab,*),astif(*)
      integer ievab,i,jst,jen,j,jpos,jevab
c---------------------------------------------------------------
      do ievab = 1,nevab
        i = ldn(ievab)
        if (i.ne.0) then
          jst = icsr(i)
          jen = icsr(i+1)-1
          do jevab = 1,nevab
            j = ldn(jevab)
            if ((i.le.j).and.(j.ne.0)) then

              jpos = jst
              do while ( (jcsr(jpos).ne.j).and.(jpos.le.jen))
                jpos = jpos+1
              enddo
              if ((jcsr(jpos).ne.j).and.(jpos.eq.jen)) then
                write(*,*)' Error in c_csr_addm2'
              endif
              astif(jpos) = astif(jpos)+estif(ievab,jevab)

            endif
          enddo
        endif
      enddo

      return
      end


c**************************************************************
      subroutine csr_addum (estif,astif,mevab,nevab,ldn,icsr,jcsr)
c     assemble unsymmetric, real  element matrix into real 
c        global array stored in csr form
      implicit none
      integer mevab,nevab
      integer ldn(*),icsr(*),jcsr(*)
      double precision estif(mevab,*),astif(*)
      integer ievab,i,jst,jen,j,jpos,jevab
c---------------------------------------------------------------
      do ievab = 1,nevab
        i = ldn(ievab)
        if (i.ne.0) then
          jst = icsr(i)
          jen = icsr(i+1)
c          l = jdiag(j) - j
          do jevab = 1,nevab
            j = ldn(jevab)
            jpos = jst
            if ( j.ne.0) then

              jpos = jst
              do while ( (jcsr(jpos).ne.j).and.(jpos.le.jen))
                jpos = jpos+1
              enddo
              if ((jcsr(jpos).ne.j).and.(jpos.eq.jen)) then
                write(*,*)' Error in c_csr_addum'
              endif
              astif(jpos) = astif(jpos)+estif(ievab,jevab)
                   
            endif
          enddo
        endif
      enddo

      return
      end


c*************************************************************
      subroutine CSR_MATLAB (atang,icsr,jcsr,aload,neqns,iflag)
C     Dump Assembled coefficient matrix A in coordinate sparse 
C       storage scheme --- (irow,icol,A(irow,icol)) --- to the
C       following files: 
C                       AssA.m : Assembled A matrix
C                       Ival.m : Row number of each element in A
C                       Jval.m : Col number of each element in A
C     Arguments:
C        c_atang: coefficient matrix in SSR format (input)
C        icsr   : row indices of entries in c_atang (input)
C        jcsr   : col indices of entries in c_atang (input)
C        neqns  : number of equations (input)
C        iflag  : iflag = 1 => CSR, 2 => CSC.
      implicit none
      double precision atang(*),aload(*)
      integer icsr(*),jcsr(*),neqns,iflag
      integer ii,jj,ibeg,iend
c--------------------------------------------------------------
      open(95,file='mat.m',status='unknown')
      open(96,file='rhs.m',status='unknown')
      write(95,*)'%****************************************'
      write(95,*)'%** File created from subroutine bsolv **'
      write(95,*)'%**     Assembled coefficient matrix A **'
      write(95,*)'%****************************************'
      write(95,*)'% Total dof      (neqns) =',neqns
      write(95,*)'% Total nonzeros (nnz)   =',ICSR(neqns+1)-1
      write(95,*)' neq =    ',neqns,' ;'
      write(95,*)' nnz =    ',ICSR(neqns+1)-1,' ;'

      do ii = 1, neqns
        write(96,11) aload(ii)
        ibeg = ICSR(ii)
        iend = ICSR(ii+1)-1
        do jj = ibeg, iend
          if (ii.eq.1.and.jj.eq.ibeg) then
            if (iflag.eq.1) then
              write(95,14)ii,JCSR(jj),atang(jj)
            else if (iflag.eq.2) then
              write(95,14)JCSR(jj),ii,atang(jj)
            endif
          else if (ii.eq.neqns.and.jj.eq.iend) then
            if (iflag.eq.1) then
              write(95,15)ii,JCSR(jj),atang(jj)
            else if (iflag.eq.2) then
              write(95,15)JCSR(jj),ii,atang(jj)
            endif
          else
            if (iflag.eq.1) then
              write(95,12)ii,JCSR(jj),atang(jj)
            else if (iflag.eq.2) then
              write(95,12)JCSR(jj),ii,atang(jj)
            endif
          endif

        end do
      end do
 11   format(e25.14)
 12   format(2i8,e25.14)
 14   format('mat = [',2i8,e25.14)
 15   format(2i8,e25.14,' ];')
      close(95)
      close(96)

      return
      end


c**************************************************************
      subroutine CSR2CSC (icsr,jcsr,atang,neqns)
c     Convert from CSR to CSC format.
c     Assumes a symmetric fill in.
c     Also can be used to construct a transpose.
      implicit none
      integer neqns, icsr(*),jcsr(*)
      double precision atang(*)
      integer ieq,ibeg,iend,jj,jeq,jbeg,jend,jpos,imatch,kk
      double precision atemp
c---------------------------------------------------------------
      do ieq = 1,neqns
         ibeg = icsr(ieq)
         iend = icsr(ieq+1)-1
         do jj = ibeg,iend
            jeq = jcsr(jj)
            if (jeq.gt.ieq) then
               jbeg = icsr(jeq)
               jend = icsr(jeq+1)-1
               jpos = -1
               imatch = 0
               do kk = jbeg,jend
                  if ( jcsr(kk).eq.ieq) then
                     imatch = imatch +1
                     jpos = kk
                  endif
               enddo
               if (imatch.ne.1) then
                  write(*,*)"Error in csr2csc: Matrix fillin not sym."
                  stop
               endif
               atemp = atang(jj)
               atang(jj) = atang(jpos)
               atang(jpos) = atemp
            endif
         enddo
      enddo
      
      return
      end
      
      
c**********************************************************
      subroutine timer(n) 
c     measure and print the time difference between an interval
      USE IOUNIT
      implicit none
      double precision time1, time2
      integer n,ti0,ti1,ti2,tir
      save time1, time2
      save ti0,ti1,ti2
c------------------------------------------------------------
      call cpu_time(time2)
      call system_clock(ti2,tir)
      if (n.eq.0) then
c       don't print anything, this is just for initializing the time
        ti0=ti2
      elseif (n.eq.1) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $  '    -- Primal problem done on master core in       ',
     $  dble(ti2-ti1)/dble(tir),'s (All core(s) ',
     $  time2-time1,'s)'
      elseif (n.eq.2) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $  '    -- Dual/adjoint problem done on master core in ',
     $  dble(ti2-ti1)/dble(tir),'s (All core(s) ',
     $  time2-time1,'s)'
      elseif (n.eq.3) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $  '    -- Gradient built on master core in            ',
     $  dble(ti2-ti1)/dble(tir),'s (All core(s) ',
     $  time2-time1,'s)'
      elseif (n.eq.4) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $  '    -- Initialize done on master core in           ',
     $  dble(ti2-ti1)/dble(tir),'s (All core(s) ',
     $  time2-time1,'s)'
      else
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $  '    -- Step done on master core in                 ',
     $  dble(ti2-ti1)/dble(tir),'s (All core(s) ',
     $  time2-time1,'s)'
      endif
      time1 = time2
      ti1=ti2
      
      return
      end


c**********************************************************
      subroutine timerMain(n)
c     measure and print the time difference for a step in the main
c      file and give the total time
      USE IOUNIT
      implicit none
      double precision timeMain1, timeMain2
      integer n,tim0,tim1,tim2,timr
      save timeMain1, timeMain2
      save tim0,tim1,tim2
c------------------------------------------------------------
      call cpu_time(timeMain2)
      call system_clock(tim2,timr)
      if (n.eq.0) then
c       initialization - no print
        tim0=tim2
      else 
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     ' ## Elapsed time for this step on the master core:',
     $     dble(tim2-tim1)/dble(timr),'s (All core(s):',
     $     timeMain2-timeMain1,'s)'
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     ' ## Total elapsed time on the master core:        ',
     $     dble(tim2-tim0)/dble(timr),'s (All core(s):',
     $     timeMain2,'s)'
      endif
      timeMain1 = timeMain2
      tim1=tim2
      return
      end
      

c**********************************************************
      subroutine timerIteration(n)
c     measure and print the time difference for every iteration
      USE IOUNIT
      implicit none
      double precision timeIter1, timeIter2
      integer n,tii1,tii2,tiir
      save timeIter1, timeIter2
      save tii1,tii2
c------------------------------------------------------------
      call cpu_time(timeIter2)
      call system_clock(tii2,tiir)
      if (n.eq.0) then
c       initialization - no print
      else
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '  ++ Elapsed time for this iteration on the master core:',
     $     dble(tii2-tii1)/dble(tiir),'s (All core(s):',
     $     timeIter2-timeIter1,'s)'
      endif
      timeIter1 = timeIter2
      tii1=tii2
      return
      end


c**********************************************************
      subroutine timerProfiling(n)
c     measure and print the time difference between an interval
      USE IOUNIT
      implicit none
      integer n,tip1,tip2,tipr
      double precision timeP1, timeP2
      save timeP1, timeP2
      save tip1,tip2
c------------------------------------------------------------
      call cpu_time(timeP2)
      call system_clock(tip2,tipr)
      if (n.eq.0) then
c       initialization - no print
      elseif (n.eq.1) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** LoadIncrement done in:   ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      elseif (n.eq.2) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** T. Stif. + RHS built in: ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      elseif (n.eq.3) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** residual(1) done in:     ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      elseif (n.eq.4) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** solve done in:           ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      elseif (n.eq.5) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** residual(2) done in:     ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      elseif (n.eq.6) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** update_primal done in:   ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      elseif (n.eq.7) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** ContinIncrement done in: ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      else
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '      ** step done in:            ',
     $     dble(tip2-tip1)/dble(tipr),'s (All core(s):',
     $     timeP2-timeP1,'s)'
      endif
      timeP1 = timeP2
      tip1 = tip2
      return
      end


c**********************************************************
      subroutine timerPardiso(n)
c     measure and print the time difference between an interval
      USE IOUNIT
      implicit none
      integer n,tia1,tia2,tiar
      double precision timePa1, timePa2
      save timePa1, timePa2
      save tia1,tia2
c------------------------------------------------------------
      call cpu_time(timePa2)
      call system_clock(tia2,tiar)
      if (n.eq.0) then
c       initialization - no print
      elseif (n.eq.1) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '        oo reordering done in:       ',
     $     dble(tia2-tia1)/dble(tiar),'s (All core(s):',
     $     timePa2-timePa1,'s)'
      elseif (n.eq.2) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '        oo factorization done in:    ',
     $     dble(tia2-tia1)/dble(tiar),'s (All core(s):',
     $     timePa2-timePa1,'s)'
      elseif (n.eq.3) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '        oo backsubstitution done in: ',
     $     dble(tia2-tia1)/dble(tiar),'s (All core(s):',
     $     timePa2-timePa1,'s)'
      elseif (n.eq.4) then
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '        oo memory released in:       ',
     $     dble(tia2-tia1)/dble(tiar),'s (All core(s):',
     $     timePa2-timePa1,'s)'
      else
        write(iwrit,'(A,1p,E10.4,A,1p,E10.4,A)')
     $     '        oo step done in:             ',
     $     dble(tia2-tia1)/dble(tiar),'s (All core(s):',
     $     timePa2-timePa1,'s)'
      endif
      timePa1 = timePa2
      tia1 = tia2
      return
      end

