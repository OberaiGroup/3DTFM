c******************************************************
      subroutine inicsr
c     Evaluate row and column pointer arrays, icsr and jcsr, 
c     for Compressed Sparse Row (CSR) storage of tangent
c     matrix

c     At the end of this routine the following arrays have been 
c     created and filled in:
c     icsr, jcsr, atang (not filled in), idnum
c     nod_el_coun, nod_el_conn
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ipoin,idofn,j,mj
      integer status
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
c------------------------------------------------------

c     set up equation numbers
      neqns = 0
      do ipoin = 1,npoin
        do idofn = 1,ldofn(ipoin)
          j = idnum(ipoin,idofn)
          if (j.eq.0) then
            neqns = neqns + 1
            idnum(ipoin,idofn) = neqns
          else if (j.lt.0) then
            mj = -j
            if (mj.lt.ipoin) then
              idnum(ipoin,idofn) =mj
            else if (mj.gt.ipoin) then
              neqns = neqns + 1
              idnum(ipoin,idofn) = neqns
              idnum(mj,idofn) = -ipoin
            endif
          else if (j.eq.1) then
            idnum(ipoin,idofn) = 0
          end if
        enddo
      enddo

      ALLOCATE (itemp(neqns+1), STAT=status)
      if (status.ne.0) then
         write(*,*)'Source/inicsr.f: memory failiure: neqns'
         stop
      endif
      
      call MEMORY_ALLOCB

c      write(iwrit,*)'  Source/inicsr.f: NEQNS ',neqns
      
      call nod_elem_conn (itemp)
      
      call csr1 (itemp)

      DEALLOCATE (itemp, STAT = status)
      
      return
      end


c******************************************************
      subroutine nod_elem_conn (itemp)
c     create the node to element connectivity array
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer itemp (*)
      integer itemp1,ii,ielem,inode,iidof,idofn,kk,istart,ieq
      integer nnot
c------------------------------------------------------

      nnzconn = 0
          
      do ii = 1,neqns
        nod_el_coun(ii) = 0
      enddo
       
      do ielem = 1,nelem
        nnot = lrefn(ielem,2)! number of nodes in the element
        do inode = 1,nnot
          ii = lnods(ielem,inode)
          
          iidof = ldofn(ii)
          do idofn = 1,iidof
            kk = idnum(ii,idofn)
            if (kk.ne.0) then
              nod_el_coun(kk) = nod_el_coun(kk) + 1
              nnzconn = nnzconn + 1
            endif
          enddo
        enddo
      enddo
       
c      write(iwrit,*)'  Source/inicsr.f: NNZCONN ',nnzconn

      call MEMORY_ALLOCC

      istart = 1
      do ieq = 1,neqns
        itemp1 = nod_el_coun(ieq)
        nod_el_coun(ieq) = istart
        itemp(ieq) = nod_el_coun(ieq)
        istart = istart + itemp1
      enddo
      nod_el_coun(neqns+1) = istart
      itemp(neqns+1) = nod_el_coun(neqns+1)
          


      do ii = 1,nnzconn
        nod_el_conn(ii) = 0
      enddo
          
       
      do ielem = 1,nelem
        nnot = lrefn(ielem,2)! number of nodes in the element
        do inode = 1,nnot
          ii = lnods(ielem,inode)
          
          iidof = ldofn(ii)
          do idofn = 1,iidof
            kk = idnum(ii,idofn)
            if (kk.ne.0) then
              nod_el_conn(itemp(kk)) = ielem
              itemp(kk) = itemp(kk)+1
            endif
          enddo
        enddo
      enddo

      return
      end


c******************************************************
      subroutine csr1 (itemp)
c     create icsr (called jdiag) and jcsr using the node to element
c     connectivty array nod_el_coun and nod_el_conn
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer itemp(*)
      integer ii,ieq,ic,iidof,idofn,kk,istart,ielem,inode,itemp1,ibeg
      integer iend,jj
      integer nnot
c------------------------------------------------------
      nnz = 0
      do ii = 1,neqns+1
        itemp(ii) = 0
      enddo
 
      do ieq = 1,neqns
        do ic = nod_el_coun(ieq),nod_el_coun(ieq+1)-1
          ielem = nod_el_conn(ic)
          nnot = lrefn(ielem,2)! number of nodes in the element
          do inode = 1,nnot
            ii = lnods(ielem,inode)

            iidof = ldofn(ii)
            do idofn = 1,iidof
              kk = idnum(ii,idofn)
              if ((kk.ne.0).and.(itemp(kk).ne.ieq)) then
                if ((.not.buildSymmetricMatrices).or.
     &              (buildSymmetricMatrices.and.(ieq.le.kk))) then
c     &              (buildSymmetricMatrices.and.(ieq.ge.kk))) then
                  nnz = nnz + 1
                  itemp(kk) = ieq
                endif
              endif
            enddo
          enddo
        enddo
      enddo

      call MEMORY_ALLOCD

      nnz = 0
      do ii = 1,neqns+1
        jdiag(ii) = 0
        itemp(ii) = 0
      enddo
 
      do ieq = 1,neqns
        do ic = nod_el_coun(ieq),nod_el_coun(ieq+1)-1
          ielem = nod_el_conn(ic)
          nnot = lrefn(ielem,2)! number of nodes in the element
          do inode = 1,nnot
            ii = lnods(ielem,inode)

            iidof = ldofn(ii)
            do idofn = 1,iidof
              kk = idnum(ii,idofn)
              if ((kk.ne.0).and.(itemp(kk).ne.ieq)) then
                if ((.not.buildSymmetricMatrices).or.
     &              (buildSymmetricMatrices.and.(ieq.le.kk))) then
c     &              (buildSymmetricMatrices.and.(ieq.ge.kk))) then
                  nnz = nnz + 1
                  jcsr(nnz) = kk
                  jdiag(ieq) = jdiag(ieq)+1
                  itemp(kk) = ieq
                endif
              endif
            enddo
          enddo
        enddo
      enddo

      istart = 1
      do ieq = 1,neqns
        itemp1 = jdiag(ieq)
        jdiag(ieq) = istart
        istart = istart + itemp1
      enddo
      jdiag(neqns+1) = istart
       
      do ieq = 1,neqns
        ibeg = jdiag(ieq)
        iend = jdiag(ieq+1)-1
        do jj = ibeg,iend
          itemp(jj-ibeg+1) = jcsr(jj)
        enddo
        kk =  iend-ibeg+1
        if (kk.gt.1) then
          call sort_shell(itemp, kk)
        endif
        do jj = ibeg,iend
          jcsr(jj) = itemp(jj-ibeg+1) 
        enddo
      enddo

      return
      end


c******************************************************
      SUBROUTINE sort_shell(arr,n)
c     bubble sort in ascending order (used for jcsr only - 2010-03-05)
c     modified from numerical recipies
      IMPLICIT NONE
      INTEGER arr(*)
      INTEGER  n
c     local variables
      INTEGER :: i,j,inc
      INTEGER :: v
c------------------------------------------------------
      inc=1
      do
        inc=3*inc+1
        if (inc > n) exit
      enddo

      do
        inc=inc/3
        do i=inc+1,n
          v=arr(i)
          j=i
          do
            if (arr(j-inc) <= v) exit
            arr(j)=arr(j-inc)
            j=j-inc
            if (j <= inc) exit
          end do
          arr(j)=v
        end do
        if (inc <= 1) exit
      end do
      END SUBROUTINE sort_shell


c******************************************************
      subroutine reduceSymm
c     reduces jdiag, jcsr and atang for symmetric matrices
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer countEntry, i, j, pos, posSt, posEn
      integer tmpi(neqns+1),tmpj(nnz)
c      double precision tmpv(nnz)
c------------------------------------------------------
      print*,"neqns=",neqns
      print*,jdiag
      tmpi(:)=0
      tmpj(:)=0
c      tmpv(:)=0.0d0
      countEntry=0
      tmpi(1)=1
      do i=1,neqns
        posSt=jdiag(i)
        posEn=jdiag(i+1)-1
        do pos=posSt,posEn
          j=jcsr(pos)
          if (j.ge.i) then
            countEntry=countEntry+1
c            tmpv(countEntry)=atang(j)
            tmpj(countEntry)=j
          endif
        enddo
        tmpi(i+1)=countEntry+1
      enddo
      do i=1,neqns
        jdiag(i+1)=tmpi(i+1)
      enddo
      do i=1,nnz
        jcsr(i)=tmpj(i)
c        atang(i)=tmpv(i)
      enddo
      print*,jdiag
      stop
      return
      end
