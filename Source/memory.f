c     ---------------------------------------------------------
      subroutine MEMORY_ALLOCA
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------

      mevab = mnode*mdofn

      ALLOCATE (ldofn(npoin), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: ldofn'
         stop
      endif
      ldofn(:) = 0

      ALLOCATE (elemvec_ndofn(mnode), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: elemvec_ndofn'
         stop
      endif
      elemvec_ndofn(:) = 0

      ALLOCATE (lnods(nelem,mnode), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: ldofn'
         stop
      endif
      lnods(:,:) = 0

      ALLOCATE (elemDataMatch(nelem), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: elemDataMatch'
         stop
      endif
      elemDataMatch(:) = 0.0d0

      ALLOCATE (elemRegul(nelem), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: elemRegul'
         stop
      endif
      elemRegul(:) = 0.0d0

      ALLOCATE (l2grad1(nelem), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: l2grad2'
         stop
      endif
      l2grad1(:) = 0.0d0

      ALLOCATE (l2grad2(nelem), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: l2grad2'
         stop
      endif
      l2grad2(:) = 0.0d0

      ALLOCATE (lrefn(nelem,3), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: lrefn'
         stop
      endif
      lrefn(:,:) = 0

      ALLOCATE (idnum(npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: idnum'
         stop
      endif
      idnum(:,:) = 0

      ALLOCATE (tvect(mdofn*npoin), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: tvect'
         stop
      endif
      tvect(:) = 0.0d0

      ALLOCATE (itvect(mdofn*npoin), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: itvect'
         stop
      endif
      itvect(:) = 0

      ALLOCATE (props(nmat,mprops), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: props'
         stop
      endif
      props(:,:) = 0.0d0

      ALLOCATE (coord(ndime,npoin), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: coord'
         stop
      endif
      coord(:,:) = 0.0d0

      ALLOCATE (coord_copy(ndime,npoin), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: coord_copy'
         stop
      endif
      coord_copy(:,:) = 0.0d0

      ALLOCATE (uelem_dual(mdofn,mnode), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: uelem_dual'
         stop
      endif
      uelem_dual(:,:) = 0.0d0

      ALLOCATE (primal(npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: primal'
         stop
      endif
      primal(:,:) = 0.0d0

      ALLOCATE (dual(npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: dual'
         stop
      endif
      dual(:,:)   = 0.0d0


      ALLOCATE (meas(nmeas,npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: meas'
         stop
      endif
      meas(:,:,:) = 0.0d0

      ALLOCATE (total_dual(nmeas,npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: total_dual'
         stop
      endif
      total_dual(:,:,:) = 0.0d0

      ALLOCATE (total_primal(nmeas,npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: total_primal'
         stop
      endif
      total_primal(:,:,:) = 0.0d0

      ALLOCATE (uelem_diff(mdofn,mnode), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: uelem_diff'
         stop
      endif
      uelem_diff(:,:) = 0.0d0

c     Defining uelem_meas is not compatible with openmp; if needed make it local
c      ALLOCATE (uelem_meas(mdofn,mnode), STAT=status)
c      if (status.ne.0) then
c         write(*,*)'memory failiure: uelem_meas'
c         stop
c      endif
c      uelem_meas(:,:) = 0.0d0

      ALLOCATE (fac_meas(nmeas,npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: fac_meas'
         stop
      endif
      fac_meas(:,:,:)=0.0d0

      ALLOCATE (fac_meas_nd(nmeas,npoin,mdofn*(mdofn-1)/2), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: fac_meas_nd'
         stop
      endif
      fac_meas_nd(:,:,:)=0.0d0

      ALLOCATE (bodyforce(npoin,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: bodyforce'
         stop
      endif
      bodyforce(:,:)=0.0d0

      ALLOCATE (icode_bc(nmeas,mpoinbc,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: icode_bc'
         stop
      endif
      icode_bc(:,:,:) = 0

      ALLOCATE (value_bc(nmeas,mpoinbc,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: value_bc'
         stop
      endif
      value_bc(:,:,:) = 0.0d0

      ALLOCATE (icode_ws(nmeas,mpoinbc,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: icode_bc'
         stop
      endif
      icode_bc(:,:,:) = 0

      ALLOCATE (value_ws(nmeas,mpoinbc,mdofn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: value_bc'
         stop
      endif
      value_bc(:,:,:) = 0.0d0


c      ALLOCATE (value_bc2(nmeas,mpoinbc,mdofn), STAT=status)
c      if (status.ne.0) then
c         write(*,*)'memory failiure: value_bc2'
c         stop
c      endif
c      value_bc2(:,:,:) = 0.0d0

      ALLOCATE (ipoin_bc(nmeas,mpoinbc), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: ipoin_bc'
         stop
      endif
      ipoin_bc(:,:) = 0

      ALLOCATE (ipoin_ws(nmeas,mpoinbc), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: ipoin_bc'
         stop
      endif
      ipoin_bc(:,:) = 0


      return
      end

c     ---------------------------------------------------------
      subroutine MEMORY_ALLOCB
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------


      ALLOCATE (nod_el_coun(neqns+1), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: nod_el_coun'
         stop
      endif
      nod_el_coun(:) = 0

      ALLOCATE (jdiag(neqns+1), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: neqns'
         stop
      endif
      jdiag(:) = 0

      ALLOCATE (jdiag_temp(neqns+1), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: neqns'
         stop
      endif
      jdiag_temp(:) = 0

      ALLOCATE (aload(neqns), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: aload'
         stop
      endif
      aload(:) = 0.0d0

      return
      end

c     ---------------------------------------------------------
      subroutine MEMORY_ALLOCC
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------

      ALLOCATE (nod_el_conn(nnzconn), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: nnzconn'
         stop
      endif
      nod_el_conn(:) = 0

      return
      end
c     ---------------------------------------------------------
      subroutine MEMORY_ALLOCD
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------
      ALLOCATE (jcsr(nnz), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: jcsr: status = ',status
         stop
      endif
      jcsr(:) = 0

      ALLOCATE (jcsr_temp(nnz), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: jcsr: status = ',status
         stop
      endif
      jcsr_temp(:) = 0

      ALLOCATE (atang(nnz), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: atang',status
         stop
      endif
      atang(:) = 0.0d0

      ALLOCATE (atang_temp(nnz), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: atang_temp',status
         stop
      endif
      atang_temp(:) = 0.0d0

      return
      end

c     ---------------------------------------------------------
      subroutine MEMORY_ALLOCE
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------
      ALLOCATE (adata_nodal(nset_nodal,npoin), STAT=status)
      if (status.ne.0) then
        write(*,*)'memory failiure: adata_nodal: status = ',status
        stop
      endif
      adata_nodal(:,:) = 0.0d0

      ALLOCATE (adata_nodal2(nset_nodal,npoin), STAT=status)
      if (status.ne.0) then
        write(*,*)'memory failiure: adata_nodal2: status = ',status
        stop
      endif
      adata_nodal2(:,:) = 0.0d0

      ALLOCATE (delta_adata_nodal(nset_nodal,npoin), STAT=status)
      if (status.ne.0) then
        write(*,*)'memory failiure: delta_adata_nodal: status = ',status
        stop
      endif
      delta_adata_nodal(:,:) = 0.0d0

c      grad_nodal is not used anymore
c       ALLOCATE (grad_nodal(nset_nodal,npoin), STAT=status)
c       if (status.ne.0) then
c          write(*,*)'memory failiure: grad_nodal: status = ',status
c          stop
c       endif
c       grad_nodal(:,:) = 0.0d0

      ALLOCATE (ielem_nod_grad(nset_nodal,mnode), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: ielem_nod_grad: status = ',status
         stop
      endif
      ielem_nod_grad(:,:)=0

      ALLOCATE (egrad(nset_nodal,mnode), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: egrad: status = ',status
         stop
      endif
      egrad(:,:) = 0.0d0

      ALLOCATE (ipoin_nod_grad(nset_nodal,npoin), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: ipoin_nod_grad: status = ',status
         stop
      endif
      ipoin_nod_grad(:,:) = 0
      return
      end


c     ---------------------------------------------------------
      subroutine MEMORY_ALLOCF
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------
      ALLOCATE (adata_elemental(nset_elemental,nelem), STAT=status)
      if (status.ne.0) then
         write(*,*)'memory failiure: adata_elemental: status = ',status
         stop
      endif
      adata_elemental(:,:) = 0.0d0

      return
      end



c     ---------------------------------------------------------
      subroutine MEMORY_DEALLOC
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------

      DEALLOCATE ( ldofn, STAT = status)      
      DEALLOCATE ( lnods, STAT = status)
      DEALLOCATE ( lrefn, STAT = status)
      DEALLOCATE ( idnum, STAT = status)
      DEALLOCATE ( props, STAT = status)
      DEALLOCATE ( tvect, STAT = status)
      DEALLOCATE ( itvect, STAT = status)
      DEALLOCATE ( coord, STAT = status)
      DEALLOCATE ( uelem_dual, STAT = status)
      DEALLOCATE ( primal, STAT = status)
      DEALLOCATE ( dual, STAT = status)
      DEALLOCATE ( meas, STAT = status)
      DEALLOCATE ( total_primal, STAT = status)
      DEALLOCATE ( total_dual, STAT = status)
      DEALLOCATE ( uelem_diff, STAT = status)
c make it local      DEALLOCATE ( uelem_meas, STAT = status)
      DEALLOCATE ( fac_meas, STAT = status)
      DEALLOCATE ( fac_meas_nd, STAT = status)
      DEALLOCATE ( bodyforce, STAT = status)
      DEALLOCATE ( icode_bc, STAT = status)
      DEALLOCATE ( value_bc, STAT = status)
      DEALLOCATE ( ipoin_bc, STAT = status)
      DEALLOCATE ( icode_ws, STAT = status)
      DEALLOCATE ( value_ws, STAT = status)
      DEALLOCATE ( ipoin_ws, STAT = status)
      DEALLOCATE ( nod_el_coun, STAT = status)
      DEALLOCATE ( jdiag, STAT = status)
      DEALLOCATE ( jdiag_temp, STAT = status)
      DEALLOCATE ( aload, STAT = status)
      DEALLOCATE ( nod_el_conn, STAT = status)
      DEALLOCATE ( jcsr, STAT = status)
      DEALLOCATE ( jcsr_temp, STAT = status)
      DEALLOCATE ( atang, STAT = status)
      DEALLOCATE ( atang_temp, STAT = status)
      DEALLOCATE ( adata_nodal, STAT = status)
      DEALLOCATE ( adata_nodal2, STAT = status)
c JFD      DEALLOCATE ( grad_nodal, STAT = status)
      DEALLOCATE ( ielem_nod_grad, STAT = status)
      DEALLOCATE ( egrad, STAT = status)
      DEALLOCATE ( ipoin_nod_grad, STAT = status)
      DEALLOCATE ( adata_elemental, STAT = status)
c     added new
      DEALLOCATE ( elemDataMatch, STAT = status)
c     end new

      return
      end





c     ---------------------------------------------------------
      subroutine MEMORY_DEALLOCBCD
c     ---------------------------------------------------------
      USE MAINMEM
      implicit none
      integer status
c     ---------------------------------------------------------

      DEALLOCATE ( nod_el_coun, STAT = status)
      DEALLOCATE ( jdiag, STAT = status)
      DEALLOCATE ( jdiag_temp, STAT = status)
      DEALLOCATE ( aload, STAT = status)
      DEALLOCATE ( nod_el_conn, STAT = status)
      DEALLOCATE ( jcsr, STAT = status)
      DEALLOCATE ( jcsr_temp, STAT = status)
      DEALLOCATE ( atang, STAT = status)
      DEALLOCATE ( atang_temp, STAT = status)

      return
      end








