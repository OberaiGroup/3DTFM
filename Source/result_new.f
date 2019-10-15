c***********************************************************
      subroutine result(ite,g)
c     print results 
c      USE MAINMEM! not needed here
      implicit none
      integer ite
      double precision g(*)
c-----------------------------------------------------------
      call paraview(ite,g)
      return
      end


c****************************************************************
      subroutine paraview(ite,g) 
c     Subroutine to ouput data for PARAVIEW
c     Coordinates, Connectivity and Results
      USE IOUNIT
      USE MAINMEM
      implicit none

      integer ipoin, ncelldat, nen, ntype, inode, ielem, idofn
      integer iset, imeas, ic, jelem, itemp
      integer nF,ite,nF3, surfElemCount
      integer active1, active2, active3
      character*12 itec
      double precision g(*)
      double precision xx,yy,zz
      character str1*1, str2*50, str3*1
      integer i
c----------------------------------------------------------------

      write(iwrit,*)'---------------------'
      write(iwrit,*)'Post processor format: paraview'
      write(iwrit,*)'---------------------'
      if (ite.le.9) then
        nF3=1
      elseif (ite.le.99) then
        nF3=2
      elseif (ite.le.999) then
        nF3=3
      elseif (ite.le.9999) then
        nF3=4
      elseif (ite.le.99999) then
        nF3=5
      else
        write(iwrit,*) 'result.f (paraview): too many iter. (>=1e6)'
      endif
      write(itec,'(I12)') ite
      nF=1
      do while (prefixOut(nF:nF).ne." ")
        nF=nF+1
      enddo
      prefixOut(nF:nF)='i'
      prefixOut((nF+1):(nF+nF3))=itec((12-nF3+1):12)
      prefixOut((nF+nF3+1):(nF+nF3+4))='.vtk'
      open(12,file=prefixOut(1:(nF+nF3+4)))
      if (problemType.eq.1) then! print only for the inverse problem
        prefixOut((nF+nF3+1):(nF+nF3+4))='.res'
        open(24,file=prefixOut(1:(nF+nF3+4)))
        write(24,'(I1)') nset_nodal
      endif! if (problemType.eq.1)
      prefixOut(nF:(nF+5))='      '

c     Header of *.vtk file
      write(12,'(A)') '# vtk DataFile Version 2.0'
      write(12,'(A)') 'Output of NLACE'!comments if any (must fit on one line)
      write(12,'(A)') 'ASCII'!begining of vtk format
      write(12,'(A)') 'DATASET UNSTRUCTURED_GRID'
c     Nodal coordinates
      write(12,'(A,I12,A)') 'POINTS ',npoin,' FLOAT'

      do ipoin = 1,npoin
        if (ndime.eq.3) then
          xx = coord(1,ipoin)
          yy = coord(2,ipoin)
          zz = coord(3,ipoin)
        else if (ndime.eq.2) then
          xx = coord(1,ipoin)
          yy = coord(2,ipoin)
          zz = 0.0d0
        else if (ndime.eq.1) then
          xx = coord(1,ipoin)
          yy = 0.0d0
          zz = 0.0d0
        endif
        write(12,101)xx,yy,zz
      enddo
     
c     Print connectivity
      ncelldat = 0
      jelem = 0
      do ielem = 1,nelem
         if ((lrefn(ielem,1).ne.306).and.(lrefn(ielem,1).ne.611)
     $           .and.(lrefn(ielem,1).ne.309)) then! skip surface elements
           nen = lrefn(ielem,2)
c         transform a biquadratic element into 4 bilinear elements
          if (nen.eq.9.and.ndime.eq.2) then
            jelem = jelem+4
            ncelldat = ncelldat + 20
          else
            jelem = jelem + 1
            ncelldat = ncelldat + nen +1
          endif
        endif
      enddo

      write(12,'(A,I12,I12)') 'CELLS ',jelem,ncelldat
      
      do ielem = 1,nelem

        nen = lrefn(ielem,2)
        if ((lrefn(ielem,1).ne.306).and.(lrefn(ielem,1).ne.611)
     $ .and.(lrefn(ielem,1).ne.309)) then! skip surface elements
         if (lrefn(ielem,1).ne.608) then! 608 has a different numbering (nodes' order differs)
          if (nen.eq.4.and.ndime.eq.3) then
c           if the element is a tetrahedron then reconcile the 
c           numbering between fem notes(TJRH,pg 170) and vtk (notes, pg9).
            if (npoin.le.9999) then
              write(12,'(I1,I5,I5,I5,I5)')nen,lnods(ielem,4)-1,
     $           lnods(ielem,1)-1,lnods(ielem,2)-1,lnods(ielem,3)-1
            elseif (npoin.le.99999) then
              write(12,'(I1,I6,I6,I6,I6)')nen,lnods(ielem,4)-1,
     $           lnods(ielem,1)-1,lnods(ielem,2)-1,lnods(ielem,3)-1
            elseif (npoin.le.999999) then
              write(12,'(I1,I7,I7,I7,I7)')nen,lnods(ielem,4)-1,
     $           lnods(ielem,1)-1,lnods(ielem,2)-1,lnods(ielem,3)-1
            elseif (npoin.le.999999) then
              write(12,'(I1,I8,I8,I8,I8)')nen,lnods(ielem,4)-1,
     $           lnods(ielem,1)-1,lnods(ielem,2)-1,lnods(ielem,3)-1
            else
              write(12,*)nen,lnods(ielem,4)-1,
     $           lnods(ielem,1)-1,lnods(ielem,2)-1,lnods(ielem,3)-1
            endif
          elseif (nen.eq.9.and.ndime.eq.2) then
c           if the element is biquadratic write as 4 bilinears
            write(12,*)4, 
     $           lnods(ielem,1)-1,lnods(ielem,5)-1,
     $           lnods(ielem,9)-1,lnods(ielem,8)-1
            write(12,*)4, 
     $           lnods(ielem,5)-1,lnods(ielem,2)-1,
     $           lnods(ielem,6)-1,lnods(ielem,9)-1
            write(12,*)4, 
     $           lnods(ielem,8)-1,lnods(ielem,9)-1,
     $           lnods(ielem,7)-1,lnods(ielem,4)-1
            write(12,*)4, 
     $           lnods(ielem,9)-1,lnods(ielem,6)-1,
     $           lnods(ielem,3)-1,lnods(ielem,7)-1
c          write(12,*)8,(lnods(ielem,inode)-1,inode=1,8)
          elseif (nen.eq.4.and.ndime.eq.2) then
            if (npoin.le.9999) then
              write(12,'(I1,I5,I5,I5,I5)')nen,
     $           (lnods(ielem,inode)-1,inode=1,nen)
            elseif (npoin.le.99999) then
              write(12,'(I1,I6,I6,I6,I6)')nen,
     $           (lnods(ielem,inode)-1,inode=1,nen)
            elseif (npoin.le.999999) then
              write(12,'(I1,I7,I7,I7,I7)')nen,
     $           (lnods(ielem,inode)-1,inode=1,nen)
            else
              write(12,*)nen,(lnods(ielem,inode)-1,inode=1,nen)
            endif
          else
            write(12,*)nen,(lnods(ielem,inode)-1,inode=1,nen)
          endif
         else!if (lrefn(ielem,1).ne.608) then
          if (nen.eq.4.and.ndime.eq.2) then
             write(12,*)4,lnods(ielem,1)-1,lnods(ielem,2)-1,
     $           lnods(ielem,4)-1,lnods(ielem,3)-1
          elseif (nen.eq.9.and.ndime.eq.2) then
             write(12,*)4,
     $           lnods(ielem,1)-1,lnods(ielem,2)-1,
     $           lnods(ielem,5)-1,lnods(ielem,4)-1
             write(12,*)4,
     $           lnods(ielem,2)-1,lnods(ielem,3)-1,
     $           lnods(ielem,6)-1,lnods(ielem,5)-1
             write(12,*)4,
     $           lnods(ielem,4)-1,lnods(ielem,5)-1,
     $           lnods(ielem,8)-1,lnods(ielem,7)-1
             write(12,*)4,
     $           lnods(ielem,5)-1,lnods(ielem,6)-1,
     $           lnods(ielem,9)-1,lnods(ielem,8)-1
          else
            print*,"result.f: paraview: case not implemented ...",
     $           " exiting"
            stop
          endif
         endif!if (lrefn(ielem,1).ne.608) then
        endif!if (lrefn(ielem,1).ne.306) then! skip the surface elements
      enddo
      
      surfElemCount=0
      write(12,'(A,I12)') 'CELL_TYPES ',jelem
      do ielem = 1,nelem
        nen = lrefn(ielem,2)
        if ((lrefn(ielem,1).eq.306).or.(lrefn(ielem,1).eq.611)
     $            .or.(lrefn(ielem,1).eq.309)) then! count the surface elements
          surfElemCount=surfElemCount+1
        else
          if (nen.eq.3.and.ndime.eq.2) then
            ntype = 5
            write(12,'(I1)')ntype
          elseif (nen.eq.4.and.ndime.eq.2) then 
            ntype = 9
            write(12,'(I1)')ntype
          elseif (nen.eq.4.and.ndime.eq.3) then
            ntype = 10
            write(12,'(I2)')ntype
          elseif (nen.eq.8.and.ndime.eq.3) then
            ntype = 12
            write(12,'(I2)')ntype
          elseif (nen.eq.9.and.ndime.eq.2) then
            ntype = 9
            write(12,'(I1)')ntype
            write(12,'(I1)')ntype
            write(12,'(I1)')ntype
            write(12,'(I1)')ntype
          endif
        endif
      enddo

c     Print results
c     WARNING: WORKS ONLY FOR THE SAME # OF DOFS/NODE
      write(12,'(A,I12)') 'POINT_DATA ',npoin

c     for each set of measured data:
      ic = 0
      do imeas = 1,nmeas
        ic = ic +1 
        if ((problemType.eq.1).or.
     &      ((problemType.eq.2).and.(lrefn(1,1).eq.31)).or.
     &      ((problemType.eq.2).and.(lrefn(1,1).eq.32))) then! write for the inverse problem or the filter
c         write the measured field (given in input file)
          do idofn = 1,mdofn
            write(unit=str1,fmt='(I1)') ic
            write(unit=str3,fmt='(I1)') idofn
            str2 ='SCALARS measuredDisp'//str1//str3//' float  1'
            write(12,*)str2
            write(12,106)
            do ipoin = 1,npoin
              write(12,112)meas(imeas,ipoin,idofn)
            enddo
          enddo
        endif! if (problemType.eq.1)

c       write the primal field 
        do idofn = 1,mdofn
          write(unit=str1,fmt='(I1)') ic
          write(unit=str3,fmt='(I1)') idofn
          str2 ='SCALARS predictedDisp'//str1//str3//' float  1'
          write(12,*)str2
          write(12,106)
          do ipoin = 1,npoin
            write(12,112)total_primal(imeas,ipoin,idofn)
          enddo
        enddo
        
        if (problemType.eq.1) then! write only for the inverse problem
c         write the dual 
          do idofn = 1,mdofn
            write(unit=str1,fmt='(I1)') ic
            write(unit=str3,fmt='(I1)') idofn
            str2 ='SCALARS dual'//str1//str3//' float  1'
            write(12,*)str2
            write(12,106)
            do ipoin = 1,npoin
              write(12,112)total_dual(imeas,ipoin,idofn)
            enddo
          enddo
        endif! if (problemType.eq.1)
      enddo

      if ((problemType.eq.1).or.
     &    ((problemType.eq.2).and.(lrefn(1,1).ne.31).and.
     &                            (lrefn(1,1).ne.32))) then
c       write the nodal data
        ic = 0
        do iset = 1,nset_nodal
          ic = ic + 1
          write(unit=str1,fmt='(I1)') ic
          if (ic.eq.1) then! NLP=non-linear parameter
            str2 ='SCALARS NLP float  1'
          elseif (ic.eq.2) then! SM=shear modulus
            str2 ='SCALARS SM float  1'
          else
            str2 ='SCALARS property'//str1//' float  1'
          endif
          write(12,*)str2
          write(12,106)
          do ipoin = 1,npoin
            write(12,112)adata_nodal(iset,ipoin)
          enddo
        enddo
      endif

      if (problemType.eq.1) then! print only for an inverse problem
        if (nset_nodal.eq.1) then
          do ipoin=1,npoin
            if (ipoin_nod_grad(1,ipoin).eq.0) then
              active1=0
            else
              active1=1
            endif
            write(24,'(I1,1p,E11.4)') active1,adata_nodal(1,ipoin)
          enddo
        elseif (nset_nodal.eq.2) then
          do ipoin=1,npoin
            if (ipoin_nod_grad(1,ipoin).eq.0) then
              active1=0
            else
              active1=1
            endif
            if (ipoin_nod_grad(2,ipoin).eq.0) then
              active2=0
            else
              active2=1
            endif
            write(24,'(I1,1p,E11.4,I2,1p,E11.4)') 
     +                  active1,adata_nodal(1,ipoin),
     +                  active2,adata_nodal(2,ipoin)
          enddo
        elseif (nset_nodal.eq.3) then
          do ipoin=1,npoin
            if (ipoin_nod_grad(1,ipoin).eq.0) then
              active1=0
            else
              active1=1
            endif
            if (ipoin_nod_grad(2,ipoin).eq.0) then
              active2=0
            else
              active2=1
            endif
            if (ipoin_nod_grad(3,ipoin).eq.0) then
              active3=0
            else
              active3=1
            endif
            write(24,'(I1,1p,E11.4,I2,1p,E11.4,I2,1p,E11.4)')
     +                  active1,adata_nodal(1,ipoin),
     +                  active2,adata_nodal(2,ipoin),
     +                  active3,adata_nodal(3,ipoin)
          enddo
        else
          write(24,*) 'result.f: nset_nodal is too large!!!'
        endif! if (nset_nodal.eq.X)
      endif! if (problemType.eq.1)

      if (problemType.eq.1) then! write only for the inverse problem
c       write the gradient for nodal data 
        ic = 0
        do iset = 1,nset_nodal
          ic = ic +1
          write(unit=str1,fmt='(I1)') ic
          str2 ='SCALARS gradient'//str1//' float  1'
          write(12,*)str2
          write(12,106)
          if (g(1).eq.-1234567.0d0) then! this is the dummy value cf lmain.f
           do ipoin = 1,npoin
             write(12,112) 0.0d0! make a default at 0.0d0
           enddo
          else
           do ipoin = 1,npoin
            itemp = ipoin_nod_grad(iset,ipoin)
            if (itemp.ne.0) then
              if (g(itemp).gt.1.0d12) then
                write(12,112) 1.0d12!grad_nodal(iset,ipoin) - grad_nodal has been removed
              elseif (g(itemp).lt.(-1.0d12)) then
                write(12,112) -1.0d12!grad_nodal(iset,ipoin) - grad_nodal has been removed
              elseif((g(itemp).le.1.0d12).and.
     &               (g(itemp).ge.(-1.0d12)))then
                write(12,112) g(itemp)!grad_nodal(iset,ipoin) - grad_nodal has been removed
              else
                print*,'warning :g(',itemp,')',g(itemp)
                write(12,112) g(itemp)!grad_nodal(iset,ipoin) - grad_nodal has been removed
              endif
            else
              write(12,112) 0.0d0
            endif
           enddo
          endif
        enddo

        write(12,'(A,I12)') 'CELL_DATA',jelem!nelem-surfElemCount - JFD
        write(12,'(A)') 'SCALARS elemDataMatch  float  1'
        write(12,106)
        do ielem = 1,nelem
          if ((lrefn(ielem,1).ne.306).and.(lrefn(ielem,1).ne.611)
     $    .and.(lrefn(ielem,1).ne.309)) then! skip surface elements                   
              write(12,112) elemDataMatch(ielem)
            if ((lrefn(ielem,1).eq.608).or.(lrefn(ielem,1).eq.606)) then! write x more times the value
              nen = lrefn(ielem,2)
              if (nen.eq.4) then
                nen=2-1
              elseif (nen.eq.9) then
                nen=4-1
              elseif (nen.eq.16) then
                nen=9-1
              elseif (nen.eq.25) then
                nen=16-1
              else
                print*,"results.f: paraview: case not implemented",
     $                 " ... exiting"
                stop
              endif
              if (nen.gt.1) then
                do i=1,nen
                  write(12,112) elemDataMatch(ielem)
                enddo
              endif
            endif
          endif
        enddo
        write(12,'(A)') 'SCALARS elemRegul  float  1'
        write(12,106)
        do ielem = 1,nelem
          if ((lrefn(ielem,1).ne.306) .and.
     $       (lrefn(ielem,1).ne.611).and.
     $        (lrefn(ielem,1).ne.309))then! skip surface elements
            write(12,112) elemRegul(ielem)
            if ((lrefn(ielem,1).eq.608).or.(lrefn(ielem,1).eq.606)) then! write x more times the value
              nen = lrefn(ielem,2)
              if (nen.eq.4) then
                nen=2-1
              elseif (nen.eq.9) then
                nen=4-1
              elseif (nen.eq.16) then
                nen=9-1
              elseif (nen.eq.25) then
                nen=16-1
              else
                print*,"results.f: paraview: case not implemented",
     $                 " ... exiting"
                stop
              endif
              if (nen.gt.1) then
                do i=1,nen
                  write(12,112) elemRegul(ielem)
                enddo
              endif
            endif
          endif
        enddo
        write(12,'(A)') 'SCALARS l2grad2  float  1'
        write(12,106)
        do ielem = 1,nelem
          if ((lrefn(ielem,1).ne.306).and.(lrefn(ielem,1).ne.611)
     $        .and.(lrefn(ielem,1).ne.309) ) then! skip surface elements
            write(12,112) l2grad2(ielem)
            if ((lrefn(ielem,1).eq.608).or.(lrefn(ielem,1).eq.606)) then! write x more times the value
              nen = lrefn(ielem,2)
              if (nen.eq.4) then
                nen=2-1
              elseif (nen.eq.9) then
                nen=4-1
              elseif (nen.eq.16) then
                nen=9-1
              elseif (nen.eq.25) then
                nen=16-1
              else
                print*,"results.f: paraview: case not implemented",
     $                 " ... exiting"
                stop
              endif
              if (nen.gt.1) then
                do i=1,nen
                  write(12,112) l2grad2(ielem)
                enddo
              endif
            endif
          endif
        enddo
      endif! if (problemType.eq.1)
            
            
 101  format(1p,e13.6,1p,e14.6,1p,e14.6)
 106  format('LOOKUP_TABLE default')
 112  format(1p,e12.5)

      close(12,status='keep')
      if (problemType.eq.1) close(24,status='keep')

      return
      end

