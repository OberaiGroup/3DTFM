c**********************************************************
      subroutine elem315 (ielem,itask,pelem,nnode,estif,eforc,elnods,    !!DS: elem315
     $     xelem,elemdata_nodal,uelem,uelem_meas) 
c      Incompressible Blatz material, modified from elem305 by Dawei Song. 
c      This aims to solve a forward problem formulated in the current configuration. 
c      3D, linear tetrahedra, finite elasticity, incompressible.
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem, s, ii
      parameter (mnode_elem = 4,mquad_elem=30)
      integer ielem, itask,k,ireg,iset
      integer l, ninte, iinte, inode, jnode,knode
      integer ievab, jevab, jdofn, idofn, i, j, q, r, t
      integer tmp6, tmp7
      double precision xjaco, wtjac, temp, temp_der
      double precision deno, powe, alpha(2), beta(2), Treg(2)! variables for the regularization
      double precision h, tmp1, tmp2, tmp3, tmp4, tmp5, tauMult ! variables for the stabilization
      double precision gamm, mu, Cdet, Inv1, Inv2, K1, K2(3,3), Fdet
      double precision shap(4,mnode_elem), pres(4)
      double precision sg(mquad_elem),tg(mquad_elem),zg(mquad_elem)
      double precision SecPK_grad(3,3), wg(mquad_elem)
      double precision Fdef(3,3), Ctens(3,3), dWdC(3,3), dJdC(3,3)
      double precision Btens(3,3),Sigma(3,3), BFtens(3,3)  !!DS: Define the left Cauchy-Green tensor, the Cauchy (true) stress.
      double precision Cinv(3,3), SecPK(3,3), ident(3,3), Finv(3,3)
	double precision Finv_det    !!DS: Define the determinant of the inverse of the deformation gradient.
      double precision Ctang(3,3,3,3),d2JdC(3,3,3,3)
c JFD for debug      double precision Ctangt(3,3,3,3),d2JdCt(3,3,3,3)
      double precision Igeo(9,9),Lmat(3,3,3,3), Dgeomat(9,9)
      double precision Tp(4), FS(3,3), bb(mnode_elem,9,3)
      double precision dWdC_grad(3,3),udiff,sigma_grad(3,3) !!DS: Change udiff(3) to a scalar udiff
      double precision temp_dual(4*mnode_elem),prop_grad(2,4) 
      double precision var_shape_prop(2,4)
      double precision dCinvdC(3,3,3,3)
      double precision dCinvdCt
c     Localized information
      integer nnode
      double precision pelem(*)
      double precision estif(mevab,mevab)    !!DS: Is mevab the total dof for an element?
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      integer ndim ! = ndime (local copy)
c-------------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask  
      
c     --------------------------------------------------------
c     Initialize the elemental parameters
c     --------------------------------------------------------
 1    elemvec_ndofn(1:nnode) = 4    !!DS: Incompressible
      buildSymmetricMatrices=.false.
      return


c     --------------------------------------------------------
c     Read and store the material properties
c     --------------------------------------------------------
c     l          : no. of gauss pts/direction
c     ireg       : regularization type (0/1/2/3/4:none/H1/TVD/newTV/power)
c     alpha(1)   : regularization weight for the non-linear parameter
c     alpha(2)   : regularization weight for the shear modulus
c     beta(1)    : extra parameter for TVD - JFD ask Sevan for a better name
c     beta(2)    : extra parameter for TVD - JFD ask Sevan for a better name
c     Treg(1)    : second extra parameter
c     Treg(2)    : second extra parameter
c     tauMult    : stabilization factor
c     s          : stabilization type (0/1:off/on) - JFD use stype instead of s?
 2    read(iread,*) l,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +     Treg(1),Treg(2),tauMult,s
      write(iwrit,205) l,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +     Treg(1),Treg(2),tauMult,s
 205  format(/' finite elasticity (3D,4DOF) '/      !!DS: each node has 4DOF for incompressible materials
     +       ' gauss pts/dir .........................',i12,/
     +       ' reg type(0/1/2/4:none/H1/TVD/log, etc).',i12,/
     +       ' regularization parameter 1 ............',1p,e16.4,/
     +       ' regularization parameter 2 ............',1p,e16.4,/
     +       ' extra parameter (TVD) .................',1p,e16.4,/! modify names here JFD
     +       ' extra parameter (TVD) .................',1p,e16.4,/! idem
     +       ' second extra parameter (TVD) ..........',1p,e16.4,/! modify names here JFD
     +       ' second extra parameter (TVD) ..........',1p,e16.4,/! idem
     +       ' stabilization factor ..................',1p,e16.4,/
     +       ' stabilization type (0/1:off/on) .......',i12)

c
      if (alpha(1).lt.0.0d0) then
        print*,'elem315.f: alpha(1) is less than 0.0d0: exiting'
        stop
      elseif (alpha(2).lt.0.0d0) then
        print*,'elem315.f: alpha(2) is less than 0.0d0: exiting'
        stop
      elseif (beta(1).le.0.0d0) then
        print*,'elem315.f: beta(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (beta(2).le.0.0d0) then
        print*,'elem315.f: beta(2) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(1).le.0.0d0) then
        print*,'elem315.f: Treg(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(2).le.0.0d0) then
        print*,'elem315.f: Treg(2) is less or equal to 0.0d0: exiting'
        stop
      endif
     
c     note: ireg is checked in itask.eq.7 for keeping things simple
c     note: s is checked in itask.eq.3 for keeping things simple
      pelem(1) = dble(l)
      pelem(2) = dble(ireg)
      pelem(3) = alpha(1)
      pelem(4) = alpha(2)
      pelem(5) = beta(1)
      pelem(7) = Treg(1)
      pelem(6) = beta(2)
      pelem(8) = Treg(2)
      pelem(9) = tauMult
      pelem(10) = dble(s)
      return


c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c      and the elemental RHS/residual (eforc)
c     --------------------------------------------------------
 3    l      = int(pelem(1))
      tauMult= pelem(9)
      s      = int(pelem(10))
      ndim = ndime

c     determine the characteristic length h of the element
c       (for tetrahedra: h is the length of the longest edge)  !!DS: Need this for stabilization terms
      tmp1=(xelem(1,1)-xelem(1,2))**2.0d0
      tmp2=(xelem(1,1)-xelem(1,3))**2.0d0
      tmp3=(xelem(1,1)-xelem(1,4))**2.0d0
      tmp4=(xelem(1,2)-xelem(1,3))**2.0d0
      tmp5=(xelem(1,2)-xelem(1,4))**2.0d0
      h=(xelem(1,3)-xelem(1,4))**2.0d0


      do i=2,ndim ! 1, ndim
        tmp1 = tmp1 + (xelem(i,1)-xelem(i,2))**2.0d0
        tmp2 = tmp2 + (xelem(i,1)-xelem(i,3))**2.0d0
        tmp3 = tmp3 + (xelem(i,1)-xelem(i,4))**2.0d0
        tmp4 = tmp4 + (xelem(i,2)-xelem(i,3))**2.0d0
        tmp5 = tmp5 + (xelem(i,2)-xelem(i,4))**2.0d0
        h = h + (xelem(i,3)-xelem(i,4))**2.0d0
      enddo


      if (tmp1.gt.h) then
        h=tmp1
      endif
      if (tmp2.gt.h)then
        h=tmp2
      endif
      if (tmp3.gt.h) then
        h=tmp3
      endif
      if (tmp4.gt.h) then
        h=tmp4
      endif
      if (tmp5.gt.h) then
        h=tmp5
      endif
      h=sqrt(h)

c     initialize variables
      estif(:,:) = 0.0d0 ! elemental consistent tangent stiffness
      eforc(:) = 0.0d0   ! elemental RHS/residual
      ident(:,:) = 0.0d0 ! identity matrix
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0
      
      call gausstet(l,ninte,sg,tg,zg,wg)     

c!$OMP PARALLEL DO
c!$OMP& PRIVATE(iinte,shap,xjaco,wtjac,mu,gamm,inode,temp,Fdef,j,
c!$OMP&  Fdet,Finv,pres,Ctens,r,Cdet,Cinv,Inv1,dJdC,tmp5,K1,tmp2,K2,tmp4,
c!$OMP&  tmp3,dWdC,SecPK,dCinvdCt,Ctang,d2JdC,Igeo,Lmat,q,k,t,Dgeomat,
c!$OMP&  bb,ievab,idofn,jevab,jnode,jdofn,Tp,FS)
c!$OMP& SHARED (ninte,ielem,sg,tg,zg,nnode,ndim,elnods,xelem,wg,
c!$OMP&  elemdata_nodal,s,tauMult,h,ident,uelem)
c!$OMP& REDUCTION(+:estif,eforc)! reduction on an array requires openmp2.5, reduction on allocatable openmp3.0
c!$acc region
      do iinte = 1,ninte ! for all Gauss integration points
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,
     $               xjaco,.false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco

        mu = 0.0d0    !!DS: value of mu at the current Gauss point
        gamm = 0.0d0  !!DS: value of gamm at the current Gauss point 
c	                 
        do inode = 1,nnode
          mu = mu + shap(4,inode)*elemdata_nodal(2,inode)
          gamm = gamm + shap(4,inode)*elemdata_nodal(1,inode)  
        enddo

c
        if (s.eq.0) then 
          temp  = 0.0d0 
        else if (s.eq.1) then
          temp = (0.5d0*tauMult*(h**2.0d0))/mu 
        else
          write(iwrit,200) 
 200      format(4x,'elem315.f: Stabilization property must be 0 or 1')  !!DS: element 315
          stop
        endif

c       compute the INVERSE of the deformation gradient at a Gauss Point, given information at 
c       the current configuration (!!DS: modified by Dawei)
        Finv(:,:) = ident(:,:)
        do inode = 1,nnode
          do j = 1,ndim
            do i = 1,ndim
              Finv(i,j)=Finv(i,j)-uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do
 
c       Finv_det is the determinant of the INVERSE of the deformation gradient
        Finv_det = Finv(1,1)*Finv(2,2)*Finv(3,3) +
     $         Finv(1,2)*Finv(2,3)*Finv(3,1) +
     $         Finv(1,3)*Finv(2,1)*Finv(3,2) -
     $         Finv(1,3)*Finv(2,2)*Finv(3,1) -
     $         Finv(1,2)*Finv(2,1)*Finv(3,3) -
     $         Finv(1,1)*Finv(2,3)*Finv(3,2)

        if (Finv_det .lt. 0.0d0) then
          Print*,"elem315.f: Negative Jacobian at element", ielem,  !!DS: elem315.f
     $           "Inverse of Jacob", Finv_det," ... exiting"


        negJac = .FALSE.
        endif

c       !!DS: Fdef is the deformation gradient at a Gauss Point
        Fdef(1,1) = Finv(2,2)*Finv(3,3)-Finv(2,3)*Finv(3,2)
        Fdef(1,2) = Finv(1,3)*Finv(3,2)-Finv(1,2)*Finv(3,3)
        Fdef(1,3) = Finv(1,2)*Finv(2,3)-Finv(1,3)*Finv(2,2)
        Fdef(2,1) = Finv(2,3)*Finv(3,1)-Finv(2,1)*Finv(3,3)
        Fdef(2,2) = Finv(1,1)*Finv(3,3)-Finv(1,3)*Finv(3,1)
        Fdef(2,3) = Finv(1,3)*Finv(2,1)-Finv(1,1)*Finv(2,3)
        Fdef(3,1) = Finv(2,1)*Finv(3,2)-Finv(2,2)*Finv(3,1)
        Fdef(3,2) = Finv(1,2)*Finv(3,1)-Finv(1,1)*Finv(3,2)
        Fdef(3,3) = Finv(1,1)*Finv(2,2)-Finv(1,2)*Finv(2,1)
        Fdef(:,:) = (1.0d0/Finv_det)*Fdef(:,:)

c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet=1.0d0/Finv_det


c       compute the spatial derivatives of the pressure and the pressure at the Gauss Point
c        (shap contain the derivative of the shape functions and its value)
        pres(1:4) = 0.0d0
        do i = 1, 4
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(4,inode)*shap(i,inode)
          end do
        end do

c       compute the left Cauchy-Green tensor      !!DS
        Btens(1:3,1:3) = 0.0d0
        do j = 1,ndim
          do i = j,ndim ! Btens is symmetric
            do r = 1,ndim
              Btens(i,j)= Btens(i,j)+Fdef(i,r)*Fdef(j,r)
            end do
          end do
        end do

c       enforce the symmetries
        Btens(1,2)=Btens(2,1)
        Btens(1,3)=Btens(3,1)
        Btens(2,3)=Btens(3,2)

c       principal invariant Inv1
        Inv1 = Btens(1,1)+Btens(2,2)+Btens(3,3)   !!DS
	  K1=Inv1/3.d0                              !!DS      
	  tmp5=mu*Fdet**(-5.d0/3.d0)*
     $	  dexp(gamm*(Fdet**(-2.d0/3.d0)*Inv1-3.d0))  !!DS  

c
        BFtens(1:3,1:3)=0.0d0                     !!DS
	  do i=1,ndim
	    do j=1,ndim
            do r=1,ndim
              BFtens(i,j)=BFtens(i,j)+Fdef(r,i)*Btens(r,j)
	      end do
	    end do
        end do


c       compute the Cauchy stress
        Sigma(1:3,1:3)=tmp5*(Btens(1:3,1:3)-K1*ident(1:3,1:3))-
     $	  pres(4)*ident(1:3,1:3)


c       compute the material tangent Ctang   !!DS: incompressible is different from compressible
        do i = 1, ndime
          do j = 1, ndime
            do q = 1, ndime
              do r = 1, ndime
 
                Ctang(i,j,q,r) =tmp5*(-5.d0/3.d0*(Btens(i,j)-
     $		  K1*ident(i,j))*Fdef(r,q)+Fdef(i,q)*Btens(j,r)+
     $          Fdef(j,q)*Btens(i,r)-2.d0/3.d0*ident(i,j)*BFtens(q,r)+
     $			  2.d0*gamm*Fdet**(-2.d0/3.d0)*(Btens(i,j)-
     $              K1*ident(i,j))*(-K1*Fdef(r,q)+BFtens(q,r)))
                         
              end do
		  end do
          end do
        end do       

c       construct the element stiffness matrix  !!DS: Here we do not use igeo and dgeomat
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,4 ! 4 = elemvec_ndofn(inode)
            ievab = ievab+1
            jevab = 0
            do jnode = 1, nnode
              do jdofn = 1,4 ! 4 = elemvec_ndofn(jnode)
                jevab = jevab+1



                if (idofn.le.3  .AND. jdofn.eq.4) then  !!DS: Compute the K12 matrix

                    estif(ievab,jevab)=estif(ievab,jevab)-
     $                shap(idofn,inode)*shap(4,jnode)*wtjac

                elseif (idofn.eq.4 .AND. jdofn.le.3) then !!DS: Compute the K21 matrix 
                 
                      do k = 1,ndim
                          estif(ievab,jevab)=estif(ievab,jevab)+  
     $                      shap(4,inode)*Fdet*Fdef(k,jdofn)* 						
     $					  shap(k,jnode)*wtjac
                      end do

                elseif (idofn.eq.4  .AND. jdofn.eq.4) then !!DS: Compute the K22 matrix

                  do i = 1,ndim
                        estif(ievab,jevab)=estif(ievab,jevab)+
     $                   temp*shap(i,inode)*shap(i,jnode)*wtjac
                  end do

                else  !!DS: Compute the K11 matrix

	            do j=1, ndim
	              do l=1,ndim
	                estif(ievab,jevab)=estif(ievab,jevab)+
     $                 shap(j,inode)*Ctang(idofn,j,jdofn,l)*
     $                 shap(l,jnode)*wtjac
                    end do
                  end do

                end if




              end do ! jdofn
            end do ! jnode
          end do ! idofn
        end do ! inode
 	  
c       create element residual for right hand side  !!DS: 
        Tp(1:4)=0.0d0

        do inode=1,4  !
          Tp(inode)=shap(4,inode)*(Fdet-1.d0)
	    do j=1,ndim
            Tp(inode)=Tp(inode)+temp*pres(j)*shap(j,inode)
	    end do
	  end do


        ievab = 0
        do inode = 1, nnode
          do i = 1,4 ! 4 = elemvec_ndofn(inode)
            ievab = ievab+1

            if (i.eq.4) then 
              eforc(ievab)=eforc(ievab)+Tp(inode)*wtjac
            else
              do j = 1, ndim
                 eforc(ievab)=eforc(ievab)+
     $                        shap(j,inode)*Sigma(i,j)*wtjac
              end do
            end if

          end do ! i over 4 dof
        end do ! inode



      end do ! iinte
c!$acc end region
c!$OMP END PARALLEL DO

      return

c     --------------------------------------------------------
c     Evaluate the elemantal RHS/residual from a source term
c     --------------------------------------------------------
 5    continue
      eforc(:) = 0.0d0 ! elemental RHS/residual
      return


c     --------------------------------------------------------
c     Build the elemental RHS/residual (eforc) for the dual problem
c     --------------------------------------------------------
 6    l  = int(pelem(1))  
      ndim = ndime

      eforc(:) = 0.0d0 ! elemental RHS/residual

      ievab = 0
      do inode = 1, nnode                             !!DS: Modified by Dawei
        do idofn = 1,4 ! 4 = elemvec_ndofn(inode)
          ievab = ievab+1
            ii = idnum(elnods(inode),idofn)
            if (ii.ne.0) then
              eforc(ievab) = -uelem_diff(idofn,inode)/
     $              float(nod_el_coun(ii+1)-nod_el_coun(ii))  ! twice jdofn ... ask Sevan JFD
            endif
        enddo ! idofn
      enddo ! inode

      return
c     --------------------------------------------------------
c     Compute the objective function (dataMatch+regularization)
c      and its gradient (egrad) on the element
c     --------------------------------------------------------
 7    l       = int(pelem(1))
      ireg    = int(pelem(2))
      alpha(1)= pelem(3)
      alpha(2)= pelem(4)
      beta(1) = pelem(5)
      beta(2) = pelem(6)
      Treg(1) = pelem(7)
      Treg(2) = pelem(8)
      tauMult = pelem(9)
      s       = int(pelem(10))

      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0
      egrad(:,:) = 0.0d0 ! initialize egrad

c     identity matrix ident
      ident(:,:) = 0.0d0
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0
      ndim = ndime

      powe = 0.5d0

      temp_dual(:) = 0.0d0 
      ii=0

      do inode = 1,nnode
        do idofn = 1,4 !4 = elemvec_ndofn(inode)  
          ii = ii+1
          temp_dual(ii)   = uelem_dual(idofn,inode)
        enddo
      enddo


c     determine the characteristic length h of the element
c       (for tetrahedra: h is the length of the longest edge)  !!DS: Need this for stabilization terms
      tmp1=(xelem(1,1)-xelem(1,2))**2.0d0
      tmp2=(xelem(1,1)-xelem(1,3))**2.0d0
      tmp3=(xelem(1,1)-xelem(1,4))**2.0d0
      tmp4=(xelem(1,2)-xelem(1,3))**2.0d0
      tmp5=(xelem(1,2)-xelem(1,4))**2.0d0
      h=(xelem(1,3)-xelem(1,4))**2.0d0


      do i=2,ndim ! 1, ndim
        tmp1 = tmp1 + (xelem(i,1)-xelem(i,2))**2.0d0
        tmp2 = tmp2 + (xelem(i,1)-xelem(i,3))**2.0d0
        tmp3 = tmp3 + (xelem(i,1)-xelem(i,4))**2.0d0
        tmp4 = tmp4 + (xelem(i,2)-xelem(i,3))**2.0d0
        tmp5 = tmp5 + (xelem(i,2)-xelem(i,4))**2.0d0
        h = h + (xelem(i,3)-xelem(i,4))**2.0d0
      enddo


      if (tmp1.gt.h) then
        h=tmp1
      endif
      if (tmp2.gt.h)then
        h=tmp2
      endif
      if (tmp3.gt.h) then
        h=tmp3
      endif
      if (tmp4.gt.h) then
        h=tmp4
      endif
      if (tmp5.gt.h) then
        h=tmp5
      endif
      h=sqrt(h)


      call gausstet(l,ninte,sg,tg,zg,wg)

      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,xjaco,
     $              .false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco

c       create prop_grad (contains prop gradient and value) 
        prop_grad(:,:)  = 0.0d0
        do iset = 1,2
          do idofn = 1,4      !three gradients+the value  !!DS: This is independent of elemvec_ndofn(inode)=3 or 4
            do inode = 1,nnode
               prop_grad(iset,idofn) = prop_grad(iset,idofn)+
     $                            shap(idofn,inode)*
     $                      elemdata_nodal(iset,inode)
            enddo
          enddo
        end do

c      compute the inverse of the deformation gradient and its determinant   !!DS: copy from itask=3
        Finv(:,:) = ident(:,:)
        do inode = 1,nnode
          do j = 1,ndim
            do i = 1,ndim
              Finv(i,j)=Finv(i,j)-uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do
 
c       Finv_det is the determinant of the INVERSE of the deformation gradient
        Finv_det = Finv(1,1)*Finv(2,2)*Finv(3,3) +
     $         Finv(1,2)*Finv(2,3)*Finv(3,1) +
     $         Finv(1,3)*Finv(2,1)*Finv(3,2) -
     $         Finv(1,3)*Finv(2,2)*Finv(3,1) -
     $         Finv(1,2)*Finv(2,1)*Finv(3,3) -
     $         Finv(1,1)*Finv(2,3)*Finv(3,2)

c       !!DS: Fdef is the deformation gradient at a Gauss Point
        Fdef(1,1) = Finv(2,2)*Finv(3,3)-Finv(2,3)*Finv(3,2)
        Fdef(1,2) = Finv(1,3)*Finv(3,2)-Finv(1,2)*Finv(3,3)
        Fdef(1,3) = Finv(1,2)*Finv(2,3)-Finv(1,3)*Finv(2,2)
        Fdef(2,1) = Finv(2,3)*Finv(3,1)-Finv(2,1)*Finv(3,3)
        Fdef(2,2) = Finv(1,1)*Finv(3,3)-Finv(1,3)*Finv(3,1)
        Fdef(2,3) = Finv(1,3)*Finv(2,1)-Finv(1,1)*Finv(2,3)
        Fdef(3,1) = Finv(2,1)*Finv(3,2)-Finv(2,2)*Finv(3,1)
        Fdef(3,2) = Finv(1,2)*Finv(3,1)-Finv(1,1)*Finv(3,2)
        Fdef(3,3) = Finv(1,1)*Finv(2,2)-Finv(1,2)*Finv(2,1)
        Fdef(:,:) = (1.0d0/Finv_det)*Fdef(:,:)

c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet=1.0d0/Finv_det

c       compute pressure and spatial derivatives at Gauss Point
        pres(1:4) = 0.0d0
        do i = 1, 4
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(4,inode)*shap(i,inode)
          end do
        end do

c       compute the left Cauchy-Green tensor      !!DS
        Btens(1:3,1:3) = 0.0d0
        do j = 1,ndim
          do i = j,ndim ! Btens is symmetric
            do r = 1,ndim
              Btens(i,j)= Btens(i,j)+Fdef(i,r)*Fdef(j,r)
            end do
          end do
        end do

c       enforce the symmetries
        Btens(1,2)=Btens(2,1)
        Btens(1,3)=Btens(3,1)
        Btens(2,3)=Btens(3,2)

c       principal invariant Inv1
        Inv1 = Btens(1,1)+Btens(2,2)+Btens(3,3)   !!DS
	  K1=Inv1/3.d0                             !!DS
	  tmp4=Fdet**(-2.d0/3.d0)*Inv1-3.d0      
	  tmp5=Fdet**(-5.d0/3.d0)*dexp(gamm*tmp4)  !!DS  

c
        do iset = 1,2
          do knode = 1,nnode

c           set variations in the shape function for the property 
            var_shape_prop(:,:) = 0.0d0
            if (ielem_nod_grad(iset,knode).eq.1) then
              do idofn = 1,4
                var_shape_prop(iset,idofn) = shap(idofn,knode)  !!DS: var_shape_prop is for a given optimization variable and node
              enddo
            endif

c
            if (s.eq.0) then 
              temp  = 0.0d0
              temp_der  = 0.0d0
            end if 
            if (s.eq.1) then 
              temp = 0.5d0*(tauMult*(h**2.0d0))/prop_grad(2,4)
              temp_der = 0.0d0
            end if
            if (s.eq.1 .AND. iset.eq.2) then 
              temp_der = -0.5d0*var_shape_prop(2,4)*
     $                (tauMult*(h**2.0d0))/(prop_grad(2,4)**2.0d0)
            end if

c
            Tp(1:4) = 0.0d0
            do i = 1, ndim
                Tp(i) = temp_der*pres(i)              
            end do                       !!DS: Tp(4)=0.0d0


c           compute derivative of Cauchy stress Sigma   
c           with respect to the material parameters    !!DS: modified by Dawei
            if (iset.eq.1) then
              sigma_grad(1:3,1:3) = mu*tmp5*tmp4*
     $			(Btens(1:3,1:3)-K1*ident(1:3,1:3))
            else if (iset.eq.2) then
              sigma_grad(1:3,1:3) =tmp5*(Btens(1:3,1:3)-
     $			K1*ident(1:3,1:3))

            end if

c       
            ievab = 0
            eforc(:) = 0.0d0
            do inode = 1,nnode
              do i = 1,4 ! 4 = elemvec_ndofn(inode)  
                ievab = ievab+1

                if (i.eq.4) then ! 4 = elemvec_ndofn(inode)
			     do j=1,3
				    eforc(ievab)=eforc(ievab)+
     $                         Tp(j)*shap(j,inode)*wtjac  
				 end do 
                else
			     do j = 1,3
                      eforc(ievab)=eforc(ievab)+
     $                         shap(j,inode)*sigma_grad(i,j)*wtjac
                   enddo
                end if

              enddo
            enddo ! inode
 
            do ii = 1, 16 !!DS: 4 node*4 dof=16   
              egrad(iset,knode) = egrad(iset,knode) +
     $                            eforc(ii)*temp_dual(ii)
            enddo

c           account for the regularization term
c JFD            if (iset.eq.1) then
              if (ireg.eq.0) then
c               there is not supposed to be a regularization term: do nothing
              elseif (ireg.eq.1) then ! H1 or Tikhonov reg.
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))
              elseif ((ireg.eq.2).or.(ireg.eq.21)) then ! TVD reg.
                deno = sqrt(beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $              +prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $              var_shape_prop(iset,2)*prop_grad(iset,2)+
     $              var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              elseif (ireg.eq.3) then ! power reg.
                deno = sqrt(beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $              +prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              ((1+deno-beta(iset))**(powe-1.0d0))*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              elseif (ireg.eq.31) then ! power reg. (proposed by Paul Barbone)
                deno = beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $             +prop_grad(iset,3)*prop_grad(iset,3)
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))*
     $              (deno**(1-0.5d0*powe)-(deno-beta(iset)*beta(iset))*
     $            (1-0.5d0*powe)*(deno**(-0.5d0*powe)))/(deno**(2-powe))
              elseif (ireg.eq.4) then ! logarithm reg.
                deno = sqrt(beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $              +prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/(deno*
     $              Treg(iset)*(1+(deno-beta(iset))/Treg(iset)))
              elseif (ireg.eq.41) then ! logarithm reg. (second implementation)
                deno = Treg(iset)*Treg(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2)+
     $              prop_grad(iset,3)*prop_grad(iset,3)
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              elseif (ireg.eq.5) then ! exponential reg. (contrast preserving)
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2)+
     $              prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))*
     $              exp((beta(iset)-deno)/Treg(iset))/(Treg(iset)*deno)
              elseif (ireg.eq.6) then ! fraction reg. (contrast preserving)
                deno = 1 + (prop_grad(iset,1)*prop_grad(iset,1)+
     $                      prop_grad(iset,2)*prop_grad(iset,2)+
     $                      prop_grad(iset,3)*prop_grad(iset,3))/
     $                     (Treg(iset)*Treg(iset))
                deno = Treg(iset)*Treg(iset)*deno*deno
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              else
                Print*,"elem315.f: ireg=",ireg,
     $                 " is not implemented: exiting"
                stop
              endif ! ireg
c            elseif (iset.eq.2) then
c              if (ireg.eq.1) then
c                egrad(iset,knode) = egrad(iset,knode) + alpha(iset)*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2)+    
c     $              var_shape_prop(iset,3)*prop_grad(iset,3))*
c     $              wtjac
c              elseif (ireg.eq.2) then
c                deno = sqrt(
c     $              beta(iset)*beta(iset)+
c     $              prop_grad(2,1)*prop_grad(2,1)+
c     $              prop_grad(2,2)*prop_grad(2,2)+
c     $              prop_grad(2,3)*prop_grad(2,3))
c                egrad(iset,knode) = egrad(iset,knode) + 
c     $              alpha(iset)*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2)+
c     $              var_shape_prop(iset,3)*prop_grad(iset,3))*
c     $              wtjac/deno
c              endif! ireg
c            endif! iset
          enddo ! knode
        enddo ! iset
      
        if (ireg.eq.1) then ! H1 or Tikhonov reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(2)*(prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3))+
     $      alpha(1)*(prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)))
        elseif (ireg.eq.2) then ! TVD reg. with offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3))+
     $      alpha(1)*sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)))
        elseif (ireg.eq.21) then ! TVD reg. without offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*(sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))-beta(2))+
     $      alpha(1)*(sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))-beta(1)))
        elseif (ireg.eq.3) then ! power reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)/powe*((1+sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))-beta(2))**(powe)-1.0d0)+
     $      alpha(1)/powe*((1+sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))-beta(1))**(powe)-1.0d0))
        elseif (ireg.eq.31) then ! power regularization (proposed by Paul Barbone)
c JFD: for powe.eq.1, this is another implementation of TVD
          deno = prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)+prop_grad(2,3)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(2)*deno/(beta(2)*beta(2)+deno)**(1-0.5d0*powe))
          deno = prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(1)*deno/(beta(1)*beta(1)+deno)**(1-0.5d0*powe))
        elseif (ireg.eq.4) then ! logarithm reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $       (alpha(2)*log(1+(sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))-beta(2))/Treg(2))+
     $       alpha(1)*log(1+(sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))-beta(1))/Treg(1)))
        elseif (ireg.eq.41) then ! logarithm reg. (second implementation to avoid using a sqrt operator and the associated beta that has a significant influence on the convergence speed)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $       (alpha(2)*log(1+(prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2)+prop_grad(2,3)*
     $           prop_grad(2,3))/(Treg(2)*Treg(2)))+
     $       alpha(1)*log(1+(prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)+prop_grad(1,3)*
     $           prop_grad(1,3))/(Treg(1)*Treg(1))))
        elseif (ireg.eq.5) then ! exponential reg. (contrast preserving: no penalty for large jumps)
c JFD: proposed by Paul Barbone
c JFD: this regularization is very unstable
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $       (alpha(2)*(1.0d0 - exp((beta(2)-sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $         prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3)))/Treg(2)))+
     $       alpha(1)*(1.0d0-exp((beta(1)-sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $         prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)))/Treg(1))))
        elseif (ireg.eq.6) then ! fraction reg. (contrast preserving: no penalty for large jumps)
c JFD: this regularization is very unstable
          deno = (prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3))/
     $           (Treg(2)*Treg(2))
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $                     (alpha(2)*deno/(1+deno))
          deno = (prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3))/
     $           (Treg(1)*Treg(1))
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $                     (alpha(1)*deno/(1+deno))
        endif
        l2grad1(ielem) = l2grad1(ielem) + sqrt(prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))
        l2grad2(ielem) = l2grad2(ielem) + sqrt(prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))
      enddo ! iinte

c     compute the gradient and the objective function
      udiff= 0.0d0
c-------------------------------------------------------------
c      do j=1,ndim
c        do inode=1,nnode
c          ii=idnum(elnods(inode),j)
c          if (ii.ne.0) then
c            udiff(j) = udiff(j)+uelem_diff(j,inode)/
c     $             sqrt(float(nod_el_coun(ii+1)-nod_el_coun(ii))) 
c          endif
c        enddo
c      enddo
c-------------------------------------------------------------
      do inode=1,nnode
        do j=1,ndim
          ii=idnum(elnods(inode),j)    
          if (ii.ne.0) then
	      udiff=udiff+1.d0/2.d0*(uelem_diff(j,inode))**2.d0/
     &            float(nod_el_coun(ii+1)-nod_el_coun(ii))           
          end if
        end do
	end do


      
      elemDataMatch(ielem) = udiff    !!DS: Modified by Dawei

      dataMatch=dataMatch+elemDataMatch(ielem)
      regularization = regularization + elemRegul(ielem)
      l2grad1(ielem) = l2grad1(ielem)/dble(ninte)
      l2grad2(ielem) = l2grad2(ielem)/dble(ninte)
      suml2grad1 = suml2grad1 + l2grad1(ielem)*l2grad1(ielem)
      suml2grad2 = suml2grad2 + l2grad2(ielem)*l2grad2(ielem)

      return



 4    return
 8    return
 9    return
 10   return
 11   return
 12   return
 13   return
 14   return
 15   return
 16   return
 17   return

      end