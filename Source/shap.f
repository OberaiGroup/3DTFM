c************************************************************
      subroutine gauss1(l,sg,wg)
c     gauss points and weights for one dimension
c     sg -> natural coordinates of gauss points
c     wg -> gauss weight
      USE IOUNIT
      implicit none
      integer l
      double precision t
      double precision sg(*),wg(*)
c------------------------------------------------------------
c     1-point integration
      if (l.eq.1) then
        sg(1) = 0.0d0
        wg(1) = 2.0d0

c     2-point integration
      else if (l.eq.2) then
        sg(1) = -1.0d0/sqrt(3.0d0)
        sg(2) = -sg(1)
        wg(1) = 1.0d0
        wg(2) = 1.0d0

c     3-point integration
      else if (l.eq.3) then
        sg(1) = -sqrt(0.6d0)
        sg(2) = 0.0d0
        sg(3) = -sg(1)
        wg(1) = 5.0d0/9.0d0
        wg(2) = 8.0d0/9.0d0
        wg(3) = wg(1)

c     4-point integration
      else if (l.eq.4) then
        t = sqrt(4.8d0)
        sg(1) = sqrt((3.0d0 + t)/7.0d0)
        sg(2) = sqrt((3.0d0 - t)/7.0d0)
        sg(3) = -sg(2)
        sg(4) = -sg(1)
        t = 1.0d0/3.0d0/t
        wg(1) = 0.5d0 - t
        wg(2) = 0.5d0 + t
        wg(3) = wg(2)
        wg(4) = wg(1)
      else
        write(iwrit,90) l
        call termin(6)
      end if
      return
90    format('** fatal error in subroutine gauss1** the set '
     +       'of coefficients for ',i1,
     +       ' quadrature points has not been implemented yet')
      end


c*************************************************************
      subroutine gauss2(l,ninte,sg,tg,wg)
c     gauss points and weights for two dimensions
c     sg,tg -> natural coordinates of gauss points
c     wg    -> gauss weight
      USE IOUNIT
      implicit none
      integer l,ninte,i,j,k
      double precision g,h
      double precision ls(9),lt(9),lw(9),g4(4),h4(4)
      double precision  h5(5),g5(5)
      double precision sg(*),tg(*),wg(*)
      data ls/-1.0d0,1.0d0,1.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0/
      data lt/-1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,0.0d0/
      data lw/4*25.0d0,4*40.0d0,64.0d0/
c------------------------------------------------------------
      ninte = l*l

c     1x1 integration
      if (l.eq.1) then
        sg(1) = 0.0d0
        tg(1) = 0.0d0
        wg(1) = 4.0d0

c     2x2 integration
      else if (l.eq.2) then
        g = 1.0d0 / sqrt(3.0d0)
        do i = 1,4
          sg(i) = g * ls(i)
          tg(i) = g * lt(i)
          wg(i) = 1.0d0
        enddo

c     3x3 integration
      else if (l.eq.3) then
        g = sqrt(0.6d0)
        h = 1.0d0/81.0d0
        do i = 1,9
          sg(i) = g * ls(i)
          tg(i) = g * lt(i)
          wg(i) = h * lw(i)
        enddo

c     4x4 integration
      else if (l.eq.4) then
        g = sqrt(4.8d0)
        h = sqrt(30.0d0)/36.0d0
        g4(1) = sqrt((3.0d0+g)/7.0d0)
        g4(4) = - g4(1)
        g4(2) = sqrt((3.0d0-g)/7.0d0)
        g4(3) = -g4(2)
        h4(1) = 0.5d0 - h
        h4(2) = 0.5d0 + h
        h4(3) = 0.5d0 + h
        h4(4) = 0.5d0 - h
        i = 0
        do j = 1,4
          do k = 1,4
            i = i + 1
            sg(i) = g4(k)
            tg(i) = g4(j)
            wg(i) = h4(j) * h4(k)
          enddo
        enddo

c     5x5 integration
      else if (l.eq.5) then
        g5(1) = -sqrt(245.0d0+14.0d0*sqrt(70.0d0))/21.0d0
        g5(2) = -sqrt(245.0d0-14.0d0*sqrt(70.0d0))/21.0d0
        g5(3) = 0.0d0
        g5(4) = -g5(2)
        g5(5) = -g5(1)
        h5(1) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
        h5(2) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
        h5(3) = 128.0d0/225.0d0
        h5(4) = h5(2)
        h5(5) = h5(1)
        i = 0
        do j = 1,5
          do k = 1,5
            i = i + 1
            sg(i) = g5(k)
            tg(i) = g5(j)
            wg(i) = h5(j) * h5(k)
          enddo
        enddo

c     3 point integration
      else if (l.lt.0) then
        ninte = 3
        g = sqrt(1.0d0/3.0d0)
        h = sqrt(2.0d0/3.0d0)
        sg(1) = -h
        sg(2) =  h
        sg(3) = 0.0d0
        tg(1) = -g
        tg(2) = -g
        tg(3) = g
        wg(1) = 1.0d0
        wg(2) = 1.0d0
        wg(3) = 2.0d0
      else
        write(iwrit,90) l,l
        call termin(6)
      end if
      return
90    format('** fatal error in subroutine gauss2 ** the set',
     +       ' of coefficients for ',i1,'x',i1,
     +       ' quadrature points has not been implemented yet')
      end


c*************************************************************
      subroutine gauss3(l,ninte,sg,tg,zg,wg)
c     gauss points and weights for three dimensions
c     sg,tg,zg -> natural coordinates of gauss points
c     wg    -> gauss weight
      USE IOUNIT
      implicit none
      integer l, ninte, i, j, k, m
      double precision g, h
      double precision ls(27),lt(27),lz(27),lw(27),g4(4),h4(4)
      double precision sg(*),tg(*),zg(*),wg(*)
      data ls/-1.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,
     +  -1.0d0,1.0d0,1.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,
     +  1.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0/
      data lt/-1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0,1.0d0,
     +  -1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,0.0d0,-1.0d0,
     +  0.0d0,1.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,0.0d0/
      data lz/4*-1.0d0,4*1.0d0,4*0.0d0,5*-1.0d0,5*1.0d0,5*0.0d0/
      data lw/8*125.0d0,8*200.0d0,320.0d0,4*200.0d0,320.0d0,4*320.0d0,
     +  512.0d0/
c------------------------------------------------------------
      ninte = l*l*l

c     1x1x1 integration
      if (l.eq.1) then
        sg(1) = 0.0d0
        tg(1) = 0.0d0
        zg(1) = 0.0d0
        wg(1) = 8.0d0

c     2x2x2 integration
      else if (l.eq.2) then
        g = 1.0d0/sqrt(3.0d0)
        do i = 1,8
          sg(i) = g * ls(i)
          tg(i) = g * lt(i)
          zg(i) = g * lz(i)
          wg(i) = 1.0d0
        enddo

c     3x3x3 integration
      else if (l.eq.3) then
        g = sqrt(0.6d0)
        h = 1.0d0/729.0d0
        do i = 1,27
          sg(i) = g * ls(i)
          tg(i) = g * lt(i)
          zg(i) = g * lz(i)
          wg(i) = h * lw(i)
        enddo

c     4x4x4 integration
      else if (l.eq.4) then
        g = sqrt(4.8d0)
        h = sqrt(30.0d0)/36.0d0
        g4(1) = sqrt((3.0d0+g)/7.0d0)
        g4(4) = - g4(1)
        g4(2) = sqrt((3.0d0-g)/7.0d0)
        g4(3) = -g4(2)
        h4(1) = 0.5d0 - h
        h4(2) = 0.5d0 + h
        h4(3) = 0.5d0 + h
        h4(4) = 0.5d0 - h
        i = 0
        do j = 1,4
          do k = 1,4
            do m = 1,4
              i = i + 1
              sg(i) = g4(m)
              tg(i) = g4(k)
              zg(i) = g4(j)
              wg(i) = h4(j) * h4(k) * h4(m)
            enddo
          enddo
        enddo

      else
        write(iwrit,90) l,l,l
        call termin(6)
      end if

      return
90    format('** fatal error in subroutine gauss3 ** the set of ',
     +       'coefficients for ',i1,'x',i1,'x',i1,
     +       ' quadrature points has not been implemented yet.')
       end


c*************************************************************
      subroutine gausstri(l,ninte,sg,tg,wg)
c     integration rules for triangles. assad oberai, taken from TJRHs book.
c     sg,tg -> natural coordinates of gauss points
c     wg    -> gauss weight
      USE IOUNIT
      implicit none
      integer l, ninte, iinte
      double precision s2(2),s3(3),w2(1),w3(1),s5(5),w5(2)
      double precision sg(*),tg(*),wg(*)
      data s2/0.666666666666667d0,0.166666666666667d0/,
     +     s3/0.659027622374092d0,0.231933368553031d0,
     +        0.109039009072877d0/,
     +     s5/0.797112651860071d0,0.165409927389841d0,
     +        0.037477420750088d0,0.124949503233232d0,
     +        0.437525248383384d0/,
     +     w2/0.333333333333333d0/,
     +     w3/0.1666666666666667d0/,
     +     w5/0.063691414286223d0,0.205950504760887d0/
c------------------------------------------------------------
c     deg. prec = 2
      if (l.eq.2) then
        ninte = 3
        sg(1) = s2(1)
        sg(2) = s2(2)
        sg(3) = s2(2)
        tg(1) = s2(2)
        tg(2) = s2(2)
        tg(3) = s2(1)
        do iinte = 1,ninte
          wg(iinte) = 0.5d0*w2(1)
        enddo
             
c     deg. prec = 3
      else if (l.eq.3) then
        ninte = 6
        sg(1) = s3(1)
        tg(1) = s3(2)
        sg(2) = s3(3)
        tg(2) = s3(1)
        sg(3) = s3(2)
        tg(3) = s3(3)
        sg(4) = s3(1)
        tg(4) = s3(3)
        sg(5) = s3(2)
        tg(5) = s3(1)
        sg(6) = s3(3)
        tg(6) = s3(2)
        do iinte = 1,ninte
          wg(iinte) =  0.5d0*w3(1)
        enddo

c     deg. prec = 5
      else if (l.eq.5) then
        ninte = 9
        sg(1) = s5(1)
        tg(1) = s5(2)
        sg(2) = s5(3)
        tg(2) = s5(1)
        sg(3) = s5(2)
        tg(3) = s5(3)
        sg(4) = s5(1)
        tg(4) = s5(3)
        sg(5) = s5(2)
        tg(5) = s5(1)
        sg(6) = s5(3)
        tg(6) = s5(2)
        sg(7) = s5(4)
        tg(7) = s5(4)
        sg(8) = s5(4)
        tg(8) = s5(5)
        sg(9) = s5(5)
        tg(9) = s5(4)
        do iinte = 1,6
          wg(iinte) =  0.5d0*w5(1)
        enddo
        do iinte = 7,ninte
          wg(iinte) =  0.5d0*w5(2)
        enddo

      else
        write(iwrit,90) l
        call termin(6)
      end if

      return
90    format('** fatal error in subroutine gausstri ** the set of ',
     +       ' coefficients for a precision of order ',i1,
     +       ' has not been implemented yet')
      end


c*************************************************************
      subroutine gausstet(l,ninte,sg,tg,zg,wg)
c     integration rules for tets
c     sg,tg,zg -> natural coordinates of integ. points
c     wg    ->  weight
      USE IOUNIT
      implicit none
      integer l, ninte, iinte
      double precision s4(3),w4(2)
      double precision sg(*),tg(*),zg(*),wg(*)
      data s4/0.25d0,0.5d0,0.166666666666667d0/,
     +     w4/-0.8d0,0.45d0/
c------------------------------------------------------------
c     deg. prec = 4
      if (l.eq.4) then
        ninte = 5
        sg(1) = s4(1)
        tg(1) = s4(1)
        zg(1) = s4(1)
        wg(1) = (1.0d0/6.0d0)*w4(1)
        sg(2) = s4(3)
        tg(2) = s4(3)
        zg(2) = s4(3)
        sg(3) = s4(2)
        tg(3) = s4(3)
        zg(3) = s4(3)
        sg(4) = s4(3)
        tg(4) = s4(2)
        zg(4) = s4(3)
        sg(5) = s4(3)
        tg(5) = s4(3)
        zg(5) = s4(2)
        do iinte = 2,ninte
          wg(iinte) = (1.0d0/6.0d0)*w4(2)
        enddo

      else
        write(iwrit,90) l
        call termin(6)
      end if

      return
90    format('** fatal error in subroutine gausstet ** the set of ',
     +       'coefficients for a precision of order ',i1,
     +       ' has not been implemented yet (only 4 works)')
      end


c JFD: the follwing function is not used ... it has been commented as
c      changes can not be tested/debugged!
c  c************************************************************
        subroutine shape1 (ss,shape,xjaco,deriv,nno,ndi,xel)
c  c     shape functions values at gauss points locations
c  c      for one-dimensional elements
c  c************************************************************
        USE IOUNIT
        USE MAINMEM! needed for nnode, ndime, xelem, ...
c        implicit double precision (a-h,o-z)
        implicit none
        double precision ss,xjaco,a,b,xa,xb
        integer nno,ndi,inode,j,idime
        double precision xel(ndi,*)! = xelem
c  c JFD ss and xjaco are not here?
c  c JFD would need to declare: a, b, xa, ...
        logical deriv
        double precision shape(2,2)
        double precision s(2)
        !dimension shape(2,*),s(2)
        !data s/-0.5d0,0.5d0/
c        !make sure all variables are declared  -TS
c  c-------------------------------------------------------------
        s(1) = -0.5d0
        s(2) =  0.5d0
c     form 2-node linear shape functions
        do inode = 1,2
          shape(2,inode) = 0.5d0 + s(inode) * ss ! function evaluation - TS
          shape(1,inode) = s(inode) ! derivative - TS
        enddo        
  

c     construct jacobian   
        a = 0.d0
        b = 0.d0
        do idime = 1,ndi
          xa = xel(idime,2) - xel(idime,1)
          a = a + xa * xa
        enddo
        a = sqrt(a)
        xjaco = a / 2.d0
  
c     form global derivatives
c        if (deriv) return
c        do inode = 1,nno
c          shape(1,inode) = shape(1,inode) / xjaco
c        enddo
  
        return
c  900   format(//2x,'** fatal error ** one-dimensional element has',
c       +    i2,' nodes.'/2x,'max nodes currentlty permissible is 3')
        end

c JFD ask sevan for this function ... as it is used only in a old version of 
c     elem707... write a better description of what it is doing
c     ***********************************************************
      subroutine shapehigher_deriv (ielem,ss,tt,shape,shapehod,
     $                             xjaco,deriv,nno,lno,ndi,xel)
c     ***********************************************************
      
c     shape function routine for two dimensional elements
      
c     shape(1,i) -> xl(1) derivative of shape function
c     shape(2,i) -> xl(2) derivative of shape function
c     shape(3,i) -> value of shape function
c     shapehod(1,i) -> derivative of shape function with respect to xx
c     shapehod(2,i) -> derivative of shape function with respect to yy
c     shapehod(3,i) -> derivative of shape function with respect to xy
c     xjaco      -> jacobian determinant
      
c******************************************************************
      USE IOUNIT
c JFD     USE MAINMEM! now passing arguments
      implicit none
      integer ielem,nno,ndi! = nnode, ndime
      double precision ss,tt,xjaco
      double precision shape(3,*)
      double precision shapehod(3,*)
      logical deriv
      integer lno(*)! = elnods
      double precision xel(ndi,*)! = xelem
      integer i,j,k,l,m
      double precision s(4),t(4)
      double precision xs(3,2),sx(2,2)
      double precision s2,t2,tp
      double precision shapeh(2,2,9)
      double precision xs2(2,2,2)
      double precision temp(3) 
      data s/-0.5d0,0.5d0,0.5d0,-0.5d0/,t/-0.5d0,-0.5d0,0.5d0,0.5d0/
c-------------------------------------------------------------------
      if (9.lt.nno) then
        print*,"shapehigher_deriv: nno is greater than 9: exiting"
        stop
      endif
c     form 4-node quadrilateral shape functions
      if (nno.eq.3) then
        shape(3,1) = ss
        shape(3,2) = tt
        shape(3,3) = 1.0d0-ss-tt
        shape(1,1) = 1.0d0
        shape(1,2) = 0.0d0
        shape(1,3) = -1.0d0
        shape(2,1) = 0.0d0
        shape(2,2) = 1.0d0
        shape(2,3) = -1.0d0
      else
        do i = 1,4
          shape(3,i) = (0.5d0+s(i)*ss)*(0.5d0+t(i)*tt)
          shape(1,i) = s(i)*(0.5d0+t(i)*tt)
          shape(2,i) = t(i)*(0.5d0+s(i)*ss)
        enddo
      endif

c     add quadratic terms if necessary
      if(nno.gt.4) then
        s2 = (1.0d0-ss*ss)/2.0d0
        t2 = (1.0d0-tt*tt)/2.0d0
        do i=5,9
          do j = 1,3
            shape(j,i) = 0.0d0
          enddo
        enddo

c       midside nodes (serendipity)
        if(lno(5).ne.0) then
          shape(1,5) = -ss*(1.0d0-tt)
          shape(2,5) = -s2
          shape(3,5) = s2*(1.0d0-tt)
        end if

        if(nno.lt.6) go to 210
        if(lno(6).ne.0) then
          shape(1,6) = t2
          shape(2,6) = -tt*(1.0d0+ss)
          shape(3,6) = t2*(1.0d0+ss)
        end if

        if(nno.lt.7) go to 210
        if(lno(7).ne.0) then
          shape(1,7) = -ss*(1.0d0+tt)
          shape(2,7) = s2
          shape(3,7) = s2*(1.0d0+tt)
        end if

        if(nno.lt.8) go to 210
        if(lno(8).ne.0) then
          shape(1,8) = -t2
          shape(2,8) = - tt*(1.0d0-ss)
          shape(3,8) = t2*(1.0d0-ss)
        end if

c       interior node (lagrangian)
        if(nno.lt.9) go to 210
        if(lno(9).ne.0) then
          shape(1,9) = -4.0d0*ss*t2
          shape(2,9) = -4.0d0*tt*s2
          shape(3,9) = 4.0d0*s2*t2

c         correct edge nodes for interior node (lagrangian)
          do j= 1,3
            do i = 1,4
              shape(j,i) = shape(j,i) - 0.25d0*shape(j,9)
            enddo
            do i = 5,8
              if(lno(i).ne.0) then
                shape(j,i) = shape(j,i) - 0.5d0*shape(j,9)
              end if
            enddo
          enddo
        end if

c       correct corner nodes for presense of midside nodes
210     k = 8
        do i = 1,4
          l = i + 4
          do j = 1,3
            shape(j,i) = shape(j,i) - 0.5d0*(shape(j,k) + shape(j,l))
          enddo
          k = l
        enddo
      endif!if(nno.gt.4) then

c     construct jacobian and its inverse
      do i = 1,ndi
        do j = 1,2
          xs(i,j) = 0.0d0
          do k = 1,nno
            xs(i,j) = xs(i,j) + xel(i,k)*shape(j,k)
          enddo
        enddo
      enddo

      if (ndi.eq.2) then
        xjaco = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
      else if (ndi.eq.3) then
        xjaco = sqrt( 
     $        (xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2))**2.0d0+
     $        (xs(3,1)*xs(1,2)-xs(1,1)*xs(3,2))**2.0d0+
     $        (xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1))**2.0d0)
      endif

      if (xjaco.le.0.0d0) then
        write(iwrit,*) '***** WARNING1 *****'
        write(iwrit,*) 'ielem=',ielem,'xjaco=',xjaco
        write(iwrit,*) '***** WARNING1 *****'! JFD make it an abort condition?
        stop
      endif

      if(deriv) return! JFD add warning message?
       
      if (ndi.eq.3) then
        write(iwrit,*)'elem77.f: global derivatives for 2D elements'
        write(iwrit,*)'in 3D not coded: exiting'
        call termin(8)
      endif

      sx(1,1) = xs(2,2)/xjaco
      sx(2,2) = xs(1,1)/xjaco
      sx(1,2) =-xs(1,2)/xjaco
      sx(2,1) =-xs(2,1)/xjaco

c     form global derivatives
      do i = 1,nno
        tp          = shape(1,i)*sx(1,1) + shape(2,i)*sx(2,1)
        shape(2,i)  = shape(1,i)*sx(1,2) + shape(2,i)*sx(2,2)
        shape(1,i)  = tp
      enddo

      do i = 1,nno
        shapeh(1,1,i) = 0.0d0
        shapeh(2,2,i) = 0.0d0
        shapeh(2,1,i) = t(i)*s(i)
        shapeh(1,2,i) = shapeh(2,1,i)
      enddo

      xs2(:,:,:) = 0.0d0
      do k = 1,nno
        xs2(1,1,1) = 0.0d0
        xs2(2,2,1) = xs2(2,2,1) + xel(2,k)*shapeh(1,2,k)
        xs2(1,2,1) = xs2(1,2,1) + xel(1,k)*shapeh(1,2,k)
        xs2(2,1,1) = 0.0d0
        xs2(1,1,2) = xs2(1,1,2) + xel(1,k)*shapeh(1,2,k)
        xs2(2,2,2) = 0.0d0
        xs2(1,2,2) = 0.0d0
        xs2(2,1,2) = xs2(2,1,2) + xel(2,k)*shapeh(1,2,k)
      enddo

      temp(1:3) = 0.0d0
      do m = 1, nno
        do j = 1, ndi
          do k = 1, ndi
            do l = 1, ndi
              shapehod(1,m)=temp(1)-shape(j,m)*xs2(j,k,l)*
     $                           sx(k,1)*sx(l,1)
              temp(1)=shapehod(1,m)
              shapehod(2,m)=temp(2)-shape(j,m)*xs2(j,k,l)*
     $                           sx(k,2)*sx(l,2)
              temp(2)=shapehod(2,m)
              shapehod(3,m)=temp(3)-shape(j,m)*xs2(j,k,l)*
     $                           sx(k,1)*sx(l,2)
              temp(3)=shapehod(3,m)
            enddo
            shapehod(1,m)=temp(1)+shapeh(k,j,m)*sx(k,1)*sx(j,1)
            temp(1)=shapehod(1,m)
            shapehod(2,m)=temp(2)+shapeh(k,j,m)*sx(k,2)*sx(j,2)
            temp(2)=shapehod(2,m)
            shapehod(3,m)=temp(3)+shapeh(k,j,m)*sx(k,1)*sx(j,2)
            temp(3)=shapehod(3,m)
          enddo
        enddo
        temp(1:3) = 0.0d0
      enddo

      return
      end
 

c****************************************************************
      subroutine shape2 (ielem,ss,tt,shape,xjaco,deriv,nno,lno,ndi,
     $    xel)
c     shape function routine for two dimensional elements
c     shape(1,i) -> xl(1) derivative of shape function
c     shape(2,i) -> xl(2) derivative of shape function
c     shape(3,i) -> value of shape function
c     xjaco      -> jacobian determinant

      USE IOUNIT
c JFD      USE MAINMEM! for lnods, ndime, xelem
      implicit none
      integer ielem
      logical deriv
      integer nno,ndi! = nnode, ndime
      integer lno(*)! = elnods
      double precision xjaco, ss, tt, s2, t2, tp
      double precision xel(ndi,*)! = xelem
      double precision shape(3,*),s(4),t(4),xs(3,2),sx(2,2)
      integer i,j,k,l
      data s/-0.5d0,0.5d0,0.5d0,-0.5d0/,t/-0.5d0,-0.5d0,0.5d0,0.5d0/
c-------------------------------------------------------------------
c     form the shape functions
      if (nno.eq.3) then
        shape(3,1) = ss
        shape(3,2) = tt
        shape(3,3) = 1.0d0-ss-tt
        shape(1,1) = 1.0d0
        shape(1,2) = 0.0d0
        shape(1,3) = -1.0d0
        shape(2,1) = 0.0d0
        shape(2,2) = 1.0d0
        shape(2,3) = -1.0d0
      else
        do i = 1,4
          shape(3,i) = (0.5d0+s(i)*ss)*(0.5d0+t(i)*tt)
          shape(1,i) = s(i)*(0.5d0+t(i)*tt)
          shape(2,i) = t(i)*(0.5d0+s(i)*ss)
        enddo
      endif

c     add quadratic terms if necessary
      if(nno.gt.4) then
        s2 = (1.0d0-ss*ss)/2.0d0
        t2 = (1.0d0-tt*tt)/2.0d0
        do i=5,9
          do j = 1,3
            shape(j,i) = 0.0d0
          enddo
        enddo

c       midside nodes (serendipity)
c JFD        if(lno(ielem,5).ne.0) then
        if(lno(5).ne.0) then
          shape(1,5) = -ss*(1.0d0-tt)
          shape(2,5) = -s2
          shape(3,5) = s2*(1.0d0-tt)
        end if

        if(nno.lt.6) go to 210

c JFD        if(lno(ielem,6).ne.0) then
        if(lno(6).ne.0) then
          shape(1,6) = t2
          shape(2,6) = -tt*(1.0d0+ss)
          shape(3,6) = t2*(1.0d0+ss)
        end if

        if(nno.lt.7) go to 210

c JFD        if(lno(ielem,7).ne.0) then
        if(lno(7).ne.0) then
          shape(1,7) = -ss*(1.0d0+tt)
          shape(2,7) = s2
          shape(3,7) = s2*(1.0d0+tt)
        end if

        if(nno.lt.8) go to 210

c JFD        if(lno(ielem,8).ne.0) then
        if(lno(8).ne.0) then
          shape(1,8) = -t2
          shape(2,8) = - tt*(1.0d0-ss)
          shape(3,8) = t2*(1.0d0-ss)
        end if

c       interior node (lagrangian)
        if(nno.lt.9) go to 210

c JFD        if(lno(ielem,9).ne.0) then
        if(lno(9).ne.0) then
          shape(1,9) = -4.0d0*ss*t2
          shape(2,9) = -4.0d0*tt*s2
          shape(3,9) = 4.0d0*s2*t2

c         correct edge nodes for interior node (lagrangian)
          do j= 1,3
            do i = 1,4
              shape(j,i) = shape(j,i) - 0.25d0*shape(j,9)
            enddo
            do i = 5,8
c JFD              if(lno(ielem,i).ne.0) then
              if(lno(i).ne.0) then
                shape(j,i) = shape(j,i) - 0.5d0*shape(j,9)
              end if
            enddo
          enddo
        end if

c       correct corner nodes for presense of midside nodes
210     k = 8
        do i = 1,4
          l = i + 4
          do j = 1,3
            shape(j,i) = shape(j,i) - 0.5d0*(shape(j,k) + shape(j,l))
          enddo
          k = l
        enddo
      end if! nno.gt.4

c     construct jacobian and its inverse
      do i = 1,ndi
        do j = 1,2
          xs(i,j) = 0.0d0
          do k = 1,nno
            xs(i,j) = xs(i,j) + xel(i,k)*shape(j,k)
          enddo
        enddo
      enddo

      if (ndi.eq.2) then
        xjaco = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
      else if (ndi.eq.3) then
        xjaco = sqrt( 
     $         (xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2))**2.0d0+
     $         (xs(3,1)*xs(1,2)-xs(1,1)*xs(3,2))**2.0d0+
     $         (xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1))**2.0d0)
      endif

      if (xjaco.le.0.d0) then
        write(iwrit,*) '***** WARNING2 *****'
        write(iwrit,*) 'ielem=',ielem,'xjaco=',xjaco
        write(iwrit,*) '***** WARNING2 *****'
        stop
      endif

      if(deriv) return
       
      if (ndi.eq.3) then
        write(iwrit,*)'Error: global derivatives for 2D elements'
        write(iwrit,*)'in 3D not coded'
        call termin(8)
      endif

      sx(1,1) = xs(2,2)/xjaco
      sx(2,2) = xs(1,1)/xjaco
      sx(1,2) =-xs(1,2)/xjaco
      sx(2,1) =-xs(2,1)/xjaco

c     form global derivatives
      do i = 1,nno
        tp          = shape(1,i)*sx(1,1) + shape(2,i)*sx(2,1)
        shape(2,i)  = shape(1,i)*sx(1,2) + shape(2,i)*sx(2,2)
        shape(1,i)  = tp
      enddo
      return
      end


c     ***********************************************************
      subroutine shape3 (ielem,ss,tt,zz,shape,xjaco,deriv,nno,ndi,
     $     lno,xel)
c     shape function routine for three dimensional elements

c     if deriv=true: 
c         shape(1,i) -> ss derivative of shape function
c         shape(2,i) -> tt derivative of shape function
c         shape(3,i) -> zz derivative of shape function
c     if deriv=false: 
c         shape(1,i) -> xl(1) derivative of shape function
c         shape(2,i) -> xl(2) derivative of shape function
c         shape(3,i) -> xl(3) derivative of shape function
c         shape(4,i) -> value of shape function
c     xjaco      -> jacobian determinant
 
c     nodal natural coordinates:
c       1(-1,-1,-1)
c       2( 1,-1,-1)
c       3( 1, 1,-1)
c       4(-1, 1,-1)
c       5(-1,-1, 1)
c       6( 1,-1, 1)
c       7( 1, 1, 1)
c       8(-1, 1, 1)
c       9( 0,-1,-1)
c      10( 1, 0,-1)
c      11( 0, 1,-1)
c      12(-1, 0,-1)
c      13(-1,-1, 0)
c      14( 1,-1, 0)
c      15( 1, 1, 0)
c      16(-1, 1, 0)
c      17( 0,-1, 1)
c      18( 1, 0, 1)
c      19( 0, 1, 1)
c      20(-1, 0, 1)
c      21( 0, 0,-1)
c      22( 0,-1, 0)
c      23( 1, 0, 0)
c      24( 0, 1, 0)
c      25(-1, 0, 0)
c      26( 0, 0, 1)
c      27( 0, 0, 0)

      USE IOUNIT
c JFD      USE MAINMEM! for ndime, lnods, xelem
      implicit none
      integer ielem
      double precision ss, tt, zz, xjaco
      logical deriv
      integer nno,ndi! = nnode, ndime, nelem
      integer lno(*)! = elnods
      double precision xel(ndi,*)! = xelem
      double precision shape(4,*)
      double precision s(27),t(27),z(27)
      double precision xs(3,3),sx(3,3)
      integer inode,idime,iii,ii,i,j,k,l,m,n
      double precision s2, t2, z2, ajaco, tp1, tp2
      data s/-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0,
     +        0.0d0, 0.5d0, 0.0d0,-0.5d0,
     +       -0.5d0, 0.5d0, 0.5d0,-0.5d0, 0.0d0, 0.5d0, 0.0d0,-0.5d0,
     +        0.0d0, 0.0d0, 0.5d0, 0.0d0,-0.5d0, 0.0d0, 0.0d0/,
     +     t/-0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0,
     +       -0.5d0, 0.0d0, 0.5d0, 0.0d0,
     +       -0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0, 0.0d0, 0.5d0, 0.0d0,
     +        0.0d0,-0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0/,
     +     z/4*-0.5d0,4*0.5d0,4*-0.5d0,4*0.0d0,4*0.5d0,-0.5d0,
     +        4*0.0d0,0.5d0,0.0d0/
c-----------------------------------------------------------------------
c     construct the linear tetrahedral element (assad, aug 97.)
      if (nno.eq.4) then
        do inode = 1,nno
          do idime = 1,(ndi+1)
            shape(idime,inode)=0.0d0
          enddo
        enddo
        do ii = 1,3
          shape(ii,ii)= 1.0d0
          shape(ii,4) = -1.0d0
        enddo
        shape(4,1) = ss
        shape(4,2) = tt
        shape(4,3) = zz
        shape(4,4) = 1.0d0-ss-tt-zz
      endif

      if (nno.lt.5) go to 800

c     form 8-node hexahedral shape functions
      do i = 1,8
        shape(4,i) = (0.5d0 + s(i)*ss)*(0.5d0 + t(i)*tt)*(0.5 + z(i)*zz)
        shape(1,i) =            s(i)*(0.5d0 + t(i)*tt)*(0.5d0 + z(i)*zz)
        shape(2,i) =            t(i)*(0.5d0 + s(i)*ss)*(0.5d0 + z(i)*zz)
        shape(3,i) =            z(i)*(0.5d0 + s(i)*ss)*(0.5d0 + t(i)*tt)
      enddo

      if (nno.lt.9) go to 800

c     form midside nodes 9,10,11,&12 on zz=-1 : 
c     12-node transition hexahedron

      s2 = (1.0d0 - ss**2)
      t2 = (1.0d0 - tt**2)
      z2 = (1.0d0 - zz**2)
c JFD      if (lno(ielem,9).ne.0
c JFD     +     .and.lno(ielem,10).ne.0
c JFD     +     .and.lno(ielem,11).ne.0
c JFD     +     .and.lno(ielem,12).ne.0) then
      if (lno(9).ne.0
     +     .and.lno(10).ne.0
     +     .and.lno(11).ne.0
     +     .and.lno(12).ne.0) then

        do i = 9,12
          do j = 1,4
            shape(j,i) = 0.0d0
          enddo
        enddo

        do i = 9,11,2
          ii = i + 1
          shape(4,i) = s2*(0.5d0 + t(i)*tt)*(0.5d0 + z(i)*zz)
          shape(1,i) = -2.0d0*ss*(0.5d0 + t(i)*tt)*(0.5d0 + z(i)*zz)
          shape(2,i) = s2*t(i)*(0.5d0 + z(i)*zz)
          shape(3,i) = s2*(0.5d0 + t(i)*tt)*z(i)
          shape(4,ii) = t2*(0.5d0 + s(ii)*ss)*(0.5d0 + z(ii)*zz)
          shape(1,ii) = t2*s(ii)*(0.5d0 + z(ii)*zz)
          shape(2,ii) = -2.0d0*tt*(0.5d0 + s(ii)*ss)*(0.5d0 + z(ii)*zz)
          shape(3,ii) = t2*(0.5d0 + s(ii)*ss)*z(ii)
        enddo

c       correct nodes 1-4 for presence of nodes 9-12
        k = 12
        do i = 1,4
          l = i + 8
          do j = 1,4
            shape(j,i) = shape(j,i) - 0.5d0*(shape(j,k)+shape(j,l))
          enddo
          k = l
        enddo

      else! lno(ielem,9).ne.0 .and. ...
        do inode = 9,12
c JFD          if (lno(ielem,inode).eq.0) write(iwrit,990)nno,
c JFD     $           ielem,inode
          if (lno(inode).eq.0) write(iwrit,990)nno,
     $           ielem,inode
        enddo
        call termin(7)

      endif

      if (nno.lt.13) go to 800

c     form interior node 13 (location of 21) on zz=-1: 
c     13-node transition hexahedron
      if (nno.eq.13) then
c JFD        if (lno(ielem,13).ne.0) then
        if (lno(13).ne.0) then
          do j= 1,4
            shape(j,13) = 0.0d0
          enddo

          shape(4,13) = s2*t2*(0.5d0 + z(21)*zz)
          shape(1,13) = -2.0d0*ss*t2*(0.5d0 + z(21)*zz)
          shape(2,13) = s2*(-2.0d0*tt)*(0.5d0 + z(21)*zz)
          shape(3,13) = s2*t2*z(21)

c         correct nodes 1-4,& 9-12 for presence of node 13
          do j= 1,4
            do i= 1,4
             shape(j,i) = shape(j,i) + 0.25d0*shape(j,13)
             shape(j,i+8) = shape(j,i+8) - 0.5d0*shape(j,13)
            enddo
          enddo
          go to 800
        else! lno(ielem,13).ne.0
          inode = 13
          write(iwrit,990) nno,ielem,inode
          call termin(7)
        endif
      endif

c     form edge nodes 13-20 : 20-node hexahedron
c JFD      if (nno.ge.20.and.lno(ielem,13).ne.0
c JFD     +             .and.lno(ielem,14).ne.0
c JFD     +             .and.lno(ielem,15).ne.0
c JFD     +             .and.lno(ielem,16).ne.0
c JFD     +             .and.lno(ielem,17).ne.0
c JFD     +             .and.lno(ielem,18).ne.0
c JFD     +             .and.lno(ielem,19).ne.0
c JFD     +             .and.lno(ielem,20).ne.0) then
      if (nno.ge.20.and.lno(13).ne.0
     +             .and.lno(14).ne.0
     +             .and.lno(15).ne.0
     +             .and.lno(16).ne.0
     +             .and.lno(17).ne.0
     +             .and.lno(18).ne.0
     +             .and.lno(19).ne.0
     +             .and.lno(20).ne.0) then

        do i = 13,20
          do j = 1,4
            shape(j,i) = 0.0d0
          enddo
        enddo

        do i = 13,16
          shape(4,i) = z2*(0.5d0 + s(i)*ss)*(0.5d0 + t(i)*tt)
          shape(1,i) = z2*s(i)*(0.5d0 + t(i)*tt)
          shape(2,i) = z2*(0.5d0 + s(i)*ss)*t(i)
          shape(3,i) = -2.0d0*zz*(0.5d0 + s(i)*ss)*(0.5d0 + t(i)*tt)
        enddo
        do i = 17,19,2
          ii = i+1
          shape(4,i) = s2*(0.5d0 + t(i)*tt)*(0.5d0 + z(i)*zz)
          shape(1,i) = -2.0d0*ss*(0.5d0 + t(i)*tt)*(0.5d0 + z(i)*zz)
          shape(2,i) = s2*t(i)*(0.5d0 + z(i)*zz)
          shape(3,i) = s2*(0.5d0 + t(i)*tt)*z(i)
          shape(4,ii) = t2*(0.5d0 + s(ii)*ss)*(0.5d0 + z(ii)*zz)
          shape(1,ii) = t2*s(ii)*(0.5d0 + z(ii)*zz)
          shape(2,ii) = -2.0d0*tt*(0.5d0 + s(ii)*ss)*(0.5d0 + z(ii)*zz)
          shape(3,ii) = t2*(0.5d0 + s(ii)*ss)*z(ii)
        enddo

c       correct nodes 5-8 for presence of nodes 17-20
c       and nodes 1-8 for presence of nodes 13-16
        k = 20
        do i = 5,8
          l = i+12
          do j = 1,4
            shape(j,i) = shape(j,i) - 0.5d0*(shape(j,k)+shape(j,l))
          enddo
          k = l
        enddo

        do i = 1,4
          ii = i + 4
          k  = i + 12
          do j = 1,4
            shape(j,i) = shape(j,i) - 0.5d0*shape(j,k)
            shape(j,ii) = shape(j,ii) - 0.5d0*shape(j,k)
          enddo
        enddo

      else! nno.ge.20.and.lno(ielem,13).ne.0 .and. ...
        do inode = 13,20
c JFD          if (lno(ielem,inode).eq.0) write(iwrit,990) nno,ielem,inode
          if (lno(inode).eq.0) write(iwrit,990) nno,ielem,inode
        enddo
        call termin(7)

      endif

      if (nno.lt.21) go to 800
      
c     form interior node 21 on zz=-1: 21-node transition hexahedron
c JFD      if (lno(ielem,21).ne.0) then
      if (lno(21).ne.0) then
        do j= 1,4
          shape(j,21) = 0.0d0
        enddo

        shape(4,21) = s2*t2*(0.5d0 + z(21)*zz)
        shape(1,21) = -2.0d0*ss*t2*(0.5d0 + z(21)*zz)
        shape(2,21) = s2*(-2.0d0*tt)*(0.5d0 + z(21)*zz)
        shape(3,21) = s2*t2*z(21)

c       correct nodes 1-4,& 9-12 for presence of node 21
        do j= 1,4
          do i= 1,4
            shape(j,i) = shape(j,i) + 0.25d0*shape(j,21)
            shape(j,i+8) = shape(j,i+8) - 0.5d0*shape(j,21)
          enddo
        enddo
      else! lno(ielem,21).ne.0
        inode = 21
        write(iwrit,990) nno,ielem,inode
        call termin(7)

      endif

      if (nno.eq.21) go to 800

c     form 26-node hexahedron
c JFD      if (lno(ielem,22).ne.0
c JFD     +     .and.lno(ielem,23).ne.0
c JFD     +     .and.lno(ielem,24).ne.0
c JFD     +     .and.lno(ielem,25).ne.0
c JFD     +     .and.lno(ielem,26).ne.0) then
      if (lno(22).ne.0
     +     .and.lno(23).ne.0
     +     .and.lno(24).ne.0
     +     .and.lno(25).ne.0
     +     .and.lno(26).ne.0) then

        do i = 22,26
          do j = 1,4
            shape(j,i) = 0.0d0
          enddo
        enddo

        do i = 22,24,2
          ii = i + 1
          shape(4,i) = s2*z2*(0.5d0 + t(i)*tt)
          shape(1,i) = -2.0d0*ss*z2*(0.5d0 + t(i)*tt)
          shape(2,i) = s2*z2*t(i)
          shape(3,i) = s2*(-2.0d0*zz)*(0.5d0 + t(i)*tt)
          shape(4,ii) = t2*z2*(0.5d0 + s(ii)*ss)
          shape(1,ii) = t2*z2*s(ii)
          shape(2,ii) = -2.0d0*tt*z2*(0.5d0 + s(ii)*ss)
          shape(3,ii) = t2*(-2.0d0*zz)*(0.5d0 + s(ii)*ss)
        enddo
        shape(4,26) = s2*t2*(0.5d0 + z(26)*zz)
        shape(1,26) = -2.0*ss*t2*(0.5d0 + z(26)*zz)
        shape(2,26) = s2*(-2.0*tt)*(0.5d0 + z(26)*zz)
        shape(3,26) = s2*t2*z(26)

c       correct nodes 1,2,5,6, and  9,13,14,17 for presence of node 22
c        and nodes 2,3,6,7, and 10,14,15,18 for presence of node 23
c        and nodes 3,4,7,8, and 11,15,16,19 for presence of node 24
c        and nodes 1,4,5,8, and 12,13,16,20 for presence of node 25
c        and nodes 5,6,7,8, and 17,18,19,20 for presence of node 26

        do j = 1,4
          k = 25
          l = 22
          do i = 1,4
            ii  = i + 4
            iii = i + 12
            shape(j,i)  = shape(j,i)  + 0.25d0*(shape(j,k)+shape(j,l))
            shape(j,ii) = shape(j,ii) + 0.25d0*(shape(j,k)+shape(j,l))
            shape(j,iii)= shape(j,iii) - 0.5d0*(shape(j,k)+shape(j,l))
            k = l
            l = k + 1
          enddo

          do i = 9,12
            ii = i + 8
            k  = i + 13
            shape(j,i)  = shape(j,i) - 0.5d0*shape(j,k)
            shape(j,ii) = shape(j,ii) - 0.5d0*shape(j,k)
          enddo
 
          do i = 5,8
            ii = i + 12
            shape(j,i)  = shape(j,i) + 0.25d0*shape(j,26)
            shape(j,ii) = shape(j,ii) - 0.5d0*shape(j,26)
          enddo
        enddo
      else! lno(ielem,22).ne.0 .and. ...
        do inode = 22,26
c JFD          if (lno(ielem,inode).eq.0) write(iwrit,990)nno,ielem,inode
          if (lno(inode).eq.0) write(iwrit,990)nno,ielem,inode
        enddo
        call termin(7)
      endif

      if (nno.lt.27) go to 800

c     form central node : 27-node hexahedron
c JFD      if (nno.eq.27.and.lno(ielem,27).ne.0) then
      if (nno.eq.27.and.lno(27).ne.0) then
        do j = 1,4
          shape(j,27) = 0.0d0
        enddo
        shape(4,27) = s2*t2*z2
        shape(1,27) = -2.0d0*ss*t2*z2
        shape(2,27) = -2.0d0*s2*tt*z2
        shape(3,27) = -2.0d0*s2*t2*zz

c       correct nodes 1-26 for presence of node 27
        do j = 1,4
          do i = 1,8
            shape(j,i) = shape(j,i) - 0.125d0*shape(j,27)
          enddo
          do i = 9,20
            shape(j,i) = shape(j,i) + 0.25d0*shape(j,27)
          enddo
          do i = 21,26
            shape(j,i) = shape(j,i) - 0.50d0*shape(j,27)
          enddo
        enddo
      else! nno.eq.27.and.lno(ielem,27).ne.0
        inode = 27
        write(iwrit,990) nno,ielem,inode
        call termin(7)
      endif

c     construct jacobian xs(j,i) = dx,i / d exi,j and its determinant
800   do i = 1,ndi
        do j = 1,ndi
          xs(j,i) = 0.0d0
          do k = 1,nno
            xs(j,i) = xs(j,i) + xel(i,k)*shape(j,k)
          enddo
        enddo
      enddo

      xjaco = xs(1,1)*xs(2,2)*xs(3,3) + xs(1,2)*xs(2,3)*xs(3,1)
     +      + xs(1,3)*xs(2,1)*xs(3,2) - xs(1,1)*xs(2,3)*xs(3,2)
     +      - xs(1,2)*xs(2,1)*xs(3,3) - xs(1,3)*xs(2,2)*xs(3,1)

      if (xjaco.le.0.d0) then
        write(iwrit,*) '***** WARNING shape3 *****'
        write(iwrit,*) 'ielem=',ielem,'xjaco=',xjaco
        print*,"nno = ",nno,", ndi = ",ndi
        print*,"xel = "
        print*,xel(1,1),", ",xel(1,2),", ",xel(1,3),", ",xel(1,4)
        print*,xel(2,1),", ",xel(2,2),", ",xel(2,3),", ",xel(2,4)
        print*,xel(3,1),", ",xel(3,2),", ",xel(3,3),", ",xel(3,4)
        write(iwrit,*) '***** WARNING shape3 *****'
        stop
      endif

      if (deriv) return

c     construct inverse jacobian sx(j,i) = d exi,i / dx,j
      do i = 1,ndi
        j = i + 1
        k = i - 1
        if (j.gt.3) j = 1
        if (k.lt.1) k = 3
        do l = 1,ndi
          m = l + 1
          n = l - 1
          if (m.gt.3) m = 1
          if (n.lt.1) n = 3
          sx(i,l) = (xs(m,j)*xs(n,k) - xs(m,k)*xs(n,j))/xjaco
        enddo
      enddo

c     xjaco = surface jacobian (instead of volume jacobian) if zz=-1
      if (zz.eq.-1.0) then
        ajaco = 0.0d0
        do i = 1,ndi
          ajaco = ajaco + sx(i,3)**2.0d0
        enddo
        xjaco = sqrt(ajaco)*xjaco
      endif

c     form global derivatives
      do i = 1,nno
        tp1       =shape(1,i)*sx(1,1) + shape(2,i)*sx(1,2) + 
     +             shape(3,i)*sx(1,3)
        tp2       =shape(1,i)*sx(2,1) + shape(2,i)*sx(2,2) + 
     +             shape(3,i)*sx(2,3)
        shape(3,i)=shape(1,i)*sx(3,1) + shape(2,i)*sx(3,2) + 
     +             shape(3,i)*sx(3,3)
        shape(1,i)=tp1
        shape(2,i)=tp2
      enddo

      return

990   format(//'ERROR FOUND IN SHAPE3 :'/'nno = ',
     +        i5,5x,'and lnods(',i5,','i5,') = 0'//)
      end


c****************************************************************
      subroutine shape2ho (ielem,ss,tt,shape,xjaco,deriv,nno,lno,ndi,
     $    xel)
c     shape function routine for two dimensional elements with higher-order structure
c     shape(1,i) -> xl(1) derivative of shape function
c     shape(2,i) -> xl(2) derivative of shape function
c     shape(3,i) -> value of shape function
c     xjaco      -> jacobian determinant

      USE IOUNIT
c JFD      USE MAINMEM! for lnods, ndime, xelem
      implicit none
      integer ielem
      logical deriv
      integer nno,ndi! = nnode, ndime
      integer lno(*)! = elnods
      double precision xjaco, ss, tt, tp
      double precision xel(ndi,*)! = xelem
      double precision shape(3,*),xs(3,2),sx(2,2)
      integer i,j,k,l,order
      double precision tmp1, tmp2
      double precision s(20)! warning 20 should be enough...
c-------------------------------------------------------------------
c     form the shape functions
      if ((nno.eq.3).or.(nno.eq.6).or.(nno.eq.10).or.(nno.eq.15)) then
        print*,'Source/shap.f: shape2ho: not implemented yet ',
     &          'for triangles ... exiting'
        stop
      endif

      if (nno.eq.4) then
        order=2
      elseif (nno.eq.9) then
        order=3
      elseif (nno.eq.16) then
        order=4
      elseif (nno.eq.25) then
        order=5
      elseif (nno.eq.36) then
        order=6
      else
        print*,'Source/shap.f: shape2ho: not implemented yet ',
     &          'for such high order ... exiting'
        stop
      endif
 
c      print*,"computing for gp: ",ss,tt
c      print*,"nno=",nno,", order=",order

      do i=1,order
        s(i)=-1.0d0+2.0d0*dble(i-1)/dble(order-1)
      enddo
c      print*,"s=",s
      
      ! use lagrange polynomials for the function and its derivatives
      do i=1,order
        do j=1,order
          shape(3,i+(j-1)*order)=1.0d0
          shape(1,i+(j-1)*order)=0.0d0
          shape(2,i+(j-1)*order)=1.0d0
          do k=1,order!i2
            if (k.ne.i) then
              shape(3,i+(j-1)*order)=shape(3,i+(j-1)*order)*
     $           (ss-s(k))/(s(i)-s(k))! assumes equidistant nodes (non spectral elem)
              shape(2,i+(j-1)*order)=shape(2,i+(j-1)*order)*
     $           (ss-s(k))/(s(i)-s(k))! assumes equidistant nodes (non spectral elem)
              tmp1=1.0d0
              do l=1,order! of the derivative
               if (l.ne.i) then
                if (l.ne.k) then
                  tmp1=tmp1*(ss-s(l))/(s(i)-s(l))
                else
                  tmp1=tmp1/(s(i)-s(l))
                endif
               endif
              enddo
              shape(1,i+(j-1)*order)=shape(1,i+(j-1)*order)+tmp1
            endif
          enddo
          tmp2=0.0d0
          do k=1,order!j2
            if (k.ne.j) then
              shape(3,i+(j-1)*order)=shape(3,i+(j-1)*order)*
     $           (tt-s(k))/(s(j)-s(k))! assumes equidistant nodes (non spectral elem)
              shape(1,i+(j-1)*order)=shape(1,i+(j-1)*order)*
     $           (tt-s(k))/(s(j)-s(k))! assumes equidistant nodes (non spectral elem)
              tmp1=1.0d0
              do l=1,order
               if (l.ne.j) then
                if (l.ne.k) then
                  tmp1=tmp1*(ss-s(l))/(s(j)-s(l))
                else
                  tmp1=tmp1/(s(j)-s(l))
                endif
               endif
              enddo
              tmp2=tmp2+tmp1
            endif
          enddo
          shape(2,i+(j-1)*order)=shape(2,i+(j-1)*order)*tmp2
        enddo
      enddo
c      print*,shape(1,1),shape(1,2),shape(1,3),shape(1,4),shape(1,5),
c     $  shape(1,6),shape(1,7),shape(1,8),shape(1,9)
c      print*,shape(2,1),shape(2,2),shape(2,3),shape(2,4),shape(2,5),
c     $  shape(2,6),shape(2,7),shape(2,8),shape(2,9)
c      print*,shape(3,1),shape(3,2),shape(3,3),shape(3,4),shape(3,5),
c     $  shape(3,6),shape(3,7),shape(3,8),shape(3,9)
c      print*,xel(1,1),xel(1,2),xel(1,3),xel(1,4),xel(1,5),xel(1,6),
c     $  xel(1,7),xel(1,8),xel(1,9)
c      print*,xel(2,1),xel(2,2),xel(2,3),xel(2,4),xel(2,5),xel(2,6),
c     $  xel(2,7),xel(2,8),xel(2,9)
c      stop

c        do i = 1,4
c          shape(3,i) = (0.5+s(i)*ss)*(0.5+t(i)*tt)
c          shape(1,i) = s(i)*(0.5+t(i)*tt)
c          shape(2,i) = t(i)*(0.5+s(i)*ss)
c        enddo

      
c     construct jacobian and its inverse
      do i = 1,ndi
        do j = 1,2
          xs(i,j) = 0.0d0
          do k = 1,nno
            xs(i,j) = xs(i,j) + xel(i,k)*shape(j,k)
          enddo
        enddo
      enddo

      if (ndi.eq.2) then
        xjaco = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
      else if (ndi.eq.3) then
        xjaco = sqrt( 
     $         (xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2))**2.0d0+
     $         (xs(3,1)*xs(1,2)-xs(1,1)*xs(3,2))**2.0d0+
     $         (xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1))**2.0d0)
      endif

      if (xjaco.le.0.d0) then
        write(iwrit,*) '***** WARNING4 *****'
        write(iwrit,*) 'ielem=',ielem,'xjaco=',xjaco
        write(iwrit,*) '***** WARNING4 *****'
        stop
      endif

      if(deriv) return
       
      if (ndi.eq.3) then
        write(iwrit,*)'Error: global derivatives for 2D elements'
        write(iwrit,*)'in 3D not coded'
        call termin(8)
      endif

      sx(1,1) = xs(2,2)/xjaco
      sx(2,2) = xs(1,1)/xjaco
      sx(1,2) =-xs(1,2)/xjaco
      sx(2,1) =-xs(2,1)/xjaco

c     form global derivatives
      do i = 1,nno
        tp          = shape(1,i)*sx(1,1) + shape(2,i)*sx(2,1)
        shape(2,i)  = shape(1,i)*sx(1,2) + shape(2,i)*sx(2,2)
        shape(1,i)  = tp
      enddo
      return
      end

