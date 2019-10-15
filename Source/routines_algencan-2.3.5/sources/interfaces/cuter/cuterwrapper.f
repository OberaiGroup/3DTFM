C     =================================================================
C     File: cuterwrapper.f
C     =================================================================

C     =================================================================
C     Module: CUTEr Wrapper
C     =================================================================

C     Last update of any component of this module:

C     September 11, 2009.

C     Users are encouraged to download periodically updated versions of
C     this code at the TANGO home page:

C     www.ime.usp.br/~egbirgin/tango/

C     ******************************************************************
C     ******************************************************************

      subroutine inip(n,x,l,u,m,lambda,equatn,linear,coded,checkder)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer m,n

C     ARRAY ARGUMENTS
      logical coded(10),equatn(*),linear(*)
      double precision l(*),lambda(*),u(*),x(*)

#include "../../algencan/dim.par"
#include "cuter.com"

C     LOCAL SCALARS
      integer cind,i,sind
      character * 10 pname
      logical efirst,lfirst,nvfirst
      double precision dum

C     LOCAL ARRAYS
      character * 10 gnames(mmax),xnames(nmax)
      double precision c(mmax),cl(mmax),cu(mmax)

C     EXTERNAL SUBROUTINES
      external csetup,cnames

C     Set problem data

      open(10,file='OUTSDIF.d',form='formatted',status='old')

      rewind 10

      efirst  =  .true.
      lfirst  = .false.
      nvfirst = .false.

      call csetup(10,20,n,m,x,l,u,nmax,equatn,linear,lambda,cl,cu,
     +mmax,efirst,lfirst,nvfirst)

      call cnames(n,m,pname,xnames,gnames)

      close(10)

      write(*,*) 'Solving problem: ',pname

C     There are two options for two-sides inequality constraints:
C
C     (a) Add a slack variable and transform it into an equality
C         constraints plus a bound constraint for the slack.
C
C     (b) Split it into two one-side inequality constraints.
C
C     There is a local variable of this subroutine called useslacks.
C     If useslacks is TRUE then we implement (a) else, we implement (b).

      useslacks = .false.

C     Compute constraints to initialize slacks

      if ( useslacks ) then
          call cfn(n,m,x,dum,mmax,c)
      end if

      nranges = 0

      do i = 1,m

C         Equality constraint.

          if ( equatn(i) ) then

              slaind(i) =     - 1
              ccor(i)   =       0
              cmap(i)   =       i
              ca(i)     =   1.0d0
              cb(i)     = - cu(i)

C         Ranged inequality constraint: add slack or split it.

          else if ( cl(i) .gt. - 1.0d+20 .and. cu(i) .lt. 1.0d+20 ) then

              nranges = nranges + 1

C             Replace by c(x) - s = 0 and l <= s <= u.

              if ( useslacks ) then

                  sind = n + nranges

                  l(sind) = cl(i)
                  u(sind) = cu(i)
                  x(sind) = max( cl(i), min( c(i), cu(i) ) )

                  slaind(i) =    sind
                  ccor(i)   =       0
                  cmap(i)   =       i
                  ca(i)     =   1.0d0
                  cb(i)     =   0.0d0

                  equatn(i) = .true.

C             Split into c(x) - u <= 0 and l - c(x) <= 0.

              else
                  cind = m + nranges

                  equatn(cind) = equatn(i)
                  linear(cind) = linear(i)

                  slaind(cind) =     - 1
                  cmap(cind)   =       i
                  ca(cind)     =   1.0d0
                  cb(cind)     = - cu(i)

                  slaind(i)    =     - 1
                  ccor(i)      =    cind
                  cmap(i)      =       i
                  ca(i)        = - 1.0d0
                  cb(i)        =   cl(i)
              end if

C         Inequality constraint of type c(x) <= u.

          else if ( cu(i) .lt. 1.0d+20 ) then

              slaind(i) =     - 1
              ccor(i)   =       0
              cmap(i)   =       i
              ca(i)     =   1.0d0
              cb(i)     = - cu(i)

C         Inequality constraint of type l <= c(x).

          else if ( cl(i) .gt. - 1.0d+20 ) then

              slaind(i) =     - 1
              ccor(i)   =       0
              cmap(i)   =       i
              ca(i)     = - 1.0d0
              cb(i)     =   cl(i)

          end if

      end do

      mcuter = m
      ncuter = n

      if ( useslacks ) then
          n = n + nranges
      else
          m = m + nranges
      end if

C     Lagrange multipliers approximation

      do i = 1,m
          lambda(i) =  0.0d0
      end do

C     In this CUTEr interface evalfc, evalgjac, evalhl and evalhlp are
C     present. evalf, evalg, evalh, evalc, evaljac and evalhc are not.

      coded( 1) = .false. ! evalf
      coded( 2) = .false. ! evalg
      coded( 3) = .false. ! evalh
      coded( 4) = .false. ! evalc
      coded( 5) = .false. ! evaljac
      coded( 6) = .false. ! evalhc
      coded( 7) = .true.  ! evalfc
      coded( 8) = .true.  ! evalgjac
      coded( 9) = .true.  ! evalhl
      coded(10) = .false. ! evalhlp (coded but not being used)
                          ! (In fact, the case sf=0.0d0 is incomplete,
                          ! September 11th, 2009)

      checkder = .false.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine endp(n,x,l,u,m,lambda,equatn,linear)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalf(n,x,f,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

      flag = - 1

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hnnz,n

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalc(n,x,ind,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,jcnnz,n

C     ARRAY ARGUMENTS
      integer jcvar(*)
      double precision jcval(*),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhc(n,x,hclin,hccol,hcval,hcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hcnnz,n

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hnnz,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

#include "../../algencan/dim.par"
#include "cuter.com"

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision v(mmax)

C     EXTERNAL SUBROUTINES
      external csh,csh1

      flag = 0

      if ( sf .eq. 0.0d0 ) then
C         flag = - 1
C         return

C         Special case for sf=0
C         It is too inneficient !!!
C         call ievalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,flag)
C         return

          do i = 1,mcuter
              v(i) = 0.0d0
          end do

          do i = 1,m
              v(cmap(i)) = v(cmap(i)) + ca(i) * lambda(i) * sc(i)
          end do

          call csh1(ncuter,mcuter,x,mmax,v,hnnz,hnnzmax,hval,hlin,hcol)

C         Interchange row and column indices

          do i = 1,hnnz
              call intswap(hlin(i),hcol(i))
          end do

          return
      end if

      do i = 1,mcuter
          v(i) = 0.0d0
      end do

      do i = 1,m
          v(cmap(i)) = v(cmap(i)) + ca(i) * lambda(i) * sc(i) / sf
      end do

      call csh(ncuter,mcuter,x,mmax,v,hnnz,hnnzmax,hval,hlin,hcol)

C     Interchange row and column indices

      do i = 1,hnnz
          call intswap(hlin(i),hcol(i))
          hval(i) = hval(i) * sf
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hnnz,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

#include "../../algencan/dim.par"
#include "cuter.com"

C     LOCAL SCALARS
      integer col,con,hnnztmp,hnnzmaxtmp,i,ind,itmp,j,lin,nextj,rnnz
      double precision val

C     LOCAL ARRAYS
      integer hcon(hnnzmax),pos(nmax),rind(nmax),stlin(nmax)
      double precision rval(nmax),v(mmax)

C     EXTERNAL SUBROUTINES
      external cish

      flag = 0

C     ==================================================================
C     COMPUTE MULTIPLIERS FOR CUTER-STYLE FUNCTIONS
C     ==================================================================

      do i = 1,mcuter
          v(i) = 0.0d0
      end do

      do i = 1,m
          v(cmap(i)) = v(cmap(i)) + ca(i) * lambda(i) * sc(i)
      end do

C     ==================================================================
C     COMPUTE HESSIANS
C     ==================================================================

C     COMPUTE HESSIAN OF THE OBJECTIVE FUNCTION

      call cish(ncuter,x,0,hnnz,hnnzmax,hval,hlin,hcol)

C     For each element of the Hessian of the objective function,
C     set constraint index as zero
      do i = 1,hnnz
          hval(i) = hval(i) * sf
          hcon(i) = 0
      end do

C     COMPUTE HESSIANS OF THE CONSTRAINTS

      ind = 0
      hnnzmaxtmp = hnnzmax - hnnz

      do j = 1,mcuter
C         Compute Hessian of constraint j
          call cish(ncuter,x,j,hnnztmp,hnnzmaxtmp,hval(hnnz+ind+1),
     +    hlin(hnnz+ind+1),hcol(hnnz+ind+1))

C         For each element of the Hessian, set constraint as j
          do i = hnnz + ind + 1,hnnz + ind + hnnztmp
              hcon(i) = j
          end do

          ind = ind + hnnztmp
          hnnzmaxtmp = hnnzmaxtmp - hnnztmp
      end do

      hnnz = hnnz + ind

C     INTERCHANGE ROW AND COLUMN INDICES

      do i = 1,hnnz
          call intswap(hlin(i),hcol(i))
      end do

      if ( ind .eq. 0 ) return

C     ==================================================================
C     SET ROW LINKED LISTS
C     ==================================================================

C     Initialize pointers to the first element of each row
      do i = 1,n
         stlin(i) = 0
      end do

C     Set row linked lists
      do i = 1,hnnz
         lin        = hlin(i)
         itmp       = stlin(lin)
         stlin(lin) = i
         hlin(i)    = itmp
      end do

C     ==================================================================
C     BUILD HESSIAN OF THE LAGRANGIAN ROW BY ROW
C     ==================================================================

C     Initialize array pos
      do i = 1,n
          pos(i) = 0
      end do

      do i = 1,n
C         Initialize the i-th row of the Hessian of the Lagrangian
          rnnz = 0

C         Process the i-th row of all the Hessians
          j = stlin(i)

 10       if ( j .ne. 0 ) then

C             Process element (i,hcol(j)) of the Hessian of constraint
C             hcon(j) (Hessian of the objective function if hcon(j)=0)

              col = hcol(j)
              con = hcon(j)
              if ( con .eq. 0 ) then
                  val = hval(j)
              else
                  val = hval(j) * v(con)
              end if

              if ( pos(col) .ne. 0 ) then
                  rval(pos(col)) = rval(pos(col)) + val

              else
                  rnnz           = rnnz + 1
                  pos(col)       = rnnz
                  rind(pos(col)) = col
                  rval(pos(col)) = val
              end if

C             Get next element in the i-th row linked list
              j = hlin(j)
              go to 10
          end if

C         Clean array pos
          do j = 1,rnnz
              pos(rind(j)) = 0
          end do

C         Set i-th row of hl (over the i-th rows of the Hessians)
C         and mark remaining elements to be deleted
          j = stlin(i)

 20       if ( j .ne. 0 ) then
              nextj = hlin(j)

              if ( rnnz .ne. 0 ) then
                  hlin(j) = i
                  hcol(j) = rind(rnnz)
                  hval(j) = rval(rnnz)
                  rnnz    = rnnz - 1
              else
                  hlin(j) = 0
              end if

              j = nextj
              go to 20
          end if

      end do

C     Eliminate remaining elements (marked with hlin(j)=0)
      j = 1

 30   if ( j .le. hnnz ) then
          if ( hlin(j) .eq. 0 ) then
              if ( j .ne. hnnz ) then
                  hlin(j) = hlin(hnnz)
                  hcol(j) = hcol(hnnz)
                  hval(j) = hval(hnnz)
              end if
              hnnz = hnnz - 1
          else
              j = j + 1
          end if

          go to 30
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

#include "../../algencan/dim.par"
#include "cuter.com"

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision v(mmax)

C     EXTERNAL SUBROUTINES
      external cprod

      flag = 0

      if ( sf .eq. 0.0d0 ) then
          flag = - 1
          return

C         Special case for sf=0
C         It is too inneficient !!!
C         call ievalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)
C         return
      end if

      do i = 1,mcuter
          v(i) = 0.0d0
      end do

      do i = 1,m
          v(cmap(i)) = v(cmap(i)) + ca(i) * lambda(i) * sc(i) / sf
      end do

      call cprod(ncuter,mcuter,goth,x,mmax,v,p,hp)

      do i = 1,ncuter
          hp(i) = hp(i) * sf
      end do

      do i = ncuter + 1,n
          hp(i) = 0.0d0
      end do

      goth = .true.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

#include "../../algencan/dim.par"
#include "../../algencan/hessdat.com"
#include "cuter.com"

C     LOCAL SCALARS
      integer col,i,lin
      double precision val

      flag = 0

      if ( .not. goth ) then
          goth = .true.
          call ievalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,flag)
      end if

      do i = 1,ncuter
          hp(i) = 0.0d0
      end do

      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          val = hval(i)

          hp(lin) = hp(lin) + p(col) * val

          if ( lin .ne. col ) then
              hp(col) = hp(col) + p(lin) * val
          end if
      end do

      do i = ncuter + 1,n
          hp(i) = 0.0d0
      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalfc(n,x,f,m,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

#include "../../algencan/dim.par"
#include "cuter.com"

C     LOCAL SCALARS
      integer i,sind

C     EXTERNAL SUBROUTINES
      external cfn

      flag = 0

      call cfn(ncuter,mcuter,x,f,mmax,c)

      do i = m,1,-1
          c(i) = ca(i) * c(cmap(i)) + cb(i)
      end do

      do i = 1,m
          sind = slaind(i)
          if ( sind .ne. - 1 ) then
              c(i) = c(i) - x(sind)
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

#include "../../algencan/dim.par"
#include "cuter.com"

C     LOCAL SCALARS
      integer dum2,fu2,fun,i,jcnnztmp,sind,var

C     LOCAL ARRAYS
      double precision dum3(1)

C     EXTERNAL SUBROUTINES
      external csgr

      flag = 0

      call csgr(ncuter,mcuter,.false.,dum2,dum3,x,jcnnz,jcnnzmax,jcval,
     +jcvar,jcfun)

C     Remove gradient from the sparse structure

      do i = 1,n
          g(i) = 0.0d0
      end do

      i = 1

 100  if ( i .le. jcnnz ) then

          fun = jcfun(i)
          var = jcvar(i)

          if ( fun .eq. 0 ) then
              g(var) = g(var) + jcval(i)

              if ( i .ne. jcnnz ) then
                  call intswap(jcfun(jcnnz),jcfun(i))
                  call intswap(jcvar(jcnnz),jcvar(i))
                  call dblswap(jcval(jcnnz),jcval(i))
              end if

              jcnnz = jcnnz - 1

          else
              i = i + 1
          end if

          go to 100
      end if

C     Duplicate entrances correspoding to splitted constraints

      jcnnztmp = jcnnz

      do i = 1,jcnnztmp

          fun = jcfun(i)
          fu2 = ccor(fun)

          if ( fu2 .ne. 0 ) then
              jcnnz = jcnnz + 1
              jcfun(jcnnz) = fu2
              jcvar(jcnnz) = jcvar(i)
              jcval(jcnnz) = ca(fu2) * jcval(i)
          end if

          jcval(i) = ca(fun) * jcval(i)

      end do

C     Add effect of slack variables

      do i = 1,m
          sind = slaind(i)
          if ( sind .ne. - 1 ) then
              jcnnz = jcnnz + 1
              jcfun(jcnnz) = i
              jcvar(jcnnz) = sind
              jcval(jcnnz) = - 1.0d0
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine intswap(i,j)

      implicit none

C     SCALAR ARGUMENTS
      integer i,j

C     LOCAL SCALARS
      integer tmp

      tmp = i
      i   = j
      j   = tmp

      end

C     ******************************************************************
C     ******************************************************************

      subroutine dblswap(a,b)

      implicit none

C     SCALAR ARGUMENTS
      double precision a,b

C     LOCAL SCALARS
      double precision tmp

      tmp = a
      a   = b
      b   = tmp

      end

