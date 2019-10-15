c**********************************************************
      subroutine data(nF,argin)
c     read input file, allocate memory and fill variables
c     the input file is read multiple times to allow memory allocation
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer nF
      character*80 argin
      character*4 th! for the environment variable OMP_NUM_THREADS
c----------------------------------------------------------
      call printheader ! print header on screen


      call readinput1(nF,argin)  ! read control data
       
c     set the number of openmp threads
      if (nCores.eq.0) then! nCores has not been defined: read the OMP_NUM_THREADS environment variable
        th='    '
        call getenv("OMP_NUM_THREADS",th)
        if (th(1:1).eq.' ') then! OMP_NUM_THREADS has not been defined
          write(iwrit,*) "##  Source/lmain.f: the simulation will run",
     $          " on 1 core (nCores is unspecified in the input file)."
          nCores=1! use 1 as the default value
        else
          write(iwrit,*) " "
          write(iwrit,*) "##  Source/lmain.f: the simulation will run",
     $          " on ",th(1:1),
     $          " core(s) (nCores is unspecified in the input file,",
     $          " using OMP_NUM_THREAD environment variable)."
          read(th,'(I4)') nCores
        endif
      endif
      call omp_set_num_threads(nCores)

      call readinput2(nF,argin)  ! read more control data
       
      call readinput3(nF,argin)  ! read everything else 

c     make a few more checks
      if (problemType.eq.2) then
        if (nmeas.ne.1) then
          write(iwrit,*) 'Warning: nmeas should be 1 for forward ',
     &                   'problems; forcing nmeas=1'
          nmeas=1
        endif
      endif

      return
      end


c**********************************************************
      subroutine printheader
c     print header
      USE IOUNIT
      USE MAINMEM
      implicit none
      character*2 titll(44)
      character*24 cdate1
c----------------------------------------------------------
      call fdate(cdate1)
      
      write(iwrit,901) version,verdate,cdate1
      write(iwrit,902)
      write(iwrit,903)
      
      return
 901  format(12x,51('*')/12x,'*',49(' '),'*'/
     +       12x,'*',49(' '),'*'/
     +       12x,'*',22(' '),'NLACE',22(' '),'*'/
     +       12x,'*',49(' '),'*'/
     +       12x,'*',12(' '),'Non-Linear Adjoint-based',13(' '),'*'/
     +       12x,'*',13(' '),'Coefficient Estimation',14(' '),'*'/
     +       12x,'*',49(' '),'*'/
     +       12x,'*',12(' '),'version     : ',a5,18(' '),'*'/
     +       12x,'*',12(' '),'version date: ',a8,15(' '),'*'/
     +       12x,'*',49(' '),'*'/12x,51('*')//2x,a24)
 902  format(/,' For more screen printing regarding reading the ',
     +     'inputfile, check variable (datpr) in ldata.f'/
     +     ' For more screen printing regarding pardiso, check',
     +     ' variable (msglvl) in pardisosolve.f',/)
 903  format(' The outputed files are: *.vtk (to load in paraview),',
     +    ' *.res to restart the simulation,',/,'  *.ite and *.cvg',
     +    ' to check the convergence of the optimization algorithm',/)
      end


c**********************************************************
      subroutine readinput1(nF,argin)
c     first reading of the input file: get the problem's size
      USE IOUNIT
      USE MAINMEM
      implicit none
      character*4 task(9)
      character*80 dtask,argin
      data task/'femo','opto','solo','outp','opts','timi',
     &   'pbty','optp','end '/
      integer ii,nF
c-----------------------------------------------------------
c     initialize parameters to default values
c      datpr = .true.! default: will not print a copy of the input file
      datpr = .false.! default: will not print a copy of the input file
      timingpr=.false.! default: will print minimal timing info (may be modified in readinput1 and data16)
      monitorRun=.false.! default: no plot of the reconstruction as the simulation goes
      nset_nodal = 0
      nset_elemental = 0
      prefixOut(1:40)='                                        '
      prefixOut(41:80)='                                        '
      scaleOptVar=1.0d0
      problemType=1! default: will solve an inverse problem (may be modified in readinput1 and data17)
      optPrecondition=.false.! default use a diagonal preconditioning for the optimization based on the gradient norm
      optPnZones1=1! default
      optPnZones2=1! default
      buildSymmetricMatrices=.true.! default

      open (unit = iread,file=argin(1:nF),status = 'unknown')

 150  read (iread,'(a80)') dtask
      
      do ii = 1,10
        if (dtask(1:4).eq.task(ii)) go to 350
      enddo
      go to 150
      
 350  if (ii.eq.1) then
        call data11
      elseif (ii.eq.2) then
        call data12
      elseif (ii.eq.3) then
        call data13
      elseif (ii.eq.4) then
        call data14
      elseif (ii.eq.5) then
        call data15
      elseif (ii.eq.6) then
        call data16
      elseif (ii.eq.7) then
        call data17
      elseif (ii.eq.8) then
        call data18
      elseif (ii.eq.9) then
        close(iread,status='keep')
        return
      endif
      go to 150
      end


c**********************************************************
      subroutine data11
c     read control data for the solver
      USE IOUNIT
      USE MAINMEM
      implicit none
c-----------------------------------------------------------
c     read control parameters (all global variables, see macom.f)
      tol=1.0d0
      read(iread,*) nelem,npoin,ndime,mnode,mdofn,nmat,mprops,
     $    nmeas,mpoinbc,tol,lsteps,ncontin
      write(iwrit,*) ' '
      write(iwrit,*) 'number of elements ....................',nelem
      write(iwrit,*) 'number of points    ...................',npoin
      write(iwrit,*) 'number of spatial dimensions...........',ndime
      write(iwrit,*) 'max. number of nodes per elem..........',mnode
      write(iwrit,*) 'max. number of dof per node............',mdofn
      write(iwrit,*) 'number of material types...............',nmat
      write(iwrit,*) 'max. number of properties per material.',mprops
      write(iwrit,*) 'number of measurements.................',nmeas
      write(iwrit,*) 'max. num of nodes per bc set...........',mpoinbc
      write(iwrit,'(A,1p,E10.4)')
     $ ' tolerance for Newton iteration.........      ',tol
      write(iwrit,*) 'number of loadings.....................',lsteps
      write(iwrit,*) 'number of incremental material steps...',ncontin
      
      CALL MEMORY_ALLOCA
      return
      end


c**********************************************************
      subroutine data12
c     read control data for the optimization (BFGS)
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer status
c-----------------------------------------------------------
      read(iread,*) iopt,niter,bfgsM,noutput,mnty
      write(iwrit,*) ' '
      write(iwrit,*) 'optimization (1/2/3:sd/bfgs/newton)....',iopt
      write(iwrit,*) 'maximum number of iterations...........',niter
      write(iwrit,*) 'value of m parameter in BFGS...........',bfgsM
      write(iwrit,*) 'results saved every x iterations ... x=',noutput
      write(iwrit,*) 'maximum number of fields to optimize...',mnty
       
      ALLOCATE (lowerb(mnty), STAT=status)
      if (status.ne.0) then
        write(*,*)'memory failiure: lowerb'
        stop
      endif
      ALLOCATE (upperb(mnty), STAT=status)
      if (status.ne.0) then
        write(*,*)'memory failiure: upperb'
        stop
      endif

      return
      end


c**********************************************************
      subroutine data13
c     read the options for the solver
      USE IOUNIT
      USE MAINMEM
      implicit none
c-----------------------------------------------------------
      nCores=0! initialize nCores
      read(iread,*) nCores, solveMethod
      write(iwrit,*) ' '
      write(iwrit,*) 'number of Cores used by the solver ....',
     $    nCores
      if (solveMethod.eq.1) then
        write(iwrit,*) 'inverting matrices with a direct method:',
     $        ' pardiso'
      elseif (solveMethod.eq.2) then
        write(iwrit,*) 'inverting matrices with an iterative method:',
     $        ' RCI FGMRES + ILU0 preconditioning'
      elseif (solveMethod.eq.3) then
        write(iwrit,*) 'inverting matrices with an iterative method:',
     $        ' RCI FGMRES + ILUT preconditioning'
      else
        write(iwrit,*)'solveMethod=',solveMethod,' is not ok...'
        stop
      endif
      return
      end


c**********************************************************
      subroutine data14
c     read the prefix for the output files (if specfied)
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer nF2
c-----------------------------------------------------------
      read(iread,*) prefixOut
      nF2=1 ! determine the size of word in prefixOut
      do while (prefixOut(nF2:nF2).ne." ")
        nF2=nF2+1
      enddo
      nF2=nF2-1
      write(iwrit,*) ' '
      write(iwrit,*) 'prefix for the output file ............   ',
     $    prefixOut(1:nF2)
      return
      end


c**********************************************************
      subroutine data15
c     read the scaling factor (if specified)
      USE IOUNIT
      USE MAINMEM
      implicit none
c-----------------------------------------------------------
      read(iread,*) scaleOptVar
      write(iwrit,*) ' '
      write(iwrit,*) 'scaling factor for the opt. variables..   ',
     $    scaleOptVar
      end


c**********************************************************
      subroutine data16
c     set the logical/boolean for printing the timers info (for profiling)
      USE IOUNIT
      USE MAINMEM
      implicit none
c-----------------------------------------------------------
      timingpr=.true.
      write(iwrit,*) ' '
      write(iwrit,*) 'will print the timers info for profinling'
      end


c**********************************************************
      subroutine data17
c     set the integer that defines the type of problem to solve: forward or inverse
      USE IOUNIT
      USE MAINMEM
      implicit none
c-----------------------------------------------------------
      ! problemType=1 inverse problem (default)
      ! problemType=2 forward problem
      read(iread,*) problemType
      write(iwrit,*) ' '
      if (problemType.eq.1) then
        write(iwrit,*) 'Will solve an inverse problem: problemType=',
     &                 problemType
      elseif (problemTYpe.eq.2) then
        write(iwrit,*) 'Will solve a forward problem: problemType=',
     &                 problemType
      else
        write(iwrit,*) 'Error: value ',problemType,' is not a valid ',
     &     'entry for a problemType; exiting'
        stop
      endif
      end


c**********************************************************
      subroutine data18
c     read the scaling factor (if specified)
      USE IOUNIT
      USE MAINMEM
      implicit none
c-----------------------------------------------------------
      optPrecondition=.true.
      optPnZones1=4
      optPnZones2=4
      read(iread,*) optPrestart
      write(iwrit,*) ' '
      write(iwrit,*) 'preconditioning the opt with restart period= ',
     $    optPrestart
      end

c**********************************************************
      subroutine readinput2(nF,argin)
c     read the input file for the first time
c      to get the size of the variables
      USE IOUNIT
      USE MAINMEM
      implicit none
      character*4 task(4)
      character*80 dtask,argin
      data task/'elem','optb','datn','end '/
      integer ii,nF
c-----------------------------------------------------------------
      open (unit = iread,file=argin(1:nF),status = 'unknown')

 160  read (iread,'(a80)') dtask
      
      do ii = 1,4
        if (dtask(1:4).eq.task(ii)) go to 360
      enddo
      go to 160
      
 360  if (ii.eq.1) then
        call data21
      elseif (ii.eq.2) then
        call data22
      elseif (ii.eq.3) then
        call data23(iread)
      elseif (ii.eq.4) then
        close(iread,status='keep')
        return
      endif
      go to 160
      end


c**********************************************************
      subroutine data21
c     read element reference data
c     lrefn (ielem,1) = element type number
c     lrefn (ielem,2) = number of nodes
c     lrefn (ielem,3) = material set number

      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem,jelembeg,jelemend,jelemtype,jelemnnode
      integer imat
c-----------------------------------------------------------------
      do imat = 1,nmat
        read(iread,*) jelembeg,jelemend,jelemtype,jelemnnode
        do ielem = jelembeg,jelemend
          lrefn(ielem,1) = jelemtype
          lrefn(ielem,2) = jelemnnode
          lrefn(ielem,3) = imat
        enddo
c       check the input file format; stop the code if there is a pb
        if (jelemtype.eq.305) then
          if (mprops.le.9) then
            print*,'mprops is less or equal to 9 for elem305: exiting'
            stop
          endif
          if (ndime.le.2) then
            print*,'ndime is less of equal to 2 for elem305: exiting'
            stop
          endif
          if (mnode.le.3) then
            print*,'mnode is less of equal to 3 for elem305: exiting'
            stop
          endif
          if (mdofn.le.3) then
            print*,'mdofn is less of equal to 3 for elem305: exiting'
            stop
          endif
        elseif (jelemtype.eq.306) then
          if (ndime.le.2) then
            print*,'ndime is less of equal to 1 for elem306: exiting'
            stop
          endif
          if (mnode.le.3) then
            print*,'mnode is less of equal to 3 for elem306: exiting'
            stop
          endif
          if (mdofn.le.3) then
            print*,'mdofn is less of equal to 2 for elem306: exiting'
            stop
          endif
        elseif (jelemtype.eq.307) then
          if (ndime.le.2) then
            print*,'ndime is less of equal to 1 for elem307: exiting'
            stop
          endif
          if (mnode.le.3) then
            print*,'mnode is less of equal to 3 for elem307: exiting'
            stop
          endif
          if (mdofn.le.3) then
            print*,'mdofn is less of equal to 2 for elem307: exiting'
            stop
          endif
        elseif (jelemtype.eq.505) then
          if (mprops.le.9) then
            print*,'mprops is less or equal to 9 for elem505: exiting'
            stop
          endif
          if (ndime.le.1) then
            print*,'ndime is less of equal to 1 for elem505: exiting'
            stop
          endif
          if (mnode.le.2) then
            print*,'mnode is less of equal to 3 for elem505: exiting'
            stop
          endif
          if (mdofn.le.2) then
            print*,'mdofn is less of equal to 2 for elem505: exiting'
            stop
          endif
        elseif (jelemtype.eq.606) then
          if (mprops.le.7) then
            print*,'mprops is less or equal to 7 for elem606: exiting'
            stop
          endif
          if (ndime.le.1) then
            print*,'ndime is less of equal to 1 for elem606: exiting'
            stop
          endif
          if (mnode.le.3) then
            print*,'mnode is less of equal to 3 for elem606: exiting'
            stop
          endif
          if (mdofn.le.1) then
            print*,'mdofn is less of equal to 1 for elem606: exiting'
            stop
          endif
        elseif (jelemtype.eq.37) then
          if (mprops.le.5) then
            print*,'mprops is less or equal to 5 for elem37: exiting'
            stop
          endif
          if (ndime.le.2) then
            print*,'ndime is less of equal to 1 for elem37: exiting'
            stop
          endif
          if (mnode.le.3) then
            print*,'mnode is less of equal to 3 for elem37: exiting'
            stop
          endif
          if (mdofn.le.3) then
            print*,'mdofn is less of equal to 1 for elem37: exiting'
            stop
          endif
        elseif ((jelemtype.eq.607).or.(jelemtype.eq.608)) then
          if (mprops.le.8) then
            print*,'mprops is less or equal to 8 for elem607: exiting'
            stop
          endif
          if (ndime.le.1) then
            print*,'ndime is less of equal to 1 for elem607: exiting'
            stop
          endif
          if (mnode.le.3) then
            print*,'mnode is less of equal to 3 for elem607: exiting'
            stop
          endif
          if (mdofn.le.1) then
            print*,'mdofn is less of equal to 1 for elem607: exiting'
            stop
          endif
        elseif(jelemtype.eq.306) then
        elseif(jelemtype.eq.611) then
        elseif(jelemtype.eq.618) then
        elseif(jelemtype.eq.66) then
        elseif(jelemtype.eq.77) then
        elseif(jelemtype.eq.707) then
        elseif(jelemtype.eq.31) then
        elseif(jelemtype.eq.32) then
        elseif(jelemtype.eq.309) then
        else
          print*,'element number is wrong or not checked yet: exiting'
          stop
        endif
      enddo
      return
      end


c**********************************************************
      subroutine data22
c     read the variable's bounds for optimization
      USE IOUNIT
      USE MAINMEM
      implicit none
      double precision upper,lower
      integer ii
c-----------------------------------------------------------------
      do ii= 1,mnty
        read(iread,*) lower,upper
        lowerb(ii)=lower
        upperb(ii)=upper
      enddo
      return
      end


c*****************************************************************
      subroutine data23(idfile)
c     read nodal data.
c     read in a list of data specified at nodes.
c     can be used for anythng.
c     user will have access to this data within each element.
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ipoin,iset
      integer idfile
      integer status
c-----------------------------------------------------------------
      read(idfile,*)nset_nodal
c     print heading
      if (datpr) then
        write(iwrit,800)
        if      (nset_nodal.eq. 1) then
          write(iwrit,901)
        else if (nset_nodal.eq. 2) then
          write(iwrit,902)
        else if (nset_nodal.eq. 3) then
          write(iwrit,903)
        else if (nset_nodal.eq. 4) then
          write(iwrit,904)
        else if (nset_nodal.eq. 5) then
          write(iwrit,905)
        else
          write(iwrit,998) nset_nodal
        endif
      endif

c     allocate memory to store nodal data (no reallocation with restart)
      if (idfile.eq.iread) then
        call MEMORY_ALLOCE
      endif

c     read nodal data
      noptvar = 0
      do ipoin = 1,npoin
        read(idfile,*)(ipoin_nod_grad(iset,ipoin),
     $       adata_nodal(iset,ipoin),
     $       iset=1,nset_nodal)
        if (datpr) write(iwrit,951) ipoin,
     $      (adata_nodal(iset,ipoin),
     $      ipoin_nod_grad(iset,ipoin),iset=1,nset_nodal)
        do iset = 1,nset_nodal
          if (ipoin_nod_grad(iset,ipoin).ne.0) then
            noptvar = noptvar +1
            ipoin_nod_grad(iset,ipoin) = noptvar
          endif
        enddo
      enddo
       
      adata_nodal2(:,:)=adata_nodal(:,:) ! initialize adata_nodal2 properly
      if (idfile.eq.iread) then
        ALLOCATE (xmem(noptvar), STAT=status)
        if (status.ne.0) then
           write(*,*)'memory failiure: xmem: status = ',status
           stop
        endif
        xmem(:) = 1.0d8

        ALLOCATE (gmem(noptvar), STAT=status)
        if (status.ne.0) then
           write(*,*)'memory failiure: gmem: status = ',status
           stop
        endif
        gmem(:) = 1.0d8

        print*,"nfStored=",nfStored

        ALLOCATE (fmem(nfStored), STAT=status)
        if (status.ne.0) then
           write(*,*)'memory failiure: fmem: status = ',status
           stop
        endif
        fmem(:) = 1.0d8
      endif

       
      if (datpr) then
        write(iwrit,*)'number of optimization variables=',noptvar
      end if

      return

 800  format(//2x,'nodal data'/2x,32('-'))
 901  format(2x,'node',3x,'value')
 902  format(2x,'node',2(2x,'value   '))
 903  format(2x,'node',3(2x,'value   '))
 904  format(2x,'node',4(2x,'value   '))
 905  format(2x,'node',5(2x,'value   '))
 951  format(i5,1x,7(f9.3,2x,i5,2x))
 998  format(/2x,
     +'                  proper heading in data5 is not available'/2x,
     +'                  for nset_nodal = ',i2,'.'/)
      end


c**********************************************************
      subroutine readinput3(nF,argin)
c     read problem data groups
      USE IOUNIT
      USE MAINMEM
      implicit none
      character*4 task(11)
      character*80 dtask,argin
      data task/'coor','conn','mate','boun','rest','date','meas',
     $     'fmnd','bf','weak','end '/
      integer ii,nF
c      integer i
c      double precision mincr,maxcr
c      double precision cr(nelem)
c-----------------------------------------------------------------
      open (unit = iread,file=argin(1:nF),status = 'unknown')

 100  read (iread,'(a80)') dtask
      
      do ii = 1,11
        if (dtask(1:4).eq.task(ii)) go to 300
      enddo
      go to 100
      
 300  if (ii.eq.1) then
        call data31 
      elseif (ii.eq.2) then
        call data32 
      elseif (ii.eq.3) then
        call data33
      elseif (ii.eq.4) then
        call data34 
      elseif (ii.eq.5) then
c       restart overwrites the data read with keyword datn
        call data35 
      elseif (ii.eq.6) then
        call data36 
      elseif (ii.eq.7) then
        call data37 
      elseif (ii.eq.8) then
        call data38 
      elseif (ii.eq.9) then
        call data39 
      elseif (ii.eq.10) then
        call data399 
      elseif (ii.eq.11) then
c        write(iwrit,901)
        close(iread,status='keep')
c        call meshQuality2Dtris(5,lnods(1:5,:),npoin,coord,cr)
c        mincr=cr(1)
c        maxcr=cr(1)
c        do i=2,nelem
c          if (cr(i).gt.maxcr) maxcr=cr(i)
c          if (cr(i).lt.mincr) mincr=cr(i)
c        enddo
c        print*,"min quality criterion = ",mincr
c        print*,"max quality criterion = ",maxcr
c        stop
        return
      endif
      go to 100
 901  format(/5x,22('-')'end of third reading',22('-'))
      end


c***********************************************************
      subroutine data31 
c     read nodal coordinates
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ipoin,idime,jpoin
c-----------------------------------------------------------
c     print heading
      if (datpr) then
        write(iwrit,900)
        if (ndime.eq.1) write(iwrit,901)
        if (ndime.eq.2) write(iwrit,902)
        if (ndime.eq.3) write(iwrit,903)
      end if

c     read data
      do ipoin = 1,npoin
        read(iread,*)jpoin,(tvect(idime),idime = 1,ndime)
        coord(1:ndime,jpoin) = tvect(1:ndime)
      enddo
      ! reinitialize tvect for fgmres+ilu0 and pardiso solvers
      tvect(:)=0.0d0
      tvect(1)=0.0001d0

c     print coordinates
      if (datpr) then
        do 10 ipoin = 1,npoin
          write(iwrit,910) ipoin,(coord(idime,ipoin),idime = 1,ndime)
10      continue
      end if

      return

900   format(//2x,'nodal coordinates'/2x,17('-'))
901   format(2x,'node',5x,'1 - coord')
902   format(2x,'node',5x,'1 - coord',8x,'2 - coord')
903   format(2x,'node',5x,'1 - coord',8x,'2 - coord',8x,'3 - coord')
910   format(i5,3(2x,e15.5))
      end


c******************************************************************
      subroutine data32 
c     (i)  read element topology
c     (ii) set the dof per node value for each node
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem,jelem,ipoin,inode,nn
      integer nno
      double precision pelem(mprops)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,mnode)
      double precision elemdata_nodal(nset_nodal,mnode)
      double precision uelem(mdofn,mnode)
      double precision uelem_meas(mdofn,mnode)
c-------------------------------------------------------------------
      nno = 0
      pelem(:) = 0.0d0
      estif(:,:) = 0.0d0
      eforc(:) = 0.0d0
      elnods(:) = 0
      xelem(:,:) = 0.0d0
      elemdata_nodal(:,:) = 0.0d0
      uelem(:,:) = 0.0d0
      uelem_meas(:,:) = 0.0d0
      
      if (datpr) write(iwrit,900)
       
      do ielem = 1,nelem
c       read and print element connectivity
        read(iread,*)jelem,
     $       (lnods(jelem,inode),inode = 1,lrefn(jelem,2))
        if (datpr)  then 
          write(iwrit,902) jelem,
     +          (lnods(jelem,inode),inode=1,lrefn(jelem,2))
        endif

c       set the dof per node for each node in the element
c       the information is in the element and is retrieved via
c       the call elmlib(ielem,1,...). 
c       note: if there is a discrepancy in dof/node from one element
c       to another, it will not be flagged. instead the dof/node 
c       value from the element with the higher number will be retained.

        nno = lrefn(jelem,2)
        call elmlib (jelem,1,pelem,nno,estif,eforc,elnods,xelem,
     $         elemdata_nodal,uelem,uelem_meas)
          
        do inode = 1,nno
          nn = lnods(jelem,inode)
          ldofn(nn) = elemvec_ndofn(inode)
        enddo
          
      enddo

      if (datpr) then
        write(iwrit,950)
        do ipoin = 1,npoin
          write(iwrit,951) ipoin,ldofn(ipoin)
        enddo
      end if

      return
900   format(//2x,'element topology'/2x,16('-')
     +        /2x,'element',7x,'global node numbers')
902   format(i6,7x,100(9(i5,2x)/13x))

950   format(//2x,'number of dof at nodes'/2x,22('-')
     +         /2x,' node   ',7x,'number of dof')
951   format(i6,11x,i8)
      end


c***********************************************************
      subroutine data33
c     read material properties
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer imat,ielem,iprop,imatp
      integer nno
      double precision pelem(mprops)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,mnode)
      double precision elemdata_nodal(nset_nodal,mnode)
      double precision uelem(mdofn,mnode)
      double precision uelem_meas(mdofn,mnode)
c-----------------------------------------------------------
      imatp = 0
      nno = 0
      pelem(:) = 0.0d0
      estif(:,:) = 0.0d0
      eforc(:) = 0.0d0
      elnods(:) = 0
      xelem(:,:) = 0.0d0
      elemdata_nodal(:,:) = 0.0d0
      uelem(:,:) = 0.0d0
      uelem_meas(:,:) = 0.0d0

      do ielem = 1,nelem
        imat = lrefn(ielem,3)
        if (imatp.ne.imat) then
          if (datpr) then
            write(iwrit,900)
          endif
          call elmlib(ielem,2,pelem,nno,estif,eforc,elnods,xelem,
     $              elemdata_nodal,uelem,uelem_meas)
c openmp             lprop(imat) = nprop! lprop was removed, now nprop=mprops for all element
          do iprop = 1,mprops! lprop(imat)
            props(imat,iprop) = pelem(iprop)
          enddo
          imatp = imat
        endif
      enddo

      if (datpr) then
        do imat = 1,nmat
          write(iwrit,900)
          write(iwrit,901)imat,mprops! lprop(imat)
          do iprop = 1,mprops! lprop(imat)
            write(iwrit,902)iprop,props(imat,iprop)
          enddo
        enddo
      end if
      return

 900  format(/,' element properties: ')
 901  format('material set number',i5,3x,'number of props',i5)
 902  format('property number    ',i5,3x,'property value ',e15.6)

      end


c*****************************************************************
      subroutine data34 
c     read essential (Dirichlet) boundary conditions
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ipoin,nodfx,idofn,imeas,npoin_bc
c-----------------------------------------------------------------
c     print heading
      if (datpr) then
        write(iwrit,800)
        if      (mdofn.eq. 1) then
          write(iwrit,901)
        else if (mdofn.eq. 2) then
          write(iwrit,902)
        else if (mdofn.eq. 3) then
          write(iwrit,903)
        else if (mdofn.eq. 4) then
          write(iwrit,904)
        else if (mdofn.eq. 5) then
          write(iwrit,905)
        else if (mdofn.eq. 6) then
          write(iwrit,906)
        else if (mdofn.eq. 7) then
          write(iwrit,907)
        else if (mdofn.eq. 9) then
          write(iwrit,909)
        else if (mdofn.eq.11) then
          write(iwrit,911)
        else if (mdofn.eq.12) then
          write(iwrit,912)
        else if (mdofn.eq.15) then
          write(iwrit,915)
        else if (mdofn.eq.19) then
          write(iwrit,919)
        else if (mdofn.eq.23) then
          write(iwrit,923)
        else
          write(iwrit,998) mdofn
        endif
      end if

c     read boundary condition data
      do imeas = 1,nmeas
        read(iread,*) npoin_bc
        if (mpoinbc.lt.npoin_bc) then! warning message for pb in input file
          print*,'Warning there are more points in the boun* section',
     $       ' than the number declared in the femo* section'
          print*,'Exiting'
          stop
        endif
        do ipoin = 1,npoin_bc
          read(iread,*)nodfx
          ipoin_bc(imeas,ipoin) = nodfx
          read(iread,*)(icode_bc(imeas,ipoin,idofn),
     $          value_bc(imeas,ipoin,idofn),idofn=1,ldofn(nodfx))
             
          if (datpr) write(iwrit,951) nodfx,(
     $         icode_bc(imeas,ipoin,idofn),
     $         value_bc(imeas,ipoin,idofn),
     $         idofn=1,ldofn(nodfx))
        enddo
      enddo

800   format(//2x,'displacement boundary conditions'/2x,32('-'))
901   format(2x,'node',2x,'code',3x,'value')
902   format(2x,'node',2(2x,'code   value   '))
903   format(2x,'node',3(2x,'code   value   '))
904   format(2x,'node',4(2x,'code   value   '))
905   format(2x,'node',5(2x,'code   value   '))
906   format(2x,'node',6(2x,'code   value   '))
907   format(2x,'node',7(2x,'code   value   '))
909   format(2x,'node',7(2x,'code   value   ')/6x,
     +                 2(2x,'code   value   '))
911   format(2x,'node',7(2x,'code   value   ')/6x,
     +                 4(2x,'code   value   '))
912   format(2x,'node',7(2x,'code   value   ')/6x,
     +                 5(2x,'code   value   '))
915   format(2x,'node',7(2x,'code   value   ')/6x,
     +       7(2x,'code   value   ')/8x,'code   value   ')
919   format(2x,'node',7(2x,'code   value   ')/6x,
     +       7(2x,'code   value   ')/6x,5(2x,'code   value   '))
923   format(2x,'node',7(2x,'code   value   ')/6x,
     +       7(2x,'code   value   ')/6x,7(2x,'code   value   '),
     +      /6x,2(2x,'code   value   '))
951   format(i5,1x,7(i5,1x,f9.3,2x)/,6x,7(i5,1x,f9.3,2x)/,
     +          6x,7(i5,1x,f9.3,2x)/,6x,2(i5,1x,f9.3,2x))
998   format(/2x,
     +'                  proper heading in data4 is not available'/2x,
     +'                  for mdofn = ',i2,'.'/)

      end


c*****************************************************************
      subroutine data35
c     restart for datn
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer nF
      character*80 namef
c-----------------------------------------------------------------
c     read the name of the input file
      namef(1:40)='                                       '
      namef(41:80)='                                       '
      read(iread,*)namef
c     get the number of letters in the name
      nF=1 ! number of non-empty characters
      do while (namef(nF:nF).ne." ")
        nF=nF+1
      enddo
      nF=nF-1
      write(iwrit,'(/,A,A)') ' Restarting with file .... ',namef(1:nF)
      open(unit=4,file=namef(1:nF),status='unknown')
c     overwrite the data read in datn
      call data23(4)
      close(4,status='keep')
      return
      end


c*****************************************************************
      subroutine data36
c     read element data.
c     read in a list of data specified at elements.
c     can be used for anything.
c     user will have access to this data within each element.
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem,iset
c-----------------------------------------------------------------
      read(iread,*)nset_elemental

c     print heading
      if (datpr) then
        write(iwrit,800)
        if      (nset_elemental.eq. 1) then
          write(iwrit,901)
        else if (nset_elemental.eq. 2) then
          write(iwrit,902)
        else if (nset_elemental.eq. 3) then
          write(iwrit,903)
        else if (nset_elemental.eq. 4) then
          write(iwrit,904)
        else if (nset_elemental.eq. 5) then
          write(iwrit,905)
        else
          write(iwrit,998) nset_elemental
        endif
      end if

c     allocate memory to store nodal data
      call MEMORY_ALLOCF

c     read nodal data
      do ielem = 1,nelem
        read(iread,*)(adata_elemental(iset,ielem),iset=1,
     $       nset_elemental)
        if (datpr) write(iwrit,951) ielem,
     $            (adata_elemental(iset,ielem),iset=1,nset_elemental)
      enddo

      return

 800  format(//2x,'elemental data'/2x,32('-'))
 901  format(2x,'element',3x,'value')
 902  format(2x,'element',2(2x,'value   '))
 903  format(2x,'element',3(2x,'value   '))
 904  format(2x,'element',4(2x,'value   '))
 905  format(2x,'element',5(2x,'value   '))
 951  format(i5,1x,7(f9.3,2x)/,6x,7(f9.3,2x)/,
     +          6x,7(f9.3,2x)/,6x,2(f9.3,2x))
 998  format(/2x,
     +'                  proper heading in data6 is not available'/2x,
     +'                  for nset_elemental = ',i2,'.'/)
      end


c*****************************************************************
      subroutine data37
c     read measured nodal data.
c     read in a list of measured data specified at nodes.
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer imeas,ipoin,idofn
c-----------------------------------------------------------------
      do imeas = 1,nmeas
        do ipoin = 1,npoin
          read(iread,*)(
     $          fac_meas(imeas,ipoin,idofn),
     $          meas(imeas,ipoin,idofn),
     $          idofn = 1,mdofn)
        enddo
      enddo

      end

c*****************************************************************
      subroutine data38
c     read non-diagonal components of fac_meas
c     read in a list specified at nodes.
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer imeas,ipoin,idofn
c-----------------------------------------------------------------
      do imeas = 1,nmeas
        do ipoin = 1,npoin
          read(iread,*)(
     $          fac_meas_nd(imeas,ipoin,idofn),
c     $          idofn = 1,3)
     $          idofn = 1,mdofn*(mdofn-1)/2)
        enddo
      enddo

      end

c*****************************************************************
      subroutine data39
c     read in body force values 
c     read in a list specified at nodes.
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer imeas,ipoin,idofn
c-----------------------------------------------------------------
        do ipoin = 1,npoin
          read(iread,*)(
     $          bodyforce(ipoin,idofn),
     $          idofn = 1,2)
        enddo
c        write(123,*) bodyforce
      end


c*****************************************************************
      subroutine data399 
c     read weak spring data
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ipoin,nodfx,idofn,imeas,npoin_ws
c-----------------------------------------------------------------

c     read weak spring data
      do imeas = 1,nmeas
        read(iread,*) npoin_ws
c        write(123,*) npoin_ws 
        read(iread,*) delta_ws
c        write(456,*) delta_ws
        do ipoin = 1,npoin_ws
          read(iread,*)nodfx
          ipoin_ws(imeas,ipoin) = nodfx
          read(iread,*)(icode_ws(imeas,ipoin,idofn),
     $          value_ws(imeas,ipoin,idofn),idofn=1,ldofn(nodfx))
             
        enddo
      enddo
c       write(222,*) ipoin_ws
c       write(333,*) value_ws
c       write(444,*) icode_ws

      end


