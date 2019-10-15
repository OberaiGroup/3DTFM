      MODULE MAINMEM
      SAVE
c     --------------------------------------------
c                include file "macom.f"
c      --------------------------------------------
c     DECLARATIONS
c     integers 
      integer mprops ,mnode, mdofn, nmat
      integer mevab,hh 
      integer :: countFGcalls = 0
      integer npoin,nelem,nmats,ndime
      integer neqns,nnz, nnzconn, lsteps
      integer nset_nodal,nset_elemental
      integer noptvar,nmeas,mpoinbc,ncontin
      integer niter,bfgsM,noutput,mnty,iopt ! optimization options
      integer nCores ! solver options
      integer problemType! know if it is a forward or an inverse solve
      integer optPrestart! period (in number of iterations) between restart of L-BFGS-B
      integer optPnZones1,optPnZones2! how to subdivide the domain for diagonal preconditioning
      integer solveMethod! use either direct or iterative solver
      integer :: nextOutput = 0! for ouputing data
      integer :: currentOptIt = 0! for counting the iteration of the optimization algorithm
      integer :: lastOptIt = 0! temporary variable for counting the it. of the opt. alg.
      integer :: nfStored = 10! number of objective function values stored for setting the stopping criterion
      integer :: ifield = 1 ; ! current measurement field, for element 611 

c     logical
      logical negJac
      logical datpr, timingpr, optPrecondition, monitorRun
      logical newIt! for printing the *.cvg file
      logical buildSymmetricMatrices! build symmetric matrices and selected a solver for symmetric matrices
c     character
      character(6) :: version = '1.01 a'
      character(8) :: verdate = 'Feb 2006'
      character(80) prefixOut
c     double precision
      double precision regularization,dataMatch,forceMatch
      double precision convRes,convSol, tol, forceM, Gxi, Gforce
      double precision scaleOptVar, pforce, moment
      double precision ForceTotal,MomentTotal,MuFactor
      double precision regularizationMem,dataMatchMem,forceMatchMem
      double precision suml2grad1,suml2grad2
      double precision bc_step ! for elem611 loading -TS
      double precision delta_ws ! for weak springs
      double precision :: gmemn = 1.0d8! store the previously computed result of gradfun
      
c     FIRST BATCH
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ldofn
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elemvec_ndofn! localized information for every element
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itvect
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lnods,lrefn,idnum
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipoin_bc
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: icode_bc
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipoin_ws
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: icode_ws
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tvect
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: coord
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: coord_copy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: props
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: uelem_dual! localized information for every element
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: uelem_diff! localized information for every element
c      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: uelem_meas! localized information for every element
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: primal,dual
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: value_bc,meas
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: value_ws
c      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: value_bc2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: fac_meas
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: fac_meas_nd
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: total_primal
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: total_dual
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bodyforce! store the body forces

c     added new
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: elemDataMatch
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: elemRegul
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: l2grad1,l2grad2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: upperb,lowerb
c     end new

c     SECOND BATCH
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jdiag,jdiag_temp,nod_el_coun
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: aload
c     THIRD BATCH
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nod_el_conn
c     FOURTH BATCH
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jcsr, jcsr_temp
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: atang, atang_temp
c     FIFTH BATCH
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipoin_nod_grad
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ielem_nod_grad
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: adata_nodal
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: adata_nodal2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: delta_adata_nodal
c      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: grad_nodal! not used anymore
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: egrad
c     SIXTH BATCH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: adata_elemental
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xmem, gmem, fmem! store the computed results of gradfun
c     
c     DESCRITION OF VARIABLES (work in progress)
c     note: when a variable is after the flag openmp it means that it has been transformed
c           into a local variable for compatibility with openmp
c     
c     mprops: maximum number of material properties for any material set
c     mnode: maximum number of nodes per element for any element
c     mdofn: maximum number of dof for any node
c     mevab: maximum number of dof for any element
c     npoin: number of points (nodes)
c     nelem: number of elements
c openmp    nnode: number of nodes for the current element
c     nmats: number of material property sets
c openmp    nevab: number of dofs for the current element
c     ndime: number of dimensions of the problem
c openmp    nprop: number of material properties for the current element
c     neqns: number of eqns (total number of unknowns)
c     nnz: number of nonzeros in the stiffness matrix
c     nnz_conn: number of nonzeros in the dof to element connectivity
c     nset_nodal: number of nodal data sets
c     nset_elemental: number of elemental data sets
c     datpr: logical var to control printing of data

c     
c     DESCRITION OF ARRAYS (work in progress)
c     
c     ldofn(npoin): stores # of dofs for each node
c openmp    ldnum(mevab): global dof number for each elemental dof
c     elemvec_dofn(mnode): number of dof for each elemental node number
c openmp    lprop(nmat): number of material props for each material set
c     itvect(mdofn*npoin): temp. array 
c     lnods(nelem,mnode): global node number for each element node 
c     lrefn(nelem,3): element parameters for each element
c     idnum(npoin,mdofn): global dof number for each point
c openmp    eforc(mevab): element rhs vector
c     tvect(mdofn*npoin): temp. array
c openmp    pelem(mprops): element array to store material properties
c     coord(ndime,npoin): coordinates for each point
c openmp    xelem(ndime,mnode): coordinates for each element node
c     props(nmat,mprops): material properties for each material set
c openmp    estif(mevab,mevab): element stiffness
c openmp    uelem(mdofn,mnode): element solution array
c     uelem_meas
c     uelem_dual
c     uelem_diff
c     tdisp(npoin,mdofn): global solution array
c     jdiag(neqns+1): icsr array
c     nod_el_coun(neqns+1): pointers to nod_el_conn array
c     aload(neqns): the right-hand side
c     nod_el_conn(nnzconn): array to store dof to elem connectivity
c     jcsr(nnz): array to stode dof to dof connectivity
c     atang(nnz): global stiffness matrix
c     adata_nodal(nset_nodal,npoin): stores user-input nodal data 
c     count1: integer that determines, if the forward problem is solved 
c             via the step by step loading method or the continuation method,
c             count1 = 1 step by step loading is used,
c             count1 > 1 continuation method is used
c     fac_meas: multiplicative factor of associated displacement in the objective function
c     primal: stores the forward solution (variables)
c     dual: stores the dual solution
c     dataMatch: the term in the objective function that measures the difference in the computed and measured data
c     regularization: the term in the objective function that accounts for the regularization
c     elemDataMatch: the elemental contribution to dataMatch for every element
c     elemRegul: the elemental contribution to regularization for every element
c     l2grad2: L2 norm of the gradient of the second material property
c     ielem_nod_grad: JFD
c     ipoin_nod_grad: JFD

      END MODULE MAINMEM
      
      MODULE IOUNIT
      SAVE
      integer :: iread = 5
      integer :: iwrit = 6
      END MODULE IOUNIT
       
       





