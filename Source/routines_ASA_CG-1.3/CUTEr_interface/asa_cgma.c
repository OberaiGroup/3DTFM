/* ====================================================
 * CUTEr interface for cg_descent     Sept. 8, 2007
 *
 * W. Hager
 *
 * (Based on gencma.c of D. Orban, Feb 3, 2003)
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ASA_CGMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "../../include/cuter.h"
#include "../../include/asa_user.h"

#ifdef Isg95
#define MAINENTRY MAIN_
#else
#define MAINENTRY main
#endif

/* prototypes */
double asa_value
(
    asa_objective *asa
) ;

void asa_grad
(
    asa_objective *asa
) ;

double asa_valgrad
(
    asa_objective *asa
) ;

/* global variables */
    integer CUTEr_nvar;        /* number of variables */
    integer CUTEr_ncon;        /* number of constraints */

/* main program */
    int MAINENTRY( void ) {

        char *fname = "OUTSDIF.d"; /* CUTEr data file */
        integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
        integer iout = 6;          /* FORTRAN unit number for error output */
        integer ierr;              /* Exit flag from OPEN and CLOSE */

        VarTypes vtypes;

        integer    *indvar = NULL, *indfun = NULL, ncon_dummy ;
        doublereal *x, *lo, *hi, *dummy1, *dummy2 ;
        doublereal *v = NULL, *cl = NULL, *cu = NULL, *c = NULL, *cjac = NULL;
        logical    *equatn = NULL, *linear = NULL;
        char       *pname, *vnames, *gnames;
        logical     efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;
        logical     constrained = FALSE_;

        real        calls[7], cpu[2];
        integer     nlin = 0, nbnds = 0, neq = 0;
        doublereal  dummy;
        integer     ExitCode;
        int         i;

        asa_parm AParm ;
        asacg_parm CParm ;
        int status ;


        /* Open problem description file OUTSDIF.d */
        ierr = 0;
        FORTRAN_OPEN( &funit, fname, &ierr );
        if( ierr ) {
            printf("Error opening file OUTSDIF.d.\nAborting.\n");
            exit(1);
        }

        /* Determine problem size */
        CDIMEN( &funit, &CUTEr_nvar, &CUTEr_ncon );

        /* Determine whether to call constrained or unconstrained tools */
        if( CUTEr_ncon ) constrained = TRUE_;

        /* Seems to be needed for some Solaris C compilers */
        ncon_dummy = CUTEr_ncon + 1;

        /* Reserve memory for variables, bounds, and multipliers */
        /* and call appropriate initialization routine for CUTEr */
        MALLOC( x,      CUTEr_nvar, doublereal );
        MALLOC( lo,     CUTEr_nvar, doublereal );
        MALLOC( hi,     CUTEr_nvar, doublereal );

        MALLOC( equatn, 1, logical    );
        MALLOC( linear, 1, logical    );
        MALLOC( cl, 1, doublereal );
        MALLOC( cu, 1, doublereal );
        USETUP( &funit, &iout, &CUTEr_nvar, x, lo, hi, &CUTEr_nvar );

        /* Free unneeded arrays */
        FREE( equatn );
        FREE( linear );

        /* Get problem name */
        MALLOC( pname, FSTRING_LEN+1, char );
        MALLOC( vnames, CUTEr_nvar*FSTRING_LEN, char );
        UNAMES( &CUTEr_nvar, pname, vnames );
        FREE( vnames );

        /* Make sure to null-terminate problem name */
        pname[FSTRING_LEN] = '\0';
        i = FSTRING_LEN - 1;
        while( i-- > 0 && pname[i] == ' ') {
            pname[i] = '\0';
        }

        /* Set parameters */
        asa_default (&AParm) ;
        asa_cg_default (&CParm) ;
/*      AParm.PrintLevel = 4 ;*/
/*      CParm.PrintLevel = 4 ;*/

        /* Call the optimizer (if only the default parameter values are
           used, then "&cg_parm" could be replaced by NULL) */

        status = asa_cg (x, lo, hi, CUTEr_nvar, NULL, &CParm, &AParm, 1.e-6,
                         asa_value, asa_grad, asa_valgrad, NULL) ;

        ExitCode = 0;

        /* Get CUTEr statistics */
        CREPRT( calls, cpu );

        printf ("status: %5i CPU: %8.2f\n", status, cpu [1]) ;
/*
        printf("\n\n ************************ CUTEr statistics ************************\n\n");
        printf(" Code used               : cg_descent\n");
        printf(" Problem                 : %-s\n", pname);
        printf(" # variables             = %-10d\n", CUTEr_nvar);
        printf(" # bound constraints     = %-10d\n", vtypes.nbnds);
        printf(" # objective functions   = %-15.7g\n", calls[0]);
        printf(" # objective gradients   = %-15.7g\n", calls[1]);
        printf(" # objective Hessians    = %-15.7g\n", calls[2]);
        printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
        printf(" Exit code               = %-10d\n", ExitCode);
        printf(" Final f                 = %-15.7g\n",dummy);
        printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
        printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
        printf(" ******************************************************************\n\n");
*/

        ierr = 0;
        FORTRAN_CLOSE( &funit, &ierr );
        if( ierr ) {
            printf( "Error closing %s on unit %d.\n", fname, funit );
            printf( "Trying not to abort.\n" );
        }

        /* Free workspace */
        FREE( pname );
        FREE( x ); FREE( lo ); FREE( hi );
        FREE( v ); FREE( cl ); FREE( cu );

        return 0;
    }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
double asa_value
(
    asa_objective *asa
)
{
    long int N ;
    double f ;
    N = (long int) asa->n ;
    UFN ( &N, asa->x, &f ) ;
    return (f) ;
}

void asa_grad
(
    asa_objective *asa
)
{
    long int N ;
    N = (long int) asa->n ;
    UGR ( &N, asa->x, asa->g ) ;
}

double asa_valgrad
(
    asa_objective *asa
)
{
    long int grad ;
    long int N ;
    double f ;
    N = (long int) asa->n ;
    grad = (long int) 1 ;
    UOFG ( &N, asa->x, &f, asa->g, &grad ) ;
    return (f) ;
}
