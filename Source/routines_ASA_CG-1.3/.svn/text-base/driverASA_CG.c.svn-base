
#include <math.h>
#include "asa_user.h" /* needed by the program which calls asa_cg*/

/* prototypes for the function and gradient evaluation routines */
double myvalue (asa_objective *asa) ;

void mygrad (asa_objective *asa) ;

double myvalgrad (asa_objective *asa) ;


void driver_asa_cg_(int *nnn,double *xx,double *low,double *hig,
                    double *mmiter)
{
    int nn=nnn[0];/* problem size*/
    double miter=mmiter[0];/* max number of iterations*/
    double gradTol=1.0e-15;/*tolerance on the gradient for stopping criterion*/

    /* if you want to change parameter value, you need the following: */
    asacg_parm cgParm ;
    asa_parm asaParm ;

    /* if you want to change parameter value, initialize strucs with default */
    asa_cg_default (&cgParm) ;
    asa_default (&asaParm) ;

    /* if you want to change parameters, change them here: */
    cgParm.totit_fac = miter ;
    cgParm.PrintParms = TRUE ;
    cgParm.PrintLevel = 0 ;
    asaParm.PrintParms = TRUE ;
    asaParm.PrintLevel = 0 ;

    /* run the code */
    asa_cg (xx, low, hig, nn, NULL, &cgParm, &asaParm,
                     gradTol, myvalue, mygrad, myvalgrad, NULL) ;


    /* if no change in parameters, you could replace Parm arguments by NULL*/
    /* reinitialize the starting guess */
    /*for (i = 0; i < nn; i++) xx [i] = 1 ;
    asa_cg (xx, low, hig, nn, NULL, NULL, NULL,
                     1.e-8, myvalue, mygrad, myvalgrad, NULL) ;*/

    /* with some loss of efficiency, you could omit the valgrad routine */
    /* reinitialize the starting guess */
    /*for (i = 0; i < nn; i++) xx [i] = 1 ;
    asa_cg (xx, low, hig, nn, NULL, NULL, NULL, 1.e-8, myvalue, mygrad, NULL, NULL);*/

    return;
    
}

