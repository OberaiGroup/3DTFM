#include <math.h>
#include "asa_user.h" /* needed by the program which calls asa_cg*/

double myvalue (asa_objective *asa)/* evaluate the objective function */
{
    double f, t, *x ;
    INT i, n ;
    x = asa->x ;
    n = asa->n ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        f += exp (x [i]) - t*x [i] ;
    }
    return (f) ;
}
void mygrad (asa_objective *asa)/* evaluate the gradient of the objective function */
{
    double t, *g, *x ;
    INT i, n ;
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        g [i] = exp (x [i]) -  t ;
    }
    return ;
}

double myvalgrad (asa_objective *asa)/* value and gradient of the objective function */
{
    double f, xi, t, *g, *x ;
    INT i, n ;
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    f = 0 ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        xi = x [i] ;
        f += exp (xi) - t*xi ;
        g [i] = exp (xi) -  t ;
    }
    return (f) ;
}
