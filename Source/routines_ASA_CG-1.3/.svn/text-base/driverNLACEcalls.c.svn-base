#include <math.h>
#include "asa_user.h" /* needed by the program which calls asa_cg*/


double myvalue (asa_objective *asa)/* evaluate the objective function */
{
    double f, *g, *x;
    int itN, n;
    f=0.0;
    /*fprintf(stdout,"going through myvalue\n");*/
    /*double t;
    INT i, n ;*/
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    /*fprintf(stdout,"n=%d\n",n);*/
    itN = asa->itNum ;
    /*fprintf(stdout,"itN=%d\n",itN);*/
    /*n = asa->n ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        f += exp (x [i]) - t*x [i] ;
    }*/
    /* call gradfun */
    fflush(stdout);
    gradfun_(x,&f,g,&itN);
    return (f) ;
}
void mygrad (asa_objective *asa)/* evaluate the gradient of the objective function */
{
    double f, *g, *x;
    int itN, n ;
    f=0.0;
    /*fprintf(stdout,"going through mygrad\n");*/
    /*double t;
    INT i, n ;*/
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    /*fprintf(stdout,"n=%d\n",n);*/
    itN = asa->itNum ;
    /*fprintf(stdout,"itN=%d\n",itN);*/
    /*n = asa->n ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        g [i] = exp (x [i]) -  t ;
    }*/
    /* call gradfun */
    fflush(stdout);
    gradfun_(x,&f,g,&itN);
    return ;
}

double myvalgrad (asa_objective *asa)/* value and gradient of the objective function */
{
    double f, *g, *x;
    int itN, n ;
    f=0.0;
    /*fprintf(stdout,"going through myvalgrad\n");*/
    /*double f, xi;
    INT i, n ;*/
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    /*fprintf(stdout,"n=%d\n",n);*/
    itN = asa->itNum ;
    /*fprintf(stdout,"x[0]=%e, x[1]=%e\n",x[0],x[1]);
    fprintf(stdout,"f=%e\n",f);
    fprintf(stdout,"g[0]=%e, g[1]=%e\n",g[0],g[1]);
    fprintf(stdout,"itN=%d\n",itN);*/
    /*n = asa->n ;
    f = 0 ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        xi = x [i] ;
        f += exp (xi) - t*xi ;
        g [i] = exp (xi) -  t ;
    }*/
    /* call gradfun */
    fflush(stdout);
    gradfun_(x,&f,g,&itN);
    return (f) ;
}
