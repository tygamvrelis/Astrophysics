/**
 * @file
 * @brief Based on the stratified atmosphere test problem
 */

#include "pluto.h"

static real acf, bcf, ccf; /*  polynomial coefficient for g when r < 1 */

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double rs, scrh;

/* ------------------------------------------
    with this choice g, g' and g''
    will be continuous at R = 1
   ----------------------------------------- */

  acf = -10.0;
  bcf =  15.0;
  ccf = -6.0;

/* -----------------------------------
     with this choice g and g' 
     will be continuous
   ---------------------------------- */

  acf = -3.0;
  bcf =  2.0;
  ccf = 0.0;

  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL 
   rs = sqrt(x1*x1 + x2*x2);
  #elif GEOMETRY == SPHERICAL
   rs = sqrt(x1*x1);
  #endif

  if (rs > 1.0){
    us[RHO] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
  }else{
    scrh   = 0.5*acf*(rs*rs - 1.0) + 1.0/3.0*bcf*(rs*rs*rs - 1.0)
              + 0.25*ccf*(rs*rs*rs*rs - 1.0);
    us[RHO] = exp(scrh*g_inputParam[ALPHA]);
  }

  us[VX1] = 0.0;
  us[VX2] = 0.0;
  #if ROTATING_FRAME == YES
   g_OmegaZ = (1.0 / CONST_period) * UNIT_LENGTH / UNIT_VELOCITY;
  #endif   
  us[VX3] = 0.0;
  us[PRS] = us[RHO]/g_inputParam[ALPHA];
  us[TRC] = 0.0;
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int   i, j, k;
  double *x1, *x2, *x3;
  double rs;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  if (side == 0) {
    double theta;
    // We scale the sin argument so that it is 0 wherever the end of our
    // elevation grid is, and it is maximal at the equator
    // 
    // Thus, m and b are derived from:
    //     1. sin(g_domBeg[KDIR]*m + b) = 0
    //     2. sin((pi/2)*m + b) = 1
    // 
    // Note that g_domBeg[KDIR] is the coordinate value the elevation grid
    // begins at.
    const double m = 1.0/(1.0 - 2.0*g_domBeg[KDIR]/CONST_PI);
    const double b = -1.0*g_domBeg[KDIR]*m;
    const double RMAX = g_domEnd[IDIR];
    TOT_LOOP(k, j, i) {
      rs = x1[i];
      if (rs >= 0.9*RMAX)
      {
        theta = x2[j];
        d->Vc[VX3][k][j][i] = exp(10.0*(rs/RMAX - 1.0))*sin(theta*m + b);
      }
    }
  } else if (side == X1_END) {
    X1_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == CYLINDRICAL
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == SPHERICAL
       rs = sqrt(x1[i]*x1[i]);
      #endif
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }

  } else if (side == X2_END) {
    #if GEOMETRY != SPHERICAL
    X2_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == CYLINDRICAL
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == SPHERICAL
       rs = sqrt(x1[i]*x1[i]);
      #endif
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }
    #endif

  } else if (side == X3_END) {   /* Only Cartesian */

    X3_END_LOOP(k,j,i){  
      rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double gs, rs;
  double acf, bcf, ccf;

  acf = -3.0;
  bcf =  2.0;
  ccf =  0.0;
  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
   rs = sqrt(x1*x1 + x2*x2);
  #elif GEOMETRY == SPHERICAL
   rs = sqrt(x1*x1);
  #endif

  if (rs > 1.0) gs = -1.0/rs/rs;
  else          gs = rs*(acf + rs*(bcf + rs*ccf));

  #if GEOMETRY == CARTESIAN
   g[IDIR] = gs*x1/rs;
   g[JDIR] = gs*x2/rs;
   g[KDIR] = gs*x3/rs;
  #elif GEOMETRY == CYLINDRICAL
   g[IDIR] = gs*x1/rs;
   g[JDIR] = gs*x2/rs;
   g[KDIR] = 0.0;
  #elif GEOMETRY == SPHERICAL
   g[IDIR] = gs;
   g[JDIR] = 0.0; // + fz; // No idea if this is where fz should go...
   g[KDIR] = 0.0;
  #endif

}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double rs, phi;
  double acf, bcf, ccf, C;

  acf = -3.0;
  bcf =  2.0;
  ccf =  0.0;

  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
   rs = sqrt(x1*x1 + x2*x2);
  #elif GEOMETRY == SPHERICAL
   rs = sqrt(x1*x1);
  #endif

  C = (0.5*acf + bcf/3.0 + ccf*0.25);  /* integration constant to make phi continuous */
  if (rs > 1.0) phi = -1.0/rs;
  else          phi = -rs*rs*(0.5*acf + rs*(bcf/3.0 + rs*ccf*0.25)) + C - 1.0;

  const double cf = 0.5*pow(x1*sin(x2), 2)*g_OmegaZ*g_OmegaZ;
  phi += cf;
  return phi;
}
#endif
