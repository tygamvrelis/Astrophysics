/**
 * @file
 * @brief Based on the stratified atmosphere test problem
 */

#include "pluto.h"

// Polynomial coefficient for g when r < 1. With this choice, g and g' will be
// continuous
static const real acf = -3.0;
static const real bcf =  2.0;
static const real ccf = 0.0;

void Init (double *us, double x1, double x2, double x3)
{
  double rs = x1;
  double scrh;
  double alpha = g_inputParam[ALPHA];
  if (rs > 1.0)
  {
    us[RHO] = exp(alpha*(1.0/rs - 1.0));
  }
  else
  {
    scrh = 0.5*acf*(rs*rs - 1.0) + 1.0/3.0*bcf*(rs*rs*rs - 1.0)
              + 0.25*ccf*(rs*rs*rs*rs - 1.0);
    us[RHO] = exp(alpha*scrh);
  }

  #if ROTATING_FRAME == YES
   g_OmegaZ = (1.0/CONST_period)*UNIT_LENGTH/UNIT_VELOCITY;
   us[RHO] *= exp(alpha*0.5*pow(g_OmegaZ,2)*pow(rs,2)*pow(sin(x2),2));
  #endif
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[PRS] = us[RHO]/alpha;
  us[TRC] = 0.0;
}

//-----------------------------------------------------------------------------

void Analysis (const Data *d, Grid *grid)
{

}

//-----------------------------------------------------------------------------

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int   i, j, k;
  double *x1, *x2, *x3;
  double rs;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  if (side == X1_END)
  {
    double alpha = g_inputParam[ALPHA];
    double theta;
    X1_END_LOOP(k,j,i)
    {
      rs = x1[i];
      theta = x2[j];
      d->Vc[RHO][k][j][i] = exp(alpha*(1.0/rs - 1.0));
      #if defined(ROTATING_FRAME)
       d->Vc[RHO][k][j][i] *= exp(alpha*0.5*pow(g_OmegaZ,2)*pow(rs,2)*pow(sin(theta),2));
      #endif
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]/alpha;
    }
  }
}

//-----------------------------------------------------------------------------

#if (BODY_FORCE & VECTOR)
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
    double gs;
    double rs = x1;
    double theta = x2;
    if (rs > 1.0)
    {
        gs = -1.0/rs/rs;
    }
    else
    {
        gs = rs*(acf + rs*(bcf + rs*ccf));
    }

    g[IDIR] = gs;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;

#if defined(WIND_EXP) || defined(WIND_QUAD)
    // We scale the sin argument so that it is 0 wherever the end of our
    // elevation grid is, and it is maximal at the equator
    // 
    // Thus, m and b are derived from:
    //     1. sin(g_domBeg[JDIR]*m + b) = 0
    //     2. sin((pi/2)*m + b) = 1
    // 
    // Note that g_domBeg[JDIR] is the coordinate value the elevation grid
    // begins at.
    const double R_MAX = g_domEnd[IDIR];
    if (rs >= 0.9*R_MAX)
    {
        const double m = 1.0/(1.0 - 2.0*g_domBeg[JDIR]/CONST_PI);
        const double b = -1.0*g_domBeg[JDIR]*m;
        const double T_RELAX = g_inputParam[TRELAX];
        const double V_MAX = g_inputParam[VMAX]; // Maximum speed in km/s
    #if defined(WIND_QUAD)
        double ampl = rs*rs/(R_MAX*R_MAX);
    #else
        const double decay_rate = g_inputParam[EXP_DECAY];
        double ampl = exp(decay_rate*(rs/R_MAX - 1.0));
    #endif
        double v_relax = V_MAX*ampl*sin(theta*m + b);
        g[KDIR] = (v_relax - v[VX3]) / T_RELAX;
    }
#endif
}
#endif

//-----------------------------------------------------------------------------

#if (BODY_FORCE & POTENTIAL)
double BodyForcePotential(double x1, double x2, double x3)
{
    // integration constant to make phi continuous
    static const double C = 0.5*acf + bcf/3.0 + ccf*0.25;

    double phi;
    double rs = x1;
    if (rs > 1.0)
    {
        phi = -1.0/rs;
    }
    else
    {
        phi = -rs*rs*(0.5*acf + rs*(bcf/3.0 + rs*ccf*0.25)) + C - 1.0;
    }

    return phi;
}
#endif
