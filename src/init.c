/**
 * @file
 * @brief Based on the stratified atmosphere test problem
 */

#include <float.h>
#include <stdbool.h>
#include "pluto.h"

// Polynomial coefficient for g when r < 1. With this choice, g and g' will be
// continuous
static const real acf = -3.0;
static const real bcf =  2.0;
static const real ccf = 0.0;

double init_density(double x1, double x2, double x3)
{
  double rs = x1;
  double scrh;
  double alpha = g_inputParam[ALPHA];
  double rho;
  if (rs > 1.0)
  {
    rho = exp(alpha*(1.0/rs - 1.0));
  }
  else
  {
    scrh = 0.5*acf*(rs*rs - 1.0) + 1.0/3.0*bcf*(rs*rs*rs - 1.0)
              + 0.25*ccf*(rs*rs*rs*rs - 1.0);
    rho = exp(alpha*scrh);
  }

  #if ROTATING_FRAME == YES
    rho *= exp(alpha*0.5*pow(g_OmegaZ,2)*pow(rs,2)*pow(sin(x2),2));
  #endif

  return rho;
}

double init_pressure(double rho)
{
  return rho / g_inputParam[ALPHA];
}

static double min_P = DBL_MAX;
void Init (double *us, double x1, double x2, double x3)
{
  #if ROTATING_FRAME == YES
    g_OmegaZ = (1.0/CONST_period)*UNIT_LENGTH/UNIT_VELOCITY;
  #endif

  us[RHO] = init_density(x1, x2, x3);

  double P = init_pressure(us[RHO]);
  us[PRS] = P;
  min_P = (P < min_P) ? P : min_P;

  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;

  us[TRC] = 0.0;

  #if PHYSICS == MHD
    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.0;

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = 0.0;
  #endif

  #if defined(B_DIPOLE) && (BACKGROUND_FIELD == NO)
    DipoleField (x1, x2, x3, &v[BX1], &v[BX2], &v[AX3]);
  #endif
}

//-----------------------------------------------------------------------------

#if defined(B_UNIFORM) || defined (B_DIPOLE)
// Sets initial curl-free magnetic field component
// Runs after init.
void BackgroundField (double x1, double x2, double x3, double *B0)
{
  #if defined(B_UNIFORM)
  {
    // Magnetic pressure: Bz << sqrt(8*pi*P)
    double Bz = sqrt(8.0 * CONST_PI * min_P) / 1000.0;
    B0[IDIR] = cos(x2) * Bz;      // r
    B0[JDIR] = -1 * sin(x2) * Bz; // theta
    B0[KDIR] = 0;                 // phi
  }
  #elif defined(B_DIPOLE)
  {
    DipoleField(x1, x2, x3, &B0[IDIR], &B0[JDIR], NULL);
  }
  #else
    B0[IDIR] = 0.0;
    B0[JDIR] = 0.0;
    B0[KDIR] = 0.0;
  #endif
}
#endif

//-----------------------------------------------------------------------------

// Based on init.c in Test_Problems/MHD/Torus
// Initializes dipole field in spherical coordinates
void DipoleField(double x1, double x2, double x3, 
                 double *Bx1, double *Bx2, double *A)
{
  static bool first = true;
  static double M;
  if (first)
  {
    first = false;
    // Convert input to code units
    double B = g_inputParam[BSURFACE] / (
      UNIT_VELOCITY * sqrt(4 * CONST_PI * UNIT_DENSITY)
    );

    // |M| = (B * r^3) / (1 + 3*(cos(theta))^2)
    // At equator (theta = 90), this reduces to |M| = B * r^3
    M = B * pow(g_domEnd[IDIR], 3);
  }

  double r_cubed = pow(x1, 3);
  *Bx1 = 2.0 * M *cos(x2) / r_cubed; // b_r
  *Bx2 = M * sin(x2) / r_cubed;      // b_phi
  if (A != NULL)
  {
    *A = *Bx2 * x1; // M * sin(x2) / pow(r, 2);
  }
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
