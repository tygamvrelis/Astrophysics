#include "pluto.h"
#include <stdbool.h>

/**
 * @brief Gets the current, background, and resistivity fields
 */
void ComputeUserVar (const Data *d, Grid *grid)
{
  double ***Jx1 = GetUserVar("Jx1");
  double ***Jx2 = GetUserVar("Jx2");
  double ***Jx3 = GetUserVar("Jx3");
  
  // Get current
  int i, j, k;
  DOM_LOOP(k,j,i){
#if RESISTIVITY != NO
  	Jx1[k][j][i] = d->J[IDIR][k][j][i];
  	Jx2[k][j][i] = d->J[JDIR][k][j][i];
  	Jx3[k][j][i] = d->J[KDIR][k][j][i];
#else
  	Jx1[k][j][i] = 0;
  	Jx2[k][j][i] = 0;
  	Jx3[k][j][i] = 0;
#endif
  }

  // Only compute background and resistivity fields the first time
  static bool first = true;
  if (first)
  {
    first = false;
    double ***B0x1  = GetUserVar("B0x1");
    double ***B0x2  = GetUserVar("B0x2");
    double ***B0x3  = GetUserVar("B0x3");
    double ***etax1 = GetUserVar("etax1");
    double ***etax2 = GetUserVar("etax2");
    double ***etax3 = GetUserVar("etax3");
    double* x1 = grid[IDIR].xgc;
    double* x2 = grid[JDIR].xgc;
    double* x3 = grid[KDIR].xgc;
    double B0[3], eta[3];
    DOM_LOOP(k,j,i){
      #if RESISTIVITY != NO
        Resistive_eta(NULL, x1[i], x2[j], x3[k], NULL, eta);
        etax1[k][j][i] = eta[IDIR];
        etax2[k][j][i] = eta[JDIR];
        etax3[k][j][i] = eta[KDIR];
      #else
        etax1[k][j][i] = 0;
        etax2[k][j][i] = 0;
        etax3[k][j][i] = 0;
      #endif
        BackgroundField(x1[i], x2[j], x3[k], B0);
        B0x1[k][j][i] = B0[IDIR];
        B0x2[k][j][i] = B0[JDIR];
        B0x3[k][j][i] = B0[KDIR];
    }
  }
}

//-----------------------------------------------------------------------------

void ChangeDumpVar ()
{ 

}
