#include "pluto.h"

/**
 * @brief Dumps the current vectors
 */
void ComputeUserVar (const Data *d, Grid *grid)
{
  int i, j, k;
  double ***Jx1 = GetUserVar("Jx1");
  double ***Jx2 = GetUserVar("Jx2");
  double ***Jx3 = GetUserVar("Jx3");
  
  DOM_LOOP(k,j,i){
  	Jx1[k][j][i] = d->J[IDIR][k][j][i];
  	Jx2[k][j][i] = d->J[JDIR][k][j][i];
  	Jx3[k][j][i] = d->J[KDIR][k][j][i];
  }
}

//-----------------------------------------------------------------------------

void ChangeDumpVar ()
{ 

}
