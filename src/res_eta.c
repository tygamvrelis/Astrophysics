/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Define the components of the diagonal resistive tensor. 

  Use this function to supply the resistivity in the three directions
  \f$ \eta_{x1}\f$, \f$ \eta_{x2}\f$ and \f$ \eta_{x3}\f$.
  
  \authors T. Matsakos \n
           A. Mignone (mignone@ph.unito.it)\n
  \date    March 22, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Resistive_eta(double *v, double x1, double x2, double x3,
                   double *J, double *eta)
/*!
 * Compute the resistive tensor components as function of the primitive
 * variables, coordinates and currents.
 *
 * \param [in]  v    array of primitive variables
 * \param [in]  x1   coordinate in the X1 direction
 * \param [in]  x2   coordinate in the X2 direction
 * \param [in]  x3   coordinate in the X3 direction
 * \param [in]  J    current components, J[IDIR], J[JDIR], J[KDIR]
 * \param [out] eta  an array containing the three components of
 *                   \f$ \tens{\eta}\f$.
 *
 *********************************************************************** */
{
#if defined(RES_INC_WITH_RADIUS)
  // Need to non-dimensionalize eta. Requires the profile to be inputted in
  // SI units
  const double LV = UNIT_LENGTH * UNIT_VELOCITY;
  const double ten_to_11_over_4_pi = pow(10, 11) / (4.0 * CONST_PI);
  const double eta_coeff = ten_to_11_over_4_pi / LV;

  double r = x1;
  const double R_MAX = g_domEnd[IDIR];
  double sigma = 1.e6 * exp(-r / (0.04 * R_MAX)) +
                 1.e-2 * exp((r - R_MAX) / (0.01 * R_MAX)); // SI
  double resistivity = eta_coeff / sigma;
  
  eta[IDIR] = resistivity;
  eta[JDIR] = resistivity;
  eta[KDIR] = resistivity;
#elif defined(RES_UNIFORM)
  eta[IDIR] = g_inputParam[ETA];
  eta[JDIR] = g_inputParam[ETA];
  eta[KDIR] = g_inputParam[ETA];
#else
  eta[IDIR] = 0.0;
  eta[JDIR] = 0.0;
  eta[KDIR] = 0.0; 
#endif
}
