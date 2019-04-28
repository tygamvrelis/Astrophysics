#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     5

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          YES

/* -- user-defined parameters (labels) -- */

#define  ALPHA                   0
#define  VMAX                    1
#define  EXP_DECAY               2
#define  TRELAX                  3
#define  BSURFACE                4

/* [Beg] user-defined constants (do not change this line) */

// ----------------------------------------------------------------------------
// Winds: Exponential or quadratic profile, or none
// ----------------------------------------------------------------------------
#define  WIND_EXP
// #define WIND_QUAD

// ----------------------------------------------------------------------------
// Magnetic field: uniform vertical field, dipole, or none
// ----------------------------------------------------------------------------
#define B_DIPOLE
// #define B_UNIFORM

// ----------------------------------------------------------------------------
// Planet parameters
// ----------------------------------------------------------------------------
#define  CONST_Rjupiter          6.9911e9
#define  CONST_Rplanet           (1.6*CONST_Rjupiter)
#define  CONST_period            (3.0*24.0*3600.0)

// ----------------------------------------------------------------------------
// Scaling
// ----------------------------------------------------------------------------
// With unit length defined this way, CONST_Rplanet will measure the radius of
// the planet all the way up to the end of the R-coordinate, as desired. At the
// same time, R=1 will not be the edge of the planet (R=g_domEnd[IDIR] will be)
// but will instead define the location where the gravity body force changes,
// as before
#define UNIT_LENGTH         (CONST_Rplanet/g_domEnd[IDIR])
#define UNIT_DENSITY        1e-9

#define ASSIGN_VECTOR_POTENTIAL YES

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 
#define  LIMITER                   VANLEER_LIM
#define  CT_EMF_AVERAGE            UCT_HLL
