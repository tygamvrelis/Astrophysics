#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     1
#define  INTERNAL_BOUNDARY       YES

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          YES

/* -- user-defined parameters (labels) -- */

#define  ALPHA                   0

/* [Beg] user-defined constants (do not change this line) */

#define CONST_Rjupiter           6.9911e9
#define CONST_Rplanet           (1.6*CONST_Rjupiter)

/** Number of seconds for planet to complete a rotation about its sun */
#define CONST_period            (3.0*24.0*3600.0)

// Enable winds
//#define WIND_EXP  // exponential profile
//#define WIND_QUAD // quadratic profile

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

// With unit length defined this way, CONST_Rplanet will measure the radius of
// the planet all the way up to the end of the R-coordinate, as desired. At the
// same time, R=1 will not be the edge of the planet (R=g_domEnd[IDIR] will be)
// but will instead define the location where the gravity body force changes,
// as before
#define UNIT_LENGTH         (CONST_Rplanet/g_domEnd[IDIR])
#define LIMITER             VANLEER_LIM
