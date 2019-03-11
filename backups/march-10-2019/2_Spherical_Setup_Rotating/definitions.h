#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     1

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

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define UNIT_LENGTH         (CONST_Rplanet/g_domEnd[IDIR])
#define LIMITER             VANLEER_LIM
