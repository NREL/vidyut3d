#ifndef Vidyut_vardefines_H_
#define Vidyut_vardefines_H_
#include <Chemistry.H>

#define EEN_ID (NUM_SPECIES)
#define ETEMP_ID (NUM_SPECIES+1)
#define POT_ID (NUM_SPECIES+2)
#define EFX_ID (NUM_SPECIES+3)
#define EFY_ID (NUM_SPECIES+4)
#define EFZ_ID (NUM_SPECIES+5)
#define EJH_ID (NUM_SPECIES+6) //electron joule heating
#define EIH_ID (NUM_SPECIES+7) //electron inelastic heating
#define EEH_ID (NUM_SPECIES+8) //electron elastic heating
#define REF_ID (NUM_SPECIES+9)
#define SRFCH_ID (NUM_SPECIES+10)
#define PHOTO_ION_SRC_ID (NUM_SPECIES+11)
#define CMASK_ID (NUM_SPECIES+12)
#define NUM_PLASMAVARS (CMASK_ID-NUM_SPECIES+1)
#define NVAR (NUM_SPECIES+NUM_PLASMAVARS)

#define PERBC 0
#define DIRCBC 1
#define HNEUBC 2
#define IHNEUBC 3
#define ROBINBC 4
#define AXISBC 5

#define POTEQ_TYPE 0
#define PHOTEQ_TYPE 1
#define SPECEQ_TYPE 2

#endif
