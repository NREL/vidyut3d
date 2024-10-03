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
#define NUM_PLASMAVARS (REF_ID-NUM_SPECIES+1)
#define NVAR (NUM_SPECIES+NUM_PLASMAVARS)

#define PERBC 0
#define DIRCBC 1
#define HNEUBC 2
#define IHNEUBC 3
#define ROBINBC 4

#endif
