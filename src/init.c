#include <stdlib.h>
#include <R_ext/Rdynload.h>

/*.C function calls */
extern void loadX_R(void *, void *, void *);
extern void unloadX_R();
extern void closest_R(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
	{"loadX_R", 	(DL_FUNC) &loadX_R, 	3},
	{"unloadX_R",	(DL_FUNC) &unloadX_R, 	0},
	{"closest_R", 	(DL_FUNC) &closest_R,	6},
	{NULL, NULL, 0}
};

void R_init_ligp(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}
