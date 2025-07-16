#ifndef _PD6D_H_
#define _PD6D_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define MAXLENLINE 255

#ifdef __cplusplus
extern "C" {
#endif

// Function declarations
int PD6DCountAtoms(char *pdbfn, bool txt);
int PD6DReadxyz(FILE *file, int n, double tr[][3], double xr[][3], bool txt);
void pd6dmmTR(double tr[3], double xr[3], int nAtom, double xyzold[][3], double xyznew[][3]);

#ifdef __cplusplus
}
#endif

#endif /* _PD6D_H_ */