#ifndef PQR_READER_H
#define PQR_READER_H

#include "saxs_calculator.h"

// Convert element string to element type
int get_element_type(const char* element);

// Read atomic structure from PQR file
Structure read_pqr_file(const char* filename);

#endif /* PQR_READER_H */