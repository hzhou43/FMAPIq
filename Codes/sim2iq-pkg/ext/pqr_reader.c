#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "pqr_reader.h"
#include "cromer_mann_4gauss.h"

// Convert element string to element type
int get_element_type(const char* element) {
    // Skip whitespace
    while (*element && isspace(*element)) element++;
    
    // Extract first letter (should be element symbol)
    char first = *element;
    if (!first) return ELEMENT_UNKNOWN;
    
    // Convert to uppercase
    first = toupper(first);
     
    // Determine element type
    switch (first) {
        case 'H':
            return ELEMENT_H;
        case 'C':
            return ELEMENT_C;
        case 'N':
            return ELEMENT_N;
        case 'O':
            return ELEMENT_O;
        case 'P':
            return ELEMENT_P;
        case 'S':
            return ELEMENT_S;
        default:
            return ELEMENT_UNKNOWN;
    }
}

// Read atomic structure from PQR file
Structure read_pqr_file(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    Structure structure = structure_create();
    char line[MAX_LINE_LENGTH];
    
    // Use a larger buffer for better performance with large files
    setvbuf(file, NULL, _IOFBF, 32768);
    
    while (fgets(line, sizeof(line), file)) {
        // Skip lines that don't start with ATOM or HETATM
        if (strncmp(line, "ATOM", 4) != 0 && strncmp(line, "HETATM", 6) != 0) {
            continue;
        }
        
        structure_ensure_capacity(&structure);
        Atom* atom = &structure.atoms[structure.count];
        
        // Parse PQR format (based on http://www.poissonboltzmann.org/docs/file-format-info/)
        // Columns:
        // 1-6: ATOM/HETATM
        // 7-11: Atom serial number
        // 13-16: Atom name
        // 17: Alternate location indicator
        // 18-20: Residue name
        // 22: Chain identifier
        // 23-26: Residue sequence number
        // 27: Code for insertion of residues
        // 30-38: X coordinate (Å)
        // 39-46: Y coordinate (Å)
        // 47-54: Z coordinate (Å)
        // 55-61: Charge
        // 62-70: Radius
        
        char element[MAX_ELEMENT_LENGTH] = "";
        double radius = 0.0;
        
        // Extract atom name which contains element
        if (strlen(line) >= 16) {
            memcpy(element, line + 12, 4);
            element[4] = '\0';
            
            // Trim element
            char* p = element;
            while (*p && isspace(*p)) p++;
            memmove(element, p, strlen(p) + 1);
            
            p = element + strlen(element) - 1;
            while (p >= element && isspace(*p)) *p-- = '\0';
        }
        
        // Copy element to atom
        strncpy(atom->element, element, MAX_ELEMENT_LENGTH - 1);
        atom->element[MAX_ELEMENT_LENGTH - 1] = '\0';
        
        // Extract coordinates from fixed positions
        char xbuf[9] = {0};
        char ybuf[9] = {0};
        char zbuf[9] = {0};
        char rbuf[9] = {0};
        
        // Only copy if line is long enough
        if (strlen(line) >= 54) {
            // X: positions 30-38 (9 characters)
            memcpy(xbuf, line + 30, 8);
            xbuf[8] = '\0';
            
            // Y: positions 39-46 (8 characters)
            memcpy(ybuf, line + 38, 8);
            ybuf[8] = '\0';
            
            // Z: positions 47-54 (8 characters)
            memcpy(zbuf, line + 46, 8);
            zbuf[8] = '\0';
            
            // Convert to doubles
            atom->x = strtof(xbuf, NULL);
            atom->y = strtof(ybuf, NULL);
            atom->z = strtof(zbuf, NULL);
        } else {
            fprintf(stderr, "Warning: Line too short for coordinates: %s", line);
            continue;
        }
        
        // Extract radius from fixed position: 62-70 (9 characters)
        if (strlen(line) >= 70) {
            memcpy(rbuf, line + 62, 8);
            rbuf[8] = '\0';
            radius = strtof(rbuf, NULL);
        }
        
        // Set volume based on radius and element type
        int element_type = get_element_type(atom->element);
        const CromerMann4Gauss* elem_params = get_cromer_mann_params(element_type);
        
        if (radius > 0.0 && elem_params) {
            // If radius is in the PQR file, use it with element scaling factor
            atom->vol = (4.0/3.0) * M_PI * pow(radius*elem_params->scl, 3);
        } else if (elem_params) {
            // Otherwise use vdW radius from element parameters
            atom->vol = (4.0/3.0) * M_PI * pow(elem_params->vdw_radius*elem_params->scl, 3);
        } else {
            // Fallback: use a reasonable default volume
            atom->vol = (4.0/3.0) * M_PI * 1.7 * 1.7 * 1.7; // Default C vdw radius
        }
        
        structure.count++;
    }
    
    fclose(file);
    printf("Read %zu atoms from %s\n", structure.count, filename);
    return structure;
}