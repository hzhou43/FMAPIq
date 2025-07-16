#include "pd6d.h"

// Utility function to check if a file exists
int FileExist(char *filename) {
    FILE *file;
    if ((file = fopen(filename, "r"))) {
        fclose(file);
        return 1;
    }
    return 0;
}

// PD6D Line Reading Function
int PD6DReadline(char* line, double tr[3], double xr[3], bool txt) {
    int id = 0;
    
    if (txt) {
        sscanf(line, "%16d%16lf%16lf%16lf%16lf%16lf%16lf", 
               &id, &tr[0], &tr[1], &tr[2], &xr[0], &xr[1], &xr[2]);
    } else {
        sscanf(line+6, "%5d", &id);
        sscanf(line+30, "%8lf%8lf%8lf%8lf%8lf%8lf", 
               &tr[0], &tr[1], &tr[2], &xr[0], &xr[1], &xr[2]);
    }
    return id;
}

// Count atoms in a MODEL
int PD6DCountAtoms(char *pdbfn, bool txt) {
    FILE *file;
    char line[MAXLENLINE];
    int num, id = 0;
    
    if (FileExist(pdbfn) == 0) {
        exit(EXIT_FAILURE);
    }
    
    file = fopen(pdbfn, "r");
    num = 0;
    int preid = 0;
    double tr[3], xr[3];
    
    while (fgets(line, MAXLENLINE, file) != NULL) {
        if (txt) {
            if (!strncmp(line, "MODEL", 5) || !strncmp(line, "ENDMDL", 6)) {
                continue;
            } else {
                id = PD6DReadline(line, tr, xr, txt);
                if (id > preid) {
                    preid = id;
                    num++;
                } else {
                    break;
                }
            }
        } else {
            if (!strncmp(line, "ATOM", 4) || !strncmp(line, "HETATM", 6)) { 
                id = PD6DReadline(line, tr, xr, txt);
                if (id > preid) {
                    preid = id;
                    num++;
                } else {
                    break;
                }
            }
        }
    }
    
    if (num != preid) {
        fprintf(stderr, "%s:%d:Number of atoms:%d != id %d \n", __FILE__, __LINE__, num, id);
    }
    
    fclose(file);
    return num;
}

// Read XYZ coordinates
int PD6DReadxyz(FILE *file, int n, double tr[][3], double xr[][3], bool txt) {
    char line[MAXLENLINE];
    int i = 0;
    int id = 0;
    
    while (fgets(line, MAXLENLINE, file) != NULL) {
        if (txt) {
            if (!strncmp(line, "MODEL", 5) || !strncmp(line, "ENDMDL", 6)) {
                continue;
            } else {
                id = PD6DReadline(line, tr[i], xr[i], txt);
                i++;
            }
        } else {
            if (!strncmp(line, "ATOM", 4) || !strncmp(line, "HETATM", 6)) {
                id = PD6DReadline(line, tr[i], xr[i], txt);
                i++;
            }
        }
        
        if (id == n) {
            break;
        }
    }
    
    i = (i == n) ? i : 0;
    return i;
}

void Euler2Rot(double psi, double theta, double phi, double rot[9]) {
    double r11, r21, r31, r12, r22, r32, r13, r23, r33;
    
    r11 = cos(psi) * cos(phi) - sin(psi) * cos(theta) * sin(phi);
    r21 = sin(psi) * cos(phi) + cos(psi) * cos(theta) * sin(phi);
    r31 = sin(theta) * sin(phi);

    r12 = -cos(psi) * sin(phi) - sin(psi) * cos(theta) * cos(phi);
    r22 = -sin(psi) * sin(phi) + cos(psi) * cos(theta) * cos(phi);
    r32 = sin(theta) * cos(phi);

    r13 = sin(psi) * sin(theta);
    r23 = -cos(psi) * sin(theta);
    r33 = cos(theta);

    rot[0] = r11; rot[1] = r12; rot[2] = r13;
    rot[3] = r21; rot[4] = r22; rot[5] = r23;
    rot[6] = r31; rot[7] = r32; rot[8] = r33;
}

// Apply transformation and rotation to coordinates using rotation matrix in format rot[9]
void xyz_t_r(double xyz[3], double tran[3], double rot[9], double trxyz[3]) {
    trxyz[0] = tran[0] + xyz[0] * rot[0] + xyz[1] * rot[1] + xyz[2] * rot[2];
    trxyz[1] = tran[1] + xyz[0] * rot[3] + xyz[1] * rot[4] + xyz[2] * rot[5];
    trxyz[2] = tran[2] + xyz[0] * rot[6] + xyz[1] * rot[7] + xyz[2] * rot[8];
}

// Apply transformation and rotation to multiple coordinates
void xyzs_t_r(double lxyz[][3], int nums, double tran[3], double rot[9], double ltrxyz[][3]) {
    int i;
    for (i = 0; i < nums; i++) {
        xyz_t_r(lxyz[i], tran, rot, ltrxyz[i]);
    }
}

// Main function for pd6dmm transformation and rotation
void pd6dmmTR(double tr[3], double xr[3], int nAtom, double xyzold[][3], double xyznew[][3]) {
    double rot[9];
    
    Euler2Rot(xr[0],xr[1],xr[2],rot);
    xyzs_t_r(xyzold, nAtom, tr, rot, xyznew);
}