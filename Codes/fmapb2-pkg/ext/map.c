#include "rdf.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

// Global lattice data - loaded once and reused
static int g_lattice_loaded = 0;
static int g_Ind[1000][10];
static double g_Wght[1000][10];

/* One dimension mapping into the PBC box */
int PBC(int x, int mx) {
    int px;
    if (x >= 0 && x < mx) {
        return x;
    } else {
        px = x % mx;
        return px >= 0 ? px : px + mx;
    }
}

int signI(double x) {
    return (x < 0.0) ? -1 : 1;
}

// Function to find lattice.txt in multiple possible locations
char* find_lattice_file() {
    static char lattice_path[512];
    const char* possible_paths[] = {
        "lattice.txt",                                    // Current directory
        "data/lattice.txt",                              // Data subdirectory
        "/usr/local/share/fmapb2/lattice.txt", // System installation
        NULL  // End marker
    };
    
    // Try each possible path
    for (int i = 0; possible_paths[i] != NULL; i++) {
        if (access(possible_paths[i], R_OK) == 0) {
            strcpy(lattice_path, possible_paths[i]);
            printf("Found lattice file at: %s\n", lattice_path);
            return lattice_path;
        }
    }
     
    // Last resort: try Python package location (set by Python wrapper)
    const char* pkg_path = getenv("FMAPB2_DATA_PATH");
    if (pkg_path != NULL) {
        snprintf(lattice_path, sizeof(lattice_path), "%s/lattice.txt", pkg_path);
        if (access(lattice_path, R_OK) == 0) {
            //printf("Found lattice file in package: %s\n", lattice_path);
            return lattice_path;
        }
    }
    
    return NULL;  // Not found
}

// Initialize lattice data - call this once at startup
int init_lattice_data() {
    if (g_lattice_loaded) {
        return 1;  // Already loaded
    }
    
    char* lattice_file = find_lattice_file();
    if (lattice_file == NULL) {
        fprintf(stderr, "ERROR: Cannot find lattice.txt file!\n");
        fprintf(stderr, "Searched in current directory, data/, ../data/, and system paths.\n");
        fprintf(stderr, "You can set FMAPB2_DATA_PATH environment variable to specify location.\n");
        return 0;
    }
    
    FILE *file = fopen(lattice_file, "r");
    if (file == NULL) {
        fprintf(stderr, "ERROR: Cannot open lattice file: %s\n", lattice_file);
        return 0;
    }
    
    int pos = 0;
    while (20 == fscanf(file, "%d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf",
                       &g_Ind[pos][0], &g_Wght[pos][0], &g_Ind[pos][1], &g_Wght[pos][1], 
                       &g_Ind[pos][2], &g_Wght[pos][2], &g_Ind[pos][3], &g_Wght[pos][3], 
                       &g_Ind[pos][4], &g_Wght[pos][4], &g_Ind[pos][5], &g_Wght[pos][5],
                       &g_Ind[pos][6], &g_Wght[pos][6], &g_Ind[pos][7], &g_Wght[pos][7], 
                       &g_Ind[pos][8], &g_Wght[pos][8], &g_Ind[pos][9], &g_Wght[pos][9])) {
        pos++;
        if (pos >= 1000) break;  // Safety check
    }
    fclose(file);
    
    if (pos != 1000) {
        fprintf(stderr, "WARNING: Expected 1000 lattice entries, got %d\n", pos);
    }
    
    g_lattice_loaded = 1;
    //printf("Loaded %d lattice entries from %s\n", pos, lattice_file);
    return 1;
}

void Map2GrdR(int l, fftw_real* grid, double xyz[3], const double coeff, const double sys[4], double (*rdf)(double,double)) {
    fftw_real (*grd)[l][2*(l/2+1)] = (void*)grid;
    int xyzmin[3], xyzmax[3];
    int i, j, k, px, py, pz;
    double x, y, z, x2, y2, z2, r2;
    const double rlow = sys[0];
    const double rup = sys[1];
    const double dx = sys[2];
    const double rup2 = rup * rup;
    const double rlow2 = rlow * rlow;
    const double scl = sys[3];
    
    // Calculate the bounding box
    for (i = 0; i < 3; i++) {
        xyzmin[i] = (int)floor(xyz[i] - rup/dx);
        xyzmax[i] = (int)ceil(xyz[i] + rup/dx);
    }
    
    // Spherical iteration with early exits
    for (i = xyzmin[0]; i <= xyzmax[0]; i++) {
        x = (i - xyz[0]) * dx;  // Distance from atom to grid point
        x2 = x * x;
        if (x2 > rup2) continue;  // Early exit for x
        px = PBC(i, l);
        
        for (j = xyzmin[1]; j <= xyzmax[1]; j++) {
            y = (j - xyz[1]) * dx;  // Distance from atom to grid point
            y2 = y * y;
            double xy2 = x2 + y2;
            if (xy2 > rup2) continue;  // Early exit for x+y
            py = PBC(j, l);
            
            for (k = xyzmin[2]; k <= xyzmax[2]; k++) {
                z = (k - xyz[2]) * dx;  // Distance from atom to grid point
                z2 = z * z;
                r2 = xy2 + z2;
                if (r2 > rup2) continue;  // Final distance check
                
                pz = PBC(k, l);
                double var;
                if (r2 > rlow2) {
                    var = coeff * rdf(r2, scl);
                } else {
                    var = coeff * rdf(rlow2, scl);
                }
                #pragma omp atomic
                grd[px][py][pz] += var;
            }
        }
    }
}

// Last two parameters are useless; make Map2GrdL same as Map2GrdR
void Map2GrdL(int l, fftw_real* grid, double xyz[3], const double coeff, const double sys[4], double (*rdf)(double,double)) {
    fftw_real (*grd)[l][2*(l/2+1)] = (void*)grid;
    int px, py, pz;
    int pxn, pyn, pzn;
    double a, b, c;
    double av, bv, cv;
    double xyzi[3], xyzf[3];
    int i;
    
    for (i = 0; i < 3; i++) {
        xyzf[i] = modf(xyz[i], &xyzi[i]);
    }
    
    px = PBC((int)xyzi[0], l);
    py = PBC((int)xyzi[1], l);
    pz = PBC((int)xyzi[2], l);
    a = (double)xyzf[0];
    b = (double)xyzf[1];
    c = (double)xyzf[2];
    
    pxn = PBC(signI(a) + px, l);
    pyn = PBC(signI(b) + py, l);
    pzn = PBC(signI(c) + pz, l);
    
    av = fabs(a);
    bv = fabs(b);
    cv = fabs(c);
    
    #pragma omp atomic
    grd[px][py][pz] += coeff * (1.0 - av) * (1.0 - bv) * (1.0 - cv);
    #pragma omp atomic
    grd[px][py][pzn] += coeff * (1.0 - av) * (1.0 - bv) * cv;
    #pragma omp atomic
    grd[px][pyn][pz] += coeff * (1.0 - av) * bv * (1.0 - cv);
    #pragma omp atomic
    grd[px][pyn][pzn] += coeff * (1.0 - av) * bv * cv;
    #pragma omp atomic
    grd[pxn][py][pz] += coeff * av * (1.0 - bv) * (1.0 - cv);
    #pragma omp atomic
    grd[pxn][py][pzn] += coeff * av * (1.0 - bv) * cv;
    #pragma omp atomic
    grd[pxn][pyn][pz] += coeff * av * bv * (1.0 - cv);
    #pragma omp atomic
    grd[pxn][pyn][pzn] += coeff * av * bv * cv;
}

// Simplified setIndWght - now just copies from global data
int setIndWght(char* fname, int n, int Ind[][10], double Wght[][10]) {
    // Ensure lattice data is loaded
    if (!init_lattice_data()) {
        exit(EXIT_FAILURE);
    }
    
    // Copy from global arrays
    for (int i = 0; i < n && i < 1000; i++) {
        for (int j = 0; j < 10; j++) {
            Ind[i][j] = g_Ind[i][j];
            Wght[i][j] = g_Wght[i][j];
        }
    }
    
    return (n < 1000) ? n : 1000;
}

void Ind2ijk(int n, int ijk[3]) {
    ijk[0] = n / 4 / 4;
    ijk[1] = (n - ijk[0] * 4 * 4) / 4;
    ijk[2] = n - ijk[0] * 4 * 4 - ijk[1] * 4;
}

// assume 0<=xyz[?]<1.0
int xyzf2Pos(double xyzf[3]) {
    int i, j, k;
    i = (int)(xyzf[0] * 10.0);
    j = (int)(xyzf[1] * 10.0);
    k = (int)(xyzf[2] * 10.0);
    return i * 100 + j * 10 + k;
}

void Map2GrdL2nd(int l, fftw_real* grid, double xyz[3], const double coeff, const double sys[4], double (*rdf)(double,double)) {
    fftw_real (*grd)[l][2*(l/2+1)] = (void*)grid;
    int px, py, pz;
    double xyzi[3], xyzf[3];
    int low[3];
    int i;
    double w;
    static int nlattice = 0;
    static int Ind[1000][10];
    static double Wght[1000][10];
    
    if (nlattice == 0) {
        nlattice = setIndWght("lattice.txt", 1000, Ind, Wght);
    }
    
    for (i = 0; i < 3; i++) {
        xyzi[i] = floor(xyz[i]);
        xyzf[i] = xyz[i] - xyzi[i];
        low[i] = (int)xyzi[i];
    }
    
    int pos = xyzf2Pos(xyzf);
    int ijk[3];
    for (i = 0; i < 10; i++) {
        Ind2ijk(Ind[pos][i], ijk);
        w = Wght[pos][i];
        px = PBC(low[0] - 1 + ijk[0], l);
        py = PBC(low[1] - 1 + ijk[1], l);
        pz = PBC(low[2] - 1 + ijk[2], l);
        #pragma omp atomic
        grd[px][py][pz] += coeff * w;
    }
}

// Alternative optimized version using global data directly
void Map2GrdL2nd_optimized(int l, fftw_real* grid, double xyz[3], const double coeff, const double sys[4], double (*rdf)(double,double)) {
    fftw_real (*grd)[l][2*(l/2+1)] = (void*)grid;
    int px, py, pz;
    double xyzi[3], xyzf[3];
    int low[3];
    int i;
    double w;
    
    // Ensure lattice data is loaded
    if (!g_lattice_loaded) {
        if (!init_lattice_data()) {
            fprintf(stderr, "FATAL: Cannot load lattice data\n");
            exit(EXIT_FAILURE);
        }
    }
    
    for (i = 0; i < 3; i++) {
        xyzi[i] = floor(xyz[i]);
        xyzf[i] = xyz[i] - xyzi[i];
        low[i] = (int)xyzi[i];
    }
    
    int pos = xyzf2Pos(xyzf);
    int ijk[3];
    
    // Use global data directly - no need for local copies
    for (i = 0; i < 10; i++) {
        Ind2ijk(g_Ind[pos][i], ijk);
        w = g_Wght[pos][i];
        px = PBC(low[0] - 1 + ijk[0], l);
        py = PBC(low[1] - 1 + ijk[1], l);
        pz = PBC(low[2] - 1 + ijk[2], l);
        #pragma omp atomic
        grd[px][py][pz] += coeff * w;
    }
}

// Volume grid functions
void volGrd(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4],
        void (*map)(int,fftw_real*,double[3],const double, const double[4],double(*)(double,double)),
        double (*rdf)(double,double)) {
    int i;

    #pragma omp parallel for 
    for (i = 0; i < nAtom; i++) { 
        double rc2 = atoms[i].r * atoms[i].r;
        double sys_local[4] = {sys[0], atoms[i].r, sys[2], rc2};
        map(l, grd, atoms[i].xyz, 1.0, sys_local, rdf);
    }
}

void volGrdR(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4]) {
    volGrd(l, grd, nAtom, atoms, sys, Map2GrdR, rdf_vol);
}

void debyeGrd(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4],
        void (*map)(int,fftw_real*,double[3],const double, const double[4],double(*)(double,double)),
        double (*rdf)(double,double)) {
    int i;

    #pragma omp parallel for 
    for (i = 0; i < nAtom; i++) {
        map(l, grd, atoms[i].xyz, atoms[i].q, sys, rdf);
    }
}

void debyeGrdR(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4]) {
    debyeGrd(l, grd, nAtom, atoms, sys, Map2GrdR, rdf_debye);
}

void debyeGrdL2nd(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4]) {
    debyeGrd(l, grd, nAtom, atoms, sys, Map2GrdL2nd, NULL);
}

// Van der Waals grid functions
void vdwGrd(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4],
        void (*map)(int,fftw_real*,double[3],const double, const double[4],double(*)(double,double)),
        double (*rdf)(double,double), double n, double sign) {
    int i;

    #pragma omp parallel for 
    for (i = 0; i < nAtom; i++) {
        double sys_local[4] = {sys[0], sys[1], sys[2], n};
        double coeff = (n < 4.0) ? atoms[i].Bsq : atoms[i].Asq;
        map(l, grd, atoms[i].xyz, sign * coeff, sys_local, rdf);
    }
}

void r6GrdR(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4]) {
    vdwGrd(l, grd, nAtom, atoms, sys, Map2GrdR, rdf_vdw, 3.0, -1.0);
}

void r12GrdR(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4]) {
    vdwGrd(l, grd, nAtom, atoms, sys, Map2GrdR, rdf_vdw, 6.0, 1.0);
}

void r6GrdL(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4]) {
    vdwGrd(l, grd, nAtom, atoms, sys, Map2GrdL, NULL, 3.0, 1.0);
}

void r12GrdL(int l, fftw_real *grd, int nAtom, ATOM atoms[], double sys[4]) {
    vdwGrd(l, grd, nAtom, atoms, sys, Map2GrdL, NULL, 6.0, 1.0);
}

// Cleanup function (optional)
void cleanup_lattice_data() {
    g_lattice_loaded = 0;
    printf("Lattice data cleaned up\n");
}
