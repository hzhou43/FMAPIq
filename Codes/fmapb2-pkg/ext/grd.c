#include "common.h" 

// Static cache structure - automatically initialized to zero/NULL
static struct {
    fftw_real *crowd_vol;
    fftw_real *crowd_ele;
    fftw_real *crowd_vdw;
    fftw_real *crowd_v12;
    bool vol_computed;
    bool ele_computed;
    bool vdw_computed;
    bool v12_computed;
    int cached_l;
} cache = {0};

// Unified grid computation function
static void compute_grid(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, 
                        fftw_real **crowd_cache, bool *computed_flag,
                        void (*crowd_func)(int, fftw_real*, int, ATOM[], double[]),
                        void (*pro_func)(int, fftw_real*, int, ATOM[], double[]),
                        double sys[4], fftw_real Pro[l][l][2*(l/2+1)]) {
    
    // Check if grid size changed and cleanup if needed
    if (cache.cached_l != 0 && cache.cached_l != l) {
        // Grid size changed - cleanup old caches to avoid memory leak
        softgrd_cleanup();
    }
    
    // Allocate or reallocate if grid size changed
    if (!*computed_flag || cache.cached_l != l) {
        if (*crowd_cache) free(*crowd_cache);
        *crowd_cache = malloc(l * l * 2*(l/2+1) * sizeof(fftw_real));
        if (!*crowd_cache) {
            fprintf(stderr, "Error: Failed to allocate memory for crowd cache\n");
            return;
        }
        
        if (cache.cached_l != l) {
            cache.cached_l = l;
            *computed_flag = false;
        }
    }
    
    fftw_real (*Crowd)[l][2*(l/2+1)] = (fftw_real (*)[l][2*(l/2+1)])*crowd_cache;
    SetZero(l,Pro);
    
    // Compute Crowd data only once
    if (!*computed_flag) {
        SetZero(l,Crowd);
        crowd_func(l, &Crowd[0][0][0], nCrd, Crds, sys);
        fft3d_r2c(l,l,l,Crowd);
        *computed_flag = true;
    }
    
    pro_func(l, &Pro[0][0][0], nPro, Pros, sys);
    fft3d_r2c(l,l,l,Pro);
    fft3d_add_inplace(l,l,l,Crowd,Pro);
    fft3d_c2r(l,l,l,Pro);
}

void grdvol(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, SOLV* sol, fftw_real vol_result[l][l][2*(l/2+1)]){
    double sys[4] = {0.0, 0.0, dx, 0.0};
    
    compute_grid(nCrd, Crds, nPro, Pros, l, 
                &cache.crowd_vol, &cache.vol_computed,
                volGrdR, volGrdR, sys, vol_result);
}

void grdele(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, SOLV* sol, fftw_real ele_result[l][l][2*(l/2+1)]){
    double sys[4] = {1.0, sol->Cut, dx, sol->kappa};
    
    compute_grid(nCrd, Crds, nPro, Pros, l, 
                &cache.crowd_ele, &cache.ele_computed,
                debyeGrdR, debyeGrdL2nd, sys, ele_result);
}

void grdvdw(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, SOLV* sol, fftw_real vdw_result[l][l][2*(l/2+1)]){
    double sys[4] = {1.0, sol->Cut, dx, 1.0};
    
    compute_grid(nCrd, Crds, nPro, Pros, l, 
                &cache.crowd_vdw, &cache.vdw_computed,
                r6GrdR, r6GrdL, sys, vdw_result);
}

void grdv12(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, SOLV* sol, fftw_real v12_result[l][l][2*(l/2+1)]){
    double sys[4] = {1.0, sol->Cut, dx, 1.0};
    
    compute_grid(nCrd, Crds, nPro, Pros, l, 
                &cache.crowd_v12, &cache.v12_computed,
                r12GrdR, r12GrdL, sys, v12_result);
}

void setdtErn(int l, tErn softfb[l][l][l], bool vol[l][l][l], fftw_real soft[l][l][2*(l/2+1)]){
    const int h=l/2;

    #pragma omp parallel for
    for (int i=0;i<l;i++){
        int ii=(i+h)%l;
        for (int j=0;j<l;j++){
            int jj=(j+h)%l;
            for (int k=0;k<l;k++){
                int kk=(k+h)%l;
                if (vol[i][j][k]){
                    softfb[ii][jj][kk]=tErn_MAX;
                }else{
                    float fern=soft[i][j][k]*tErn_SCL;
                    fern=round(fern);
                    tErn ern=(tErn)(fern);
                    softfb[ii][jj][kk]=ern;
                }
            }
        }
    }
}

void softgrd(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, SOLV* sol, tErn softfb[l][l][l]){
    bool vol[l][l][l];
    fftw_real working[l][l][2*(l/2+1)];
    fftw_real lcEle[l][l][2*(l/2+1)];
    fftw_real lcVdw[l][l][2*(l/2+1)];
    fftw_real combined[l][l][2*(l/2+1)];
    
    static int init_lattice=0;
    if (init_lattice==0) init_lattice_data();

    SetZero(l,lcEle);
    SetZero(l,lcVdw);
    SetZero(l,combined);
    
    // Volume calculation and reporting
    grdvol(nCrd,Crds,nPro,Pros,l,dx,sol,working);
    vRep(l,vol,0.001,working,1.0,"vol",sol->kBT);
    fflush(stdout);
    
    // Electrostatic calculation and reporting
    grdele(nCrd,Crds,nPro,Pros,l,dx,sol,working);
    AddTo1filt(l,combined,1.0,working,sol->escl*332.0/sol->sdie,vol);
    AddTo1filt(l,lcEle,1.0,working,sol->escl*332.0/sol->sdie,vol);
    vRep(l,vol,0.001,working,sol->escl*332.0/sol->sdie,"ele",sol->kBT);
    fflush(stdout);
    
    // Van der Waals calculation and reporting
    grdvdw(nCrd,Crds,nPro,Pros,l,dx,sol,working);
    AddTo1filt(l,combined,1.0,working,sol->vscl,vol);
    AddTo1filt(l,lcVdw,1.0,working,sol->vscl,vol);
    vRep(l,vol,0.001,working,sol->vscl,"vdw",sol->kBT);
    fflush(stdout);
    
    // v12 calculation and reporting
    grdv12(nCrd,Crds,nPro,Pros,l,dx,sol,working);
    AddTo1filt(l,combined,1.0,working,sol->vscl,vol);
    AddTo1filt(l,lcVdw,1.0,working,sol->vscl,vol);
    vRep(l,vol,0.001,working,sol->vscl,"v12",sol->kBT);
    vRep(l,vol,0.001,combined,1.0,"v+e",sol->kBT);
    fflush(stdout);
    
    // Set energy array data into softfb
    setdtErn(l,softfb,vol,combined); 
}

void softgrd_cleanup(void) {
    if (cache.crowd_vol) { free(cache.crowd_vol); cache.crowd_vol = NULL; }
    if (cache.crowd_ele) { free(cache.crowd_ele); cache.crowd_ele = NULL; }
    if (cache.crowd_vdw) { free(cache.crowd_vdw); cache.crowd_vdw = NULL; }
    if (cache.crowd_v12) { free(cache.crowd_v12); cache.crowd_v12 = NULL; }
    cache.vol_computed = cache.ele_computed = cache.vdw_computed = cache.v12_computed = false;
    cache.cached_l = 0;
}
