// python3 json-to-header-converter.py cromer_mann_4gauss.json vdW.json electrons.json cromer_mann_4gauss.h

#ifndef CROMER_MANN_4GAUSS_H
#define CROMER_MANN_4GAUSS_H
#include <stddef.h>

// Cromer-Mann 4-Gauss parameters structure
// Formula: f(s) = a[0]*exp(-b[0]*s^2) + a[1]*exp(-b[1]*s^2) + a[2]*exp(-b[2]*s^2) + a[3]*exp(-b[3]*s^2)
typedef struct {
    double a[4];      // Amplitude parameters
    double b[4];      // Decay parameters
    double vdw_radius; // van der Waals radius (Angstroms)
    double scl;        //scaling for uniqe volume
    int electrons;    // Number of electrons
} CromerMann4Gauss;

// Element/Molecule identifiers
enum ElementType {
    ELEMENT_UNKNOWN = 0,
    ELEMENT_HOH = 1,
    ELEMENT_H = 2,
    ELEMENT_C = 3,
    ELEMENT_N = 4,
    ELEMENT_O = 5,
    ELEMENT_P = 6,
    ELEMENT_S = 7,
};

// Static parameter definitions
static const CromerMann4Gauss CROMER_MANN_HOH = {
    .a = {5.9565474326684456, -1.0902225877160387, 4.455816562289012, -9.253357957505079},
    .b = {18.395696887124288, 1.4191113040744165, 1.4200667487158216, 69.96135092707375},
    .vdw_radius = 0.0,
    .scl = 0.0,
    .electrons = 0
};

static const CromerMann4Gauss CROMER_MANN_H = {
    .a = {0.48964434865164597, 0.25795506418162684, 0.19862808650933614, 0.05398838736648019},
    .b = {20.659384838283998, 7.740482800259867, 49.551905131184874, 2.2010683273966882},
    .vdw_radius = 1.2,
    .scl = 1.33418,
    .electrons = 1
};

static const CromerMann4Gauss CROMER_MANN_C = {
    .a = {3.049502465329157, 0.9469490189353504, 0.9205181585530644, 1.0895088896753868},
    .b = {16.122830991045564, 0.879328843206143, 0.2466316869059288, 51.88436420676904},
    .vdw_radius = 1.7,
    .scl =1.18354,
    .electrons = 6
};

static const CromerMann4Gauss CROMER_MANN_N = {
    .a = {5.551459140691302, 2.480555683008789, 2.773120764746902, -3.8187500178972695},
    .b = {0.3145270000398675, 7.917446957141198, 24.43261890668131, 0.3145232450476127},
    .vdw_radius = 1.55,
    .scl = 1.04985,
    .electrons = 7
};

static const CromerMann4Gauss CROMER_MANN_O = {
    .a = {3.2844165713667874, 2.1863494590744375, 1.734158144526276, 0.7918335350536696},
    .b = {13.427747243117944, 5.212997411416355, 0.24282907358098274, 32.87575302824983},
    .vdw_radius = 1.52,
    .scl = 9.93450e-01,
    .electrons = 8
};

static const CromerMann4Gauss CROMER_MANN_P = {
    .a = {7.138789042249852, 4.142606814307699, 2.128293691849062, 1.5929377824602673},
    .b = {1.7548326959642166, 26.16163488546785, 0.12230615790259426, 68.26204711582082},
    .vdw_radius = 1.9,
    .scl = 1.0,
    .electrons = 15
};

static const CromerMann4Gauss CROMER_MANN_S = {
    .a = {7.1088513543231375, 5.215028849933757, 2.091806378705691, 1.583611951864076},
    .b = {1.4377552638137123, 22.18043778213814, 0.10418335440325559, 56.17384058069765},
    .vdw_radius = 1.8,
    .scl = 1.0,
    .electrons = 16
};

// Get Cromer-Mann parameters for a specific element type
static inline const CromerMann4Gauss* get_cromer_mann_params(enum ElementType element_type) {
    switch (element_type) {
    case ELEMENT_HOH:
        return &CROMER_MANN_HOH;
    case ELEMENT_H:
        return &CROMER_MANN_H;
    case ELEMENT_C:
        return &CROMER_MANN_C;
    case ELEMENT_N:
        return &CROMER_MANN_N;
    case ELEMENT_O:
        return &CROMER_MANN_O;
    case ELEMENT_P:
        return &CROMER_MANN_P;
    case ELEMENT_S:
        return &CROMER_MANN_S;
    default:
        return NULL;
    }
}

#endif /* CROMER_MANN_4GAUSS_H */
