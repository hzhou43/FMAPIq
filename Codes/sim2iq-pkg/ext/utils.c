#include <stdio.h>
#include <stdlib.h>
#include "saxs_calculator.h"

void* check_malloc(size_t size) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

QIData qi_data_create(size_t initial_capacity) {
    return (QIData){
        .data = check_malloc(initial_capacity * sizeof(QIPoint)),
        .count = 0,
        .capacity = initial_capacity
    };
}

Structure structure_create(void) {
    return (Structure){
        .atoms = check_malloc(INITIAL_CAPACITY * sizeof(Atom)),
        .count = 0,
        .capacity = INITIAL_CAPACITY
    };
}

void structure_ensure_capacity(Structure* structure) {
    if (structure->count >= structure->capacity) {
        structure->capacity *= 2;
        Atom* new_atoms = realloc(structure->atoms, structure->capacity * sizeof(Atom));
        if (!new_atoms) {
            free(structure->atoms);
            fprintf(stderr, "Memory reallocation failed\n");
            exit(EXIT_FAILURE);
        }
        structure->atoms = new_atoms;
    }
}

void structure_free(Structure* structure) {
    if (structure) {
        free(structure->atoms);
        structure->atoms = NULL;
        structure->count = 0;
        structure->capacity = 0;
    }
}

void qidata_free(QIData* data) {
    if (data) {
        free(data->data);
        data->data = NULL;
        data->count = 0;
        data->capacity = 0;
    }
}