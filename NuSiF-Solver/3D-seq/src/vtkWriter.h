/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __VTKWRITER_H_
#define __VTKWRITER_H_
#include <stdio.h>

#include "grid.h"

typedef enum VtkFormat { ASCII = 0, BINARY } VtkFormat;

typedef struct VtkOptions {
    VtkFormat fmt;
    Grid grid;
    FILE* fh;
} VtkOptions;

typedef struct VtkVector {
    double *u, *v, *w;
} VtkVector;

extern void vtkOpen(VtkOptions* opts, char* filename);
extern void vtkVector(VtkOptions* opts, char* name, VtkVector vec);
extern void vtkScalar(VtkOptions* opts, char* name, double* p);
extern void vtkClose(VtkOptions* opts);
#endif // __VTKWRITER_H_
