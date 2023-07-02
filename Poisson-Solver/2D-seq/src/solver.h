/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __SOLVER_H_
#define __SOLVER_H_
#include "parameter.h"

typedef struct {
    double dx, dy;
    int imax, jmax;
    double *p, *rhs;
    double eps, omega, rho;
    int itermax;
} Solver;

extern void initSolver(Solver*, Parameter*, int problem);
extern void writeResult(Solver*);
extern void solve(Solver*);
extern void solveRB(Solver*);
extern void solveRBA(Solver*);
#endif
