/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "likwid-marker.h"
#include "parameter.h"
#include "solver.h"
#include "timing.h"

#define LIKWID_PROFILE(tag, call)                                                        \
    startTime = getTimeStamp();                                                          \
    LIKWID_MARKER_START(#tag);                                                           \
    call(&solver);                                                                       \
    LIKWID_MARKER_STOP(#tag);                                                            \
    endTime = getTimeStamp();

enum VARIANT { SOR = 1, RB, RBA };

int main(int argc, char** argv)
{
    int variant = SOR;
    double startTime, endTime;
    Parameter params;
    Solver solver;
    initParameter(&params);
    LIKWID_MARKER_INIT;

    if (argc < 2) {
        printf("Usage: %s <configFile>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }

    readParameter(&params, argv[1]);
    // printParameter(&params);
    if (argc == 3) {
        variant = atoi(argv[2]);
    }
    if (argc == 4) {
        sscanf("%lf", argv[2], &params.omg);
        sscanf("%lf", argv[2], &params.rho);
    }

    initSolver(&solver, &params, 2);
    switch (variant) {
    case SOR:
        printf("Plain SOR\n");
        LIKWID_PROFILE("SOR", solve);
        break;
    case RB:
        printf("Red-black SOR\n");
        LIKWID_PROFILE("RB", solveRB);
        break;
    case RBA:
        printf("Red-black SOR with acceleration\n");
        LIKWID_PROFILE("RBA", solveRBA);
        break;
    }
    printf(" %.2fs\n", endTime - startTime);
    // Commented out since results are not required for the benchmark.
    //writeResult(&solver);

    LIKWID_MARKER_CLOSE;
    return EXIT_SUCCESS;
}
