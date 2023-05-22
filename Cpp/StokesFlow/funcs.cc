#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdbool>
#include "funcs.h"


// Function that builds the matrix for the single-layer potential
// See Reference https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.044305 for details

double ComputeEuclideanNorm(int* vector_x, int* vector_y) {
    
    int DIM = 3;
    double* tmp = new double[DIM];
    for (int i=0;i<DIM;i++) {
        *(tmp+i)  = *(vector_x+i) - *(vector_y+i);
        *(tmp+i) *= *(tmp+i);
    }
    
    double norm_squared = 0;
    for (int i=0;i<DIM;i++) {
        norm_squared += *(tmp+i);
    }

    return norm_squared;
}

void BuildMatrixForSingleLayerPotential(int* center_of_face, int kNumberOfFaces, int* evaluation_point, int kDimension, int normal_direction, double* SingleLayerMatrix) {

    //double* current_position = new double[kDimension];
    //for (int i=0;i,kDimension;i++) *(current_position+i) = *(evaluation_point+i) - *(center_of_face+i);

    double norm_squared = ComputeEuclideanNorm(evaluation_point,center_of_face);

    double* constant_ij      = new double[kDimension*kDimension];
    double* xx_ij            = new double[kDimension*kDimension];
    for (int i=0;i<kDimension;i++) {
        for (int j=0;j<kDimension;j++) {
             if (norm_squared==0.) {
                if (i==j) {
                    *(constant_ij+i*kDimension+j) = 8.*asinh(1.);
                    *(xx_ij+i*kDimension+j) = (i!=normal_direction) ? 4.*asinh(1.) : 0.;
                } else {
                    *(constant_ij+i*kDimension+j) = 0.;
                }
            }
        }
    }
    
    for (int i=0;i<kDimension;i++) {
        for (int j=0;j<kDimension;j++) {
            *(SingleLayerMatrix+i*kDimension+j) = *(constant_ij+i*kDimension+j) + *(xx_ij+i*kDimension+j);
        }
    }

/*
    for (int i=0;i<kDimension;i++) {
        for (int j=0;j<kDimension;j++) {
            *(SingleLayerMatrix+i*kDimension+j) = 0.;
        }
    }
    */
    return;
}
