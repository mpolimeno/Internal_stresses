#ifndef FUNCS_H_INCLUDED
#define FUNCS_H_INCLUDED

extern void BuildMatrixForSingleLayerPotential(int* center_of_face, int kNumberOfFaces, int* evaluation_point, int kDimension, int normal_direction, double* SingleLayerMatrix);

extern double ComputeEuclideanNorm(int* vector_x, int* vector_y);

#endif
