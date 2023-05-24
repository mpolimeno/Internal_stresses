#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdbool>
#include "funcs.h"


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
    delete[] tmp;

    return norm_squared;
}

// Function that builds the matrix for the single-layer potential
// See Reference https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.044305 for details

void BuildMatrixForSingleLayerPotential(int* center_of_face, int kNumberOfFaces, int* evaluation_point, int kDimension, int normal_direction, double* SingleLayerMatrix) {

    double* current_position = new double[kDimension];
    for (int i=0;i<kDimension;i++) *(current_position+i) = *(evaluation_point+i) - *(center_of_face+i);

    double norm_squared = ComputeEuclideanNorm(evaluation_point,center_of_face);

    std::cout << "Squared Norm is: " << norm_squared << std::endl;
    std::cout << "\n";

    double* constant_ij = new double[kDimension*kDimension];
    double* xx_ij       = new double[kDimension*kDimension];
    for (int i=0;i<kDimension;i++) {
        for (int j=0;j<kDimension;j++) {
             if (rint(norm_squared)==0) { // We are at the singularity: center_of_face = evaluation_point
                if (i==j) {
                    *(constant_ij+i*kDimension+j) = 8.*asinh(1.);
                    *(xx_ij+i*kDimension+j)       = (i!=(normal_direction-1)) ? 4.*asinh(1.) : 0.;
                } else {
                    *(constant_ij+i*kDimension+j) = 0.;
                    *(xx_ij+i*kDimension+j)       = 0.;
                }
            } else { // We are not at the singularity
                if (i==j) {
                    int xs_index = ((normal_direction+1) % 3);
                    int ys_index = ((normal_direction+2) % 3);
                    xs_index = (xs_index==0) ? 3 : xs_index;
                    ys_index = (ys_index==0) ? 3 : ys_index;
                    
                    xs_index = xs_index - 1;
                    ys_index = ys_index - 1;

                    double x_s = *(current_position+xs_index);
                    double y_s = *(current_position+ys_index);
                    int n_dir = normal_direction - 1; // Indexing starts at 0
                    double z_s = *(current_position+n_dir);

                    std::cout << x_s << std::endl;
                    std::cout << y_s << std::endl;
                    std::cout << z_s << std::endl;

                    double z = 0.;

                    // Declaring them here does not seem like a good idea, but for now it is okay
                    double p1;
                    double p2;
                    double p3;
                    double p4;

                    int flag = 0;
                    if (rint(z)==rint(z_s) && rint(abs(x_s))==1 && rint(abs(y_s))==1) {
                        flag = 1;

                        p1 = 0.;
                        p2 = 0.;
                        p3 = 0.;
                        p4 = 0.;
                        
                        *(constant_ij+i*kDimension+j) = 4.*asinh(1.);

                    } else {
                        if (rint(abs(x_s))!=1 && rint(abs(y_s))!=1) {
                            if (rint(z)==rint(z_s)) {
                                p1 = (1-y_s)  * log(sqrt((1-x_s)*(1-x_s)   + (1-y_s)*(1-y_s))   + (1-x_s))  + (1-x_s)  * log(sqrt((1-x_s)*(1-x_s)   + (1-y_s)*(1-y_s))   + (1-y_s))  -  (1-y_s);
                                p2 = (-1-y_s) * log(sqrt((1-x_s)*(1-x_s)   + (-1-y_s)*(-1-y_s)) + (1-x_s))  + (1-x_s)  * log(sqrt((1-x_s)*(1-x_s)   + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s);
                                p3 = (1-y_s)  * log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s))   + (-1-x_s)) + (-1-x_s) * log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s))   + (1-y_s))  -  (1-y_s);
                                p4 = (-1-y_s) * log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-x_s)) + (-1-x_s) * log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s);
                            } else {
                                p1 = - (z-z_s)*atan(((1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)    + (1-x_s)*(1-x_s))     + (1-y_s)*(1-y_s)))) 
                                     + (z-z_s)*atan((1-y_s)/(z-z_s))  
                                     + (1-y_s)*log((1-x_s)   + sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)   + (1-y_s)*(1-y_s)))    + (1-x_s)*log(sqrt((z-z_s)*(z-z_s)  + (1-x_s)*(1-x_s)   + (1-y_s)*(1-y_s))   + (1-y_s))  -  (1-y_s);
                                p2 = - (z-z_s)*atan(((1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)   + (1-x_s)*(1-x_s))     + (-1-y_s)*(-1-y_s)))) 
                                     + (z-z_s)*atan((-1-y_s)/(z-z_s)) 
                                     + (-1-y_s)*log((1-x_s)  + sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)   + (-1-y_s)*(-1-y_s)))  + (1-x_s)*log(sqrt((z-z_s)*(z-z_s)  + (1-x_s)*(1-x_s)   + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s);
                                p3 = - (z-z_s)*atan(((-1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)   + (-1-x_s)*(-1-x_s))   + (1-y_s)*(1-y_s))))   
                                     + (z-z_s)*atan((1-y_s)/(z-z_s))  
                                     + (1-y_s)*log((-1-x_s)  + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)))    + (-1-x_s)*log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s))   + (1-y_s))  -  (1-y_s);
                                p4 = - (z-z_s)*atan(((-1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)  + (-1-x_s)*(-1-x_s))   + (-1-y_s)*(-1-y_s)))) 
                                     + (z-z_s)*atan((-1-y_s)/(z-z_s)) 
                                     + (-1-y_s)*log((-1-x_s) + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)))  + (-1-x_s)*log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s);
                            }
                        } else {
                            if (rint(z)!=rint(z_s)) {
                                p1 = -(z-z_s)*atan(((1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)   + (1-x_s)*(1-x_s))   + (1-y_s)*(1-y_s))))   + (z-z_s)*atan((1-y_s)/(z-z_s))  + (1-y_s)*log((1-x_s)   + sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)   + (1-y_s)*(1-y_s)))   + (1-x_s)*log(sqrt((z-z_s)*(z-z_s)  + (1-x_s)*(1-x_s)   + (1-y_s)*(1-y_s))   + (1-y_s))  -  (1-y_s);
                                p2 = -(z-z_s)*atan(((1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)  + (1-x_s)*(1-x_s))   + (-1-y_s)*(-1-y_s)))) + (z-z_s)*atan((-1-y_s)/(z-z_s)) + (-1-y_s)*log((1-x_s)  + sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)   + (-1-y_s)*(-1-y_s))) + (1-x_s)*log(sqrt((z-z_s)*(z-z_s)  + (1-x_s)*(1-x_s)   + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s);
                                p3 = -(z-z_s)*atan(((-1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)  + (-1-x_s)*(-1-x_s)) + (1-y_s)*(1-y_s))))   + (z-z_s)*atan((1-y_s)/(z-z_s))  + (1-y_s)*log((-1-x_s)  + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)))   + (-1-x_s)*log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s))   + (1-y_s))  -  (1-y_s);
                                p4 = -(z-z_s)*atan(((-1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s)) + (-1-y_s)*(-1-y_s)))) + (z-z_s)*atan((-1-y_s)/(z-z_s)) + (-1-y_s)*log((-1-x_s) + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s))) + (-1-x_s)*log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s);
                            } else {
                                if (rint(abs(x_s))==1 && rint(abs(y_s))!=1) {
                                    p1 = -((1./2.)*(1-y_s)*(log((1-y_s)*(1-y_s)) - 2));
                                    p2 = -((1-y_s)*(log(sqrt((1-y_s)*(1-y_s) + 4) + 2) - 1) + 2*asinh((1-y_s)/2.));
                                    p3 = -((1./2.) *(-1-y_s)*(log((-1-y_s)*(-1-y_s)) - 2));
                                    p4 = -((-1-y_s)*(log(sqrt((-1-y_s)*(-1-y_s) + 4) + 2) - 1) + 2*asinh((-1-y_s)/2.));
                                }
                                if (rint(abs(x_s))!=1 && rint(abs(y_s))==1) {
                                    double r_s = x_s; // swap them
                                    x_s = y_s;
                                    y_s = r_s;
                                    
                                    p1 = -((1./2.)*(1-y_s)*(log((1-y_s)*(1-y_s)) - 2));
                                    p2 = -((1-y_s)*(log(sqrt((1-y_s)*(1-y_s) + 4) + 2) - 1) + 2*asinh((1-y_s)/2.));
                                    p3 = -((1./2.) *(-1-y_s)*(log((-1-y_s)*(-1-y_s)) - 2));
                                    p4 = -((-1-y_s)*(log(sqrt((-1-y_s)*(-1-y_s) + 4) + 2) - 1) + 2*asinh((-1-y_s)/2.));
                                    
                                    r_s = x_s; // unswap them
                                    x_s = y_s;
                                    y_s = r_s;
                                }
                            }
                        }

                        *(constant_ij+i*kDimension+j) = p1 - p2 - p3 + p4;
                    }

                } else { // i!=j, i.e. we are off the diagonal
                    *(constant_ij+i*kDimension+j) = 0.;
                }
            }
            // #################################################### //
            // Now we assign the analytical values to the xx terms
            double px1;
            double px2;
            double px3;
            double px4;
            
            if (i==j && i==(normal_direction-1)) {
                int xs_index = ((normal_direction+1) % 3);
                int ys_index = ((normal_direction+2) % 3);
                xs_index = (xs_index==0) ? 3 : xs_index;
                ys_index = (ys_index==0) ? 3 : ys_index;
                
                xs_index = xs_index - 1;
                ys_index = ys_index - 1;

                double x_s = *(current_position+xs_index);
                double y_s = *(current_position+ys_index);
                int n_dir = normal_direction - 1; // Indexing starts at 0
                double z_s = *(current_position+n_dir);

                double z = 0.;

                if (rint(z)!=rint(z_s)) {
                    if (rint(x_s)!=1) {
                        px1 = (1-x_s)*(atan(((1-y_s)*(1-x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)))))/((z-z_s)*(1-x_s));
                        px2 = (1-x_s)*(atan(((1+y_s)*(1-x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1+y_s)*(1+y_s)))))/((z-z_s)*(1-x_s));    
                    } else {
                        px1 = 0.;
                        px2 = 0.;
                    }
                    if (rint(x_s)!=(-1)) {
                        px3 = (1+x_s)*(atan(((1-y_s)*(1+x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1-y_s)*(1-y_s)))))/((z-z_s)*(1+x_s));
                        px4 = (1+x_s)*(atan(((1+y_s)*(1+x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1+y_s)*(1+y_s)))))/((z-z_s)*(1+x_s));
                    } else {
                        px3 = 0.;
                        px4 = 0.;
                    }
                } else {
                    px1 = 0.;
                    px2 = 0.;
                    px3 = 0.;
                    px4 = 0.;
                }
                *(xx_ij+i*kDimension+j) = (z-z_s)*(z-z_s)*(px1+px2+px3+px4); 
 
            }
            //std::cout << "Current Position is: " << "\n";
            //for (int k=0;k<kDimension;k++) std::cout << *(current_position+k) << " ";
            //std::cout << "\n";
            if (i==j && i!=(normal_direction-1)) {
                double x_s = *(current_position+i);
                int ys_index = 5-normal_direction-1-i;
                double y_s = *(current_position+ys_index);
                int n_dir = normal_direction - 1;
                double z_s = *(current_position+n_dir);

                double z = 0.;

                if (rint(z)!=rint(z_s)) {
                px1 = - (z-z_s)*atan(((1-x_s)*(1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) +(1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)))) 
                      + (z-z_s)*atan((1-y_s)/(z-z_s)) + (1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)) +  (1-x_s)) - 1);
                px2 = - (z-z_s)*atan(((1-x_s)*(-1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)))) 
                      + (z-z_s)*atan((-1-y_s)/(z-z_s)) + (-1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) +  (1-x_s)) - 1);
                px3 = - (z-z_s)*atan(((-1-x_s)*(1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)))) 
                      + (z-z_s)*atan((1-y_s)/(z-z_s)) + (1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) +  (-1-x_s)) - 1);
                px4 = - (z-z_s)*atan(((-1-x_s)*(-1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)))) 
                      + (z-z_s)*atan((-1-y_s)/(z-z_s)) + (-1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) +  (-1-x_s)) - 1);
                } else {
                    if (rint(x_s)==1) {
                        if (rint(y_s)!=1) {
                            px1 = 0.5*(1-y_s)*(log((1-y_s)*(1-y_s)) - 2);
                        } else {
                            px1 = 0.;
                        }
                        if (rint(y_s)!=(-1)) {
                            px2 = 0.5*(-1-y_s)*(log((-1-y_s)*(-1-y_s)) - 2);
                        } else {
                            px2 = 0.;
                        }
                    } else {
                        if (rint(y_s)!=1) {
                            px1 = (1-y_s)*(log(sqrt( (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)) +  (1-x_s)) - 1);
                        } else {
                            px1 = 0.;
                        }
                        if (rint(y_s)!=(-1)) {
                            px2 =  (-1-y_s)*(log(sqrt((1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) +  (1-x_s)) - 1);
                        } else {
                            px2 = 0.;
                        }
                    }
                    if (rint(x_s)==(-1)) {
                        if (rint(y_s)!=1) {
                            px3 = (1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) +  (-1-x_s)) - 1);
                        } else {
                            px3 = 0.;
                        }
                        if (rint(y_s)!=(-1)) {
                            px4 = (-1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) +  (-1-x_s)) - 1);
                        } else {
                            px4 = 0.;
                        }
                    } else {
                        if (rint(y_s)!=1) {
                            px3 = (1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) +  (-1-x_s)) - 1);
                        } else {
                            px3 = 0.;
                        }
                        if (rint(y_s)!=(-1)) {
                            px4 = (-1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) +  (-1-x_s)) - 1);
                        } else {
                            px4 = 0.;
                        }
                    }
                }
                *(xx_ij+i*kDimension+j) = px1 - px2 - px3 + px4;
            }
            if (i!=j && i!=(normal_direction-1) && j==(normal_direction-1)) {
                double x_s = *(current_position+i);
                int ys_index = 5-normal_direction-1-i;
                double y_s = *(current_position+ys_index);
                int n_dir = normal_direction - 1;
                double z_s = *(current_position+n_dir);

                double z = 0.;

                if (rint(z)!=rint(z_s)) {
                px1 =  -log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (1-y_s)*(1-y_s)) + (1-y_s));
                px2 =  log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (1-y_s)*(1-y_s)) + (1-y_s));
                px3 =  log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s));
                px4 =  -log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s));

                } else {
                    px1 = 0.;
                    px2 = 0.;
                    px3 = 0.;
                    px4 = 0.;
                }
                *(xx_ij+i*kDimension+j) = (z-z_s)*(px1 + px2 + px3 + px4);
            }
            if (i!=j && j!=(normal_direction-1) && i==(normal_direction-1)) {
                double x_s = *(current_position+j);
                int ys_index = 5-normal_direction-1-j;
                double y_s = *(current_position+ys_index);
                int n_dir = normal_direction - 1;
                double z_s = *(current_position+n_dir);

                double z = 0.;

                if (rint(z)!=rint(z_s)) {
                px1 =  -log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (1-y_s)*(1-y_s)) + (1-y_s));
                px2 =  log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (1-y_s)*(1-y_s)) + (1-y_s));
                px3 =  log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s));
                px4 =  -log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s)); 
                } else {
                    px1 = 0.;
                    px2 = 0.;
                    px3 = 0.;
                    px4 = 0.;
                }
                *(xx_ij+i*kDimension+j) = (z-z_s)*(px1 + px2 + px3 + px4);
            }
            if (i!=j && j!=(normal_direction-1) && i!=(normal_direction-1)) {
                double x_s = *(current_position+i);
                double y_s = *(current_position+j);
                int n_dir = normal_direction - 1;
                double z_s = *(current_position+n_dir);

                double z = 0.;

                px1 = -sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s));
                px2 = sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1+y_s)*(1+y_s));
                px3 = sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1-y_s)*(1-y_s));
                px4 = -sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1+y_s)*(1+y_s));
                
                *(xx_ij+i*kDimension+j) = px1 + px2 + px3 + px4;
            }
        }
    }
    
    for (int i=0;i<kDimension;i++) {
        for (int j=0;j<kDimension;j++) {
            *(SingleLayerMatrix+i*kDimension+j) = *(constant_ij+i*kDimension+j) + *(xx_ij+i*kDimension+j);
        }
    }

    return;
}
