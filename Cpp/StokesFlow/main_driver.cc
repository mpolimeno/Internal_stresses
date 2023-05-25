#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdbool>
#include <fstream> // To use ifstream and read files
#include <chrono>
#include "funcs.h"

using namespace std::chrono;

// DRIVER
int main(int argc, char* argv[]) {
    
    int kDimension = 3;
    // Number of external faces
    int kNumberOfFaces = 10;

    int* center_of_face      = new int[kNumberOfFaces*kDimension];
    int* evaluation_point    = new int[kDimension];
    int* direction_of_normal = new int[kNumberOfFaces];
    
    // For faces center array
    int kFaceCentersArraySize = kNumberOfFaces*kDimension; 
    int count_face = 0;
    std::ifstream FaceCenters;
    FaceCenters.open("FaceCenters.txt");

    // Read the numbers from the file into the array
    while (count_face<kFaceCentersArraySize && FaceCenters >> *(center_of_face+count_face)) count_face++;
    FaceCenters.close();

    // For integration point array
    int count_integration = 0;
    std::ifstream IntegrationPoint;
    IntegrationPoint.open("IntegrationPoint.txt");
    
    while (count_integration<kDimension && IntegrationPoint >> *(evaluation_point+count_integration)) count_integration++;
    IntegrationPoint.close();

    // For direction_of_normal array
    int count_direction = 0;
    std::ifstream NormalDirection;
    NormalDirection.open("NormalDirection.txt");
    
    while (count_direction<kNumberOfFaces && NormalDirection >> *(direction_of_normal+count_direction)) count_direction++;
    NormalDirection.close();

    int face_index = kNumberOfFaces-1;
    
    // Build LHS of Linear System
    auto start = high_resolution_clock::now();
    double* LeftHandSide = new double[kDimension*kNumberOfFaces*kDimension*kNumberOfFaces]; // LHS is 3*Number_of_Faces x 3*Number_of_Faces
    for (int i=0;i<kNumberOfFaces;i++) {
        for (int j=0;j<kNumberOfFaces;j++) {
            for (int m=0;m<kDimension;m++) {
                for (int n=0;n<kDimension;n++) {
                    *(LeftHandSide+i*kNumberOfFaces*face_index+j*face_index+m*kDimension+n) = 0.;
                }
            }
        }
    }

    int count = 1;
    while (count<=10) {
        for (int i=0;i<kNumberOfFaces;i++) {
            int* x_s = new int[kDimension];
            for (int j=0;j<kDimension;j++) {
                *(x_s+j) = *(center_of_face+i*kDimension+j); // evaluation point
            }
        
            for (int f=0;f<kNumberOfFaces;f++) {

                int* current_face = new int[kDimension];
                for (int g=0;g<kDimension;g++) *(current_face+g) = *(center_of_face+f*kDimension+g);
                int n_hat = *(direction_of_normal+f);
        
                // Use Analytical Expressions for Single Layer Potential to build the matrix entries
                double* SingleLayerMatrix = new double[kDimension*kDimension];
                for (int j=0;j<kDimension;j++) {
                    for (int k=0;k<kDimension;k++) *(SingleLayerMatrix+j*kDimension+k) = 0.;
                }
                BuildMatrixForSingleLayerPotential(current_face,kNumberOfFaces,x_s,kDimension,n_hat,SingleLayerMatrix);
        
                delete[] current_face;
                 
                for (int j=0;j<kDimension;j++) {
                    for (int k=0;k<kDimension;k++) {
                        *(LeftHandSide+i*kNumberOfFaces*face_index+f*face_index+j*kDimension+k) = *(SingleLayerMatrix+j*kDimension+k);
                    }
                }

                delete[] SingleLayerMatrix;
            }
            delete[] x_s;
        }
        count++;
    }
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds> (stop - start);

    std::cout << "Time to build Single Layer Matrices: " << duration.count()  << " microseconds" << std::endl;
    
    
    std::cout << "Left-Hand Side is " << std::endl;
    std::cout << "\n";
    for (int i=0;i<kNumberOfFaces;i++) {
        for (int j=0;j<kNumberOfFaces;j++) {
            for (int m=0;m<kDimension;m++) {
                for (int n=0;n<kDimension;n++) {
                    std::cout << *(LeftHandSide+i*kNumberOfFaces*9+j*9+m*kDimension+n) << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Deallocate memory
    delete[] LeftHandSide;
    delete[] center_of_face;
    delete[] evaluation_point;
    delete[] direction_of_normal;

    return 0;
}
