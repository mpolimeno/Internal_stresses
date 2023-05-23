#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdbool>
#include <fstream> // To use ifstream and read files
#include "funcs.h"

// DRIVER
int main(int argc, char* argv[]) {
    
    // Test for a single cube
    int kDimension = 3;
    int kNumberOfFaces = 6;

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

    // Build SL Matrix
    double* SingleLayerMatrix = new double[kDimension*kDimension];
    for (int i=0;i<kNumberOfFaces;i++) {
        int* x_s          = new int[kDimension];
        int* current_face = new int[kDimension];
        
        for (int j=0;j<kDimension;j++) {
            *(x_s+j)          = *(evaluation_point+j);
            *(current_face+j) = *(center_of_face+i*kDimension+j);
        }

        std::cout << "Evaluation point is: " << "\n";
        for (int j=0;j<kDimension;j++) std::cout << *(x_s+j) << " ";
        std::cout << "\n";
        
        std::cout << "Current face is: " << "\n";
        for (int j=0;j<kDimension;j++) std::cout << *(current_face+j) << " ";
        std::cout << "\n";
        
        int n_hat = *(direction_of_normal+i);
        
        for (int j=0;j<kDimension;j++) {
            for (int k=0;k<kDimension;k++) *(SingleLayerMatrix+j*kDimension+k) = 0.;
        }
        BuildMatrixForSingleLayerPotential(current_face,kNumberOfFaces,x_s,kDimension,n_hat,SingleLayerMatrix);
        
        for (int m=0;m<kDimension;m++) {
            for (int n=0;n<kDimension;n++) std::cout << *(SingleLayerMatrix+m*kDimension+n) << " ";
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    /*
    // Display Face Centers
    for (int i=0;i<kNumberOfFaces;i++) {
        for (int j=0;j<kDimension;j++) std::cout << *(center_of_face+i*kDimension+j) << " ";
        std::cout << "\n";
    }
    std::cout << "\n";

    // Display Evaluation Point
    for (int i=0;i<kDimension;i++) std::cout << *(evaluation_point+i) << " ";
    std::cout << "\n";
    
    std::cout << "\n";
    
    // Display Normal Directions
    for (int i=0;i<kNumberOfFaces;i++) std::cout << *(direction_of_normal+i) << "\n";
    std::cout << "\n";
    */

    // Deallocate memory
    delete[] center_of_face;
    delete[] evaluation_point;
    delete[] direction_of_normal;
    delete[] SingleLayerMatrix;

    return 0;
}
