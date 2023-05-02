#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdbool>
#include "randomwalk.h"

// This class contains relevant parameters and member functions to build the fractal aggregate
class Aggregate {
    private:
        const int kDimension_;               // This is the dimensionality of the system. It does not change, so its type is const int
        const int kTotalNumberOfCubes_;      // This is the user-input total number of cubes. Does not change, so its type is const int
        int number_of_cubes_in_aggregate_;   // This is the number of cubes in the aggregate at a specific moment. It changes but cannot be greater than the total number of cubes in the system
        double* center_of_mass_;            // Center of Mass of the aggregate

    // Define some member functions to access the data from the class and set parameters and variables
    public:
        // Constructor
        // Note that 'const int' types must be initialized using constructor initializer list
        Aggregate(const int N, int M, const int DIM, double* com) : kTotalNumberOfCubes_(N), kDimension_(DIM) {
            number_of_cubes_in_aggregate_ = M;
            center_of_mass_ = com;
        };
        // Destructor
        ~Aggregate() {};
         
        // Initialize Center of Mass to [0,0,0]
        void SetCenterOfMass(const int dimensionality, double* com) {
            for (int d=0;d<dimensionality;d++) *(com+d) = 0.;
        }

        // Declaration of function to compute final center of mass of aggregate
        void ComputeFinalCenterOfMass(const int kDimension_, const int kTotalNumberOfCubes_, int* aggregate_position_, double* center_of_mass_);
        
        void ComputeCenterOfMass(const int kDimension_, int number_of_cubes_in_aggregate_, int* aggregate_position_, double* center_of_mass_);
        // For later
        // Declaration of member functions for radius of gyration of aggregate
        //void ComputeRadiusOfGyration(int DIM, int total_number_of_cubes, int* center_of_mass, int* aggregate_position);

        // Declaration of member functions for maximum radius of aggregate
        //void ComputeMaximumRadius(int DIM, int total_number_of_cubes, int* center_of_mass, int* aggregate_position);
        
};

void Aggregate::ComputeFinalCenterOfMass(const int kDimension, const int N, int* aggregate_position, double* center_of_mass) {
    for (int i=0;i<N;i++) {
        for (int d=0;d<kDimension;d++) *(center_of_mass+d) += *(aggregate_position+i*kDimension+d);
    }

    for (int d=0;d<kDimension;d++) *(center_of_mass+d) /= N;  
}

void Aggregate::ComputeCenterOfMass(const int kDimension, int M, int* aggregate_position, double* center_of_mass) {
    for (int i=0;i<M;i++) {
        for (int d=0;d<kDimension;d++) *(center_of_mass+d) += *(aggregate_position+i*kDimension+d);
    }

    for (int d=0;d<kDimension;d++) *(center_of_mass+d) /= M;  
}

// DRIVER
int main(int argc,char* argv[]) {
    
    // Check step-size restrictions
    if (DX!=DY||DX!=DZ||DY!=DZ) {
        printf("ERROR: currently only DX=DY=DZ supported\n");
        exit(0);
    }
    
    // Set parameters
    const int kDimension = 3;
    // Check dimensionality restrictions
    if (kDimension!=3) {
        printf("ERROR: currently only kDimension=3 supported\n");
        exit(0);
    }
    const int kNumberOfCubes = 15;
    int CubesInOneAggregate;
    
    // Define variable for center of mass of aggregate
    double* com = new double[kDimension];
    
    // First instance of the class
    Aggregate agg(kNumberOfCubes,CubesInOneAggregate,kDimension,com);

    // Set center of mass of aggregate
    agg.SetCenterOfMass(kDimension,com);
    
    // Deal with args
    for (int i=0;i<argc;i++) printf("** argv[%d]: %s\n",i,argv[i]);

    if (argc!=3) {
        printf("USAGE: (executable from this program) (seed number) (output filename)\n");
        exit(0);
    }

    int seed = atoi(argv[1]);
    char outfile[128];

    strcpy(outfile,argv[2]);

    printf("** seed = %d\n",seed);
    printf("** output = %s\n",outfile);

    // Intialize random number generator
    srand(seed);

    // Note: all arrays are allocated dynamically using row-wise mapping (standard in C++)

    // Array containing the position of aggregate (integer values)
    int* aggregate_position = new int[kNumberOfCubes*kDimension];

    // Initialize position of the first cube (first member of the aggregate)
    // We set the first cube to be at the origin and it is not allowed to move
    for (int d=0;d<kDimension;d++) *(aggregate_position+0+d) = 0;

    // Call function for DLA routine to build aggregate
    BuildAggregate(kNumberOfCubes,kDimension,aggregate_position); 
    //BuildAggregate(kNumberOfCubes,kDimension,aggregate_position,com); 

    // Compute center of mass of aggregate and print it out
    agg.ComputeFinalCenterOfMass(kDimension,kNumberOfCubes,aggregate_position,com);
    for (int d=0;d<kDimension;d++) std::cout << *(com+d) << " ";
    std::cout << "\n";

    // Print position of each cube in the aggregate to a file
    FILE *fp = fopen(outfile,"w+");
    PrintAggregatePosition(fp,kNumberOfCubes,kDimension,aggregate_position);
    fclose(fp);

    // Deallocate memory of array that stores final position of cubes in aggregate
    delete[] aggregate_position;

    return 0;
}
