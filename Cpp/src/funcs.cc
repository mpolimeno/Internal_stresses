#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdbool>
#include "randomwalk.h"

// This code generates a fractal aggregate whose individual particles are cubes
// The aggregate is generated using the Diffusion-Limited Aggregation (DLA) model, following the Individually-Added Aggregate (IAA) approach
// See Reference https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.044305 for details

// This function builds an aggregate using the DLA routine
void BuildAggregate(const int kTotalNumberOfCubes, const int kDimension, int* aggregate_position) {
//void BuildAggregate(const int kTotalNumberOfCubes, const int kDimension, int* aggregate_position, double* center_of_mass) {
    
    for (int j=1;j<kTotalNumberOfCubes;j++) {
        
        // First, we compute the radius of the sphere onto which each random walker will be generated at the initial time t0 

        // Update center of mass of the aggregate as it gets built
        double* center_of_mass = new double[kDimension];
        // The first cube is centered at the origin [0,0,0]
        for (int d=0;d<kDimension;d++) *(center_of_mass+d) = 0.;

        for (int i=0;i<j;i++) {
            for (int d=0;d<kDimension;d++) *(center_of_mass+d) += *(aggregate_position+i*kDimension+d);
        }
        
        for (int d=0;d<kDimension;d++) *(center_of_mass+d) /= j;

        //agg.ComputeCenterOfMass(kDimension,j,aggregate_position,center_of_mass);
        // Find maximum distance b/w current random walker and the aggregate
        double radius = 0;
        double rho; // This is the radius of the sphere on which walkers will be generated
        for (int i=0;i<j;i++) {
            double tmp = 0.;
            for (int d=0;d<kDimension;d++) tmp += (*(aggregate_position+i*kDimension+d) - *(center_of_mass+d)) * (*(aggregate_position+i*kDimension+d) - *(center_of_mass+d));

            tmp = sqrt(tmp);

            radius = (tmp>radius) ? tmp : radius;
        }
        rho = radius + DR1;

        // Compute the final position of a new cube

        int* position_random_walker = new int[kDimension]; // New walker position on the grid. It has to be an integer

        bool no_contact = true;
        while (no_contact) {

            // Generate random points on surface of sphere uniformly
            // See here https://mathworld.wolfram.com/SpherePointPicking.html
            double mu = rand()/(RAND_MAX+1.); // taking out 2PI as cos(0)=cos(2pi) and sin(0)=sin(2pi)
            double nu = rand()/(RAND_MAX+1.); // in order to be able to get phi in [0,PI]
            // Follow typical math convention for angles: theta is the azimuthal angle, while phi is the polar angle
            double theta = 2*PI*mu;
            double phi = acos(2*nu-1);

            // Initial position of the walker is chosen randomly on the surface of a sphere centered at center of mass and of radius rho 
            double* position_random_walker_at_time_t0 = new double[kDimension];
            
            // Conversion from cartesian to spherical coordinates
            *(position_random_walker_at_time_t0+0) = *(center_of_mass+0) + rho*cos(theta)*sin(phi);
            *(position_random_walker_at_time_t0+1) = *(center_of_mass+1) + rho*sin(theta)*sin(phi);
            *(position_random_walker_at_time_t0+2) = *(center_of_mass+2) + rho*cos(phi);
            
            // Placing the random walker on a grid, since it takes only integer steps
            for (int d=0;d<kDimension;d++) *(position_random_walker+d) = DX*round(*(position_random_walker_at_time_t0+d)/DX);

            // Deallocate memory of initial walker position
            delete[] position_random_walker_at_time_t0;

            // Perform random walk
            bool keep_walking = true;
            while (keep_walking) {

                // Generate uniform random number between 0 and 1
                // 'ra' is a typical naming convention for it
                double ra = rand()/(RAND_MAX+1.);
                
                // Walker takes a step in one of the six direction based on the value of random number ra
                if      ( ra<1./6. ) *(position_random_walker+0) -= DX;
                else if ( ra<1./3. ) *(position_random_walker+1) -= DY;
                else if ( ra<0.5   ) *(position_random_walker+2) -= DZ;
                else if ( ra<2./3. ) *(position_random_walker+0) += DX;
                else if ( ra<5./6. ) *(position_random_walker+1) += DY;
                else                 *(position_random_walker+2) += DZ;

                // After the walker has taken a step, check distances b/w it and center of mass of the aggregate
                // and b/w it and each cube in the aggregate

                // Find distance b/w current random walker and center of mass of the aggregate
                // to check how far they are one from the other
                double distance = 0.;
                for (int d=0;d<kDimension;d++) distance += (*(position_random_walker+d) - *(center_of_mass+d)) * (*(position_random_walker+d) - *(center_of_mass+d));

                distance = sqrt(distance);

                // Check if walker has gone too far from the aggregate
                // and if so, dismiss walker
                if (distance >= rho+DR2) {
                    no_contact = true;    // Walker has not attached to the aggregate
                    keep_walking = false; // but it has gone too far away from the aggregate and we dismiss it. New walker gets introduced.

                    break;
                }

                // If the walker has not gone too far away then
                // find distance between each cube in the aggregate and the current random walker
                // to determine the shortest distance between each member of the aggregate and the walker
                double shortest_distance = 2.*rho+DR2;    // Starting with a big number

                for (int i=0;i<j;i++) {
                    double tmp = 0.;

                    for (int d=0; d<kDimension;d++) tmp += (*(position_random_walker+d) - *(aggregate_position+i*kDimension+d)) * (*(position_random_walker+d) - *(aggregate_position+i*kDimension+d));

                    tmp = sqrt(tmp);

                    shortest_distance =  (tmp<shortest_distance) ? tmp : shortest_distance;

                    // Check if the shortest distance b/w aggregate and random walker is in attaching range
                    // and if so attach walker to aggregate
                    if (shortest_distance <= ATTACH) {
                        no_contact = false;   // Walker has attached to the aggregate
                        keep_walking = false; // and thus it stops walking. New walker gets initialized
                    }

                    break;
                }
            }
        }

        for (int d=0;d<kDimension;d++) *(aggregate_position+j*kDimension+d) = *(position_random_walker+d); // Assign walker's final position to new member of the aggregate

        // Deallocate memory of arrays of center of mass and current random walker position
        delete[] center_of_mass;
        delete[] position_random_walker;
    }
}

// This function prints the position of each cube in the final aggregate 
void PrintAggregatePosition(FILE* fp, const int kTotalNumberOfCubes, const int kDimension, int* aggregate_position) {
    
    for (int j=0;j<kTotalNumberOfCubes;j++) {
        for (int d=0;d<kDimension;d++) fprintf(fp,"%d\t",*(aggregate_position+j*kDimension+d));
        fprintf(fp,"\n");
    }
}
