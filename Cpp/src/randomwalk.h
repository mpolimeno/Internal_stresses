#ifndef RANDOMWALK_H_INCLUDED
#define RANDOMWALK_H_INCLUDED

// define the value of the number PI
#define PI       (4.*atan(1.))

// parameters for the random walk
#define DX      2           // stepsize in x // at this stage this code is set up only for DX,DY,DZ positive integers and they have to be equal
#define DY      2           // stepsize in y == stepsize in x
#define DZ      2           // stepsize in z == stepsize in x == stepsize in y
#define DR1     10.         // delta radius to increase/decrease size of the sphere on whose surface each walker is generated
#define DR2     10.         // how far away from the sphere we allow each new walker to be w/o being dismissed from the aggregation process
#define CASE    1           // Types of attachments allowed: if CASE=1 -> face-to-face; if CASE=2-> edge-to-edge; if CASE=3 -> corner-to-corner
#define EPS     1.e-4       // small number to avoid issues with the face-to-face attachment

// currently this code assumes DX=DY=DZ
#if CASE==1
#define ATTACH  (2.+EPS)            // cubes can attach only face to face to any member of the aggregate
#elif CASE==2
#define ATTACH  (DX*sqrt(2.))       // cubes can attach at least edge to edge to any member of the aggregate
#elif CASE==3
#define ATTACH  (DX*sqrt(3.))       // cubes can attach at least corner to corner to any member of the aggregate
#endif



extern void BuildAggregate(const int kTotalNumberOfCubes, const int kDimension, int* aggregate_postion);
//extern void BuildAggregate(const int kTotalNumberOfCubes, const int kDimension, int* aggregate_postion, double* center_of_mass);
extern void PrintAggregatePosition(FILE* fp, const int kTotalNumberOfCubes, const int kDimension, int* aggregate_position);


#endif // RANDOMWALK_H_INCLUDED
