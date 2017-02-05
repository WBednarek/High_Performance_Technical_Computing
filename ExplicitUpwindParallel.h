#pragma once

#include "GeneralScheme.h"


/**
* Class created for calculating Explicit Upwind Scheme
* Class inherits methods from GeneralScheme class
*/

class ExplicitUpwindParallel : public GeneralScheme {

    std::string methodName;
    Matrix explicitParallelResults;
    int myRank;
    int numOfProc;
    double lastProc;
    double workAdditional;
    double tmp;
    double localLimit;
    std::vector<double> gatherResults;
public:
    const std::vector<double> &getGatherResults() const;


private:
    MPI_Status status;


public:

    ExplicitUpwindParallel(double xMin,
                           double xMax,
                           double time,
                           double numberOfSpacePoints,
                           double CFL);

    ~ExplicitUpwindParallel();

    /**
    @brief Virtual method which solves Explicit Upwind Scheme
    @param Variable for inidicating boundary initialization set type, 1 or 2 respectively for:  sign type and exponential type
    */
    virtual void solve(int setNumber);

    /**
    @brief Method returns Matrix where Explicit Upwind Scheme solution is stored
    @return Matrix with calculated  Explicit Upwind Scheme
    */
    /* Matrix getUpwindMatrix();*/

    /**
    @brief Virtual method to returns name of ExplicitUpwindScheme class
    @return name of ExplicitUpwindScheme class
    */
    virtual std::string getName();

    std::vector<double> getLastExplicitParallelMatrixColumn();




};


