#pragma once

#include "GeneralScheme.h"


/**
* Class created for calculating Explicit Upwind Scheme
* Class inherits methods from GeneralScheme class

*/

class ImplicitParallel : public GeneralScheme
{

    std::string methodName;
    Matrix implicitResults;
    int myRank;
    int numOfProc;
    double lastNode;
    double localLimit;
    std::vector<double> gatherResults;


private:
    MPI_Status status;


public:

    ImplicitParallel(double xMin,
                     double xMax,
                     double time,
                     double numberOfSpacePoints,
                     double CFL);

    ~ImplicitParallel();

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

    virtual const std::vector<double> &getResults() const;


};



