#pragma once

#include "GeneralScheme.h"

class ImplicitUpwindScheme :
        public GeneralScheme {

    std::string methodName;
    Matrix implicitResults;
public:
    ImplicitUpwindScheme();

    ImplicitUpwindScheme(double xMin,
                         double xMax,
                         double time,
                         double numberOfSpacePoints,
                         double CFL);

    virtual ~ImplicitUpwindScheme();

    /**
    @brief Virtual method which solves Implicit Upwind Scheme

    @param Variable for inidicating boundary initialization set type, 1 or 2 respectively for:  sign type and exponential type

    */
    virtual void solve(int setNumber);

    /**
    @brief Method returns Matrix where Implicit Upwind Scheme solution is stored

    @return Matrix with calculated Implicit Upwind Schem

    */
    Matrix getImplicitUpwindMatrix();

    /**
    @brief Virtual method to returns name of  Implicit Upwind Scheme class

    @return name of  Implicit Upwind Scheme class
    */
    virtual std::string getName();


    std::vector<double> getLastImplicitMatrixColumn();
};

