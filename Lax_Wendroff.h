#pragma once

#include "GeneralScheme.h"

class Lax_Wendroff :
        public GeneralScheme {

    std::string methodName;
    Matrix Lax_WendroffResutls;

public:
    Lax_Wendroff();

    Lax_Wendroff(double xMin,
                 double xMax,
                 double time,
                 double numberOfSpacePoints,
                 double CFL);

    virtual ~Lax_Wendroff();

    /**
    @brief Virtual method which solves Lax-Wendroff

    @param Variable for inidicating boundary initialization set type, 1 or 2 respectively for:  sign type and exponential type

    */
    virtual void solve(int setNumber);
    //double solutionFunctionExplicitScheme(int numberOfSet, Matrix toUpwindSchemeCalculations);
    /**
    @brief Method returns Matrix where Lax-Wendroff solution is stored

    @return Matrix with calculated Lax-Wendroff method

    */
    Matrix getLax_WendroffdMatrix();

    /**
    @brief Virtual method to returns name of Lax-Wendroff class

    @return name of Lax-Wendroff class
    */
    virtual std::string getName();

};

