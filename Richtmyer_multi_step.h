#pragma once

#include "GeneralScheme.h"

class Richtmyer_multi_step :
        public GeneralScheme {
    std::string methodName;
    Matrix RichtmyerResutls;
    Matrix halfStepRichtmyer;

public:
    Richtmyer_multi_step();

    Richtmyer_multi_step(double xMin,
                         double xMax,
                         double time,
                         double numberOfSpacePoints,
                         double CFL);

    virtual ~Richtmyer_multi_step();

    /**
    @brief Virtual method which solves Richtmyer_multi_step method

    @param Variable for inidicating boundary initialization set type, 1 or 2 respectively for:  sign type and exponential type

    */

    //Virtual solve method as in the rest classes. Overide is an "tip" for compiler in cases then method could be by an accident not overrided but created as new one
    virtual void solve(int setNumber);

    /**
    @brief Method returns Matrix where Richtmyer_multi_step solution is stored

    @return Matrix with calculated Richtmyer_multi_step

    */
    Matrix getRichtmyer_multi_stepdMatrix();

    /**
    @brief Virtual method to returns name of Richtmyer_multi_ste class

    @return name of Richtmyer_multi_step class
    */
    virtual std::string getName();

};

