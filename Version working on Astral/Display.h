#pragma once

#include <vector>
#include <iostream>


/**
@brief Own namespace to display vectors and names
*/

namespace display
{

    /**
    @brief Fucntion to display given input vector to console

    @param Vector to display
    */
    void displayVector(std::vector<double> vector);

    /**
    @brief Display name of setected boundary condition. This function contains two sets of boundary condition: exponential and sign
    @param Number of Boundary condition to dsplay. Takes values 1 or 2
    */
    std::string getInitialBoundaryConditionName(int &numberOfBoundaryCondition);
}
