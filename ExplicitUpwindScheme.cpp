#include "ExplicitUpwindScheme.h"


ExplicitUpwindScheme::ExplicitUpwindScheme(double xMin,
                                           double xMax,
                                           double time,
                                           double numberOfSpacePoints,
                                           double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                                      numberOfSpacePoints, CFL),
                                                         methodName("ExplicitUpwindScheme") {

}


ExplicitUpwindScheme::~ExplicitUpwindScheme() {
}


void ExplicitUpwindScheme::solve(int setNumber) {
    try {

        std::cout << "Explicit upwid scheme solution runs and matrix is initialised\n";

        (*this).initializeSet(setNumber);
        explicitResutls = Matrix((*this).getMatrix());

        for (int j = 0; j < numberOfTimePoints - 1; ++j) {
            for (int i = 1; i < numberOfSpacePoints; ++i) {
                explicitResutls[i][j + 1] = (explicitResutls[i][j] -
                                             CFL * (explicitResutls[i][j] - explicitResutls[i - 1][j]));
            }

        }
        GeneralScheme::solve(setNumber);
        calculateNorms((*this).explicitResutls);

    }

    catch (std::exception &e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }
}

Matrix ExplicitUpwindScheme::getUpwindMatrix() {
    return explicitResutls;
}

std::string ExplicitUpwindScheme::getName() {
    return methodName;
}


