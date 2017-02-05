#include "ImplicitUpwindScheme.h"


ImplicitUpwindScheme::ImplicitUpwindScheme() {
}

ImplicitUpwindScheme::ImplicitUpwindScheme(double xMin,
                                           double xMax,
                                           double time,
                                           double numberOfSpacePoints,
                                           double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                                      numberOfSpacePoints, CFL),
                                                         methodName("ImplicitUpwindScheme") {
}


ImplicitUpwindScheme::~ImplicitUpwindScheme() {
}

void ImplicitUpwindScheme::solve(int setNumber) {
    try {

        std::cout << methodName << " scheme solution runs and matrix is initialised\n";


        (*this).initializeSet(setNumber);
        implicitResults = Matrix((*this).getMatrix());


        for (int j = 0; j < numberOfTimePoints - 1; ++j) {
            for (int i = 1; i < numberOfSpacePoints; ++i) {

                implicitResults[i][j + 1] = (-1.0 * CFL) * (implicitResults[i][j + 1] - implicitResults[i - 1][j + 1]) +
                                            implicitResults[i][j];

            }

        }
        GeneralScheme::solve(setNumber);
        calculateNorms((*this).implicitResults);

    }

    catch (std::exception &e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }
}

Matrix ImplicitUpwindScheme::getImplicitUpwindMatrix() {
    return implicitResults;
}

std::string ImplicitUpwindScheme::getName() {
    return methodName;
}

std::vector<double> ImplicitUpwindScheme::getLastImplicitMatrixColumn() {
    return implicitResults.getColumn(numberOfTimePoints - 1);
}

