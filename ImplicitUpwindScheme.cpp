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


        double timeOfStart = MPI_Wtime();
        std::cout << "Explicit upwind scheme solution runs and matrix is initialised\n";

        implicitResults = Matrix(numberOfSpacePoints, numberOfTimePoints);
        double actualValue = xMin;

        for (int i = 0; i < numberOfSpacePoints; ++i)
        {
            actualValue = (i * (*this).dx) + xMin;
            implicitResults[i][0] = (1.0 / 2.0) * (*this).initializationFunction(setNumber, actualValue);

        }


        for (int i = 0; i < numberOfTimePoints; ++i)
        {
            implicitResults[0][i] = 0;
            implicitResults[numberOfSpacePoints - 1][i] = 2;

        }

        //(*this).initializeSet(setNumber);
        //implicitResults = Matrix((*this).getMatrix());


        for (int j = 0; j < numberOfTimePoints - 1; ++j) {
            for (int i = 1; i < numberOfSpacePoints; ++i) {

                implicitResults[i][j + 1] = (-implicitResults[i][j] -
                                             CFL * implicitResults[i - 1][j]) / -(1 + CFL);

            }

        }

        double timeOfEnd = MPI_Wtime();
        double totalTime = timeOfEnd - timeOfStart;
        std::cout << methodName << " time is: " << totalTime << std::endl;
        //GeneralScheme::solve(setNumber);
        //calculateNorms((*this).implicitResults);

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

std::vector<double> ImplicitUpwindScheme::getResults()
{
    return implicitResults.getColumn(numberOfTimePoints - 1);
}

