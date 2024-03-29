
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

        double timeOfStart = MPI_Wtime();
        std::cout << "Explicit upwind scheme solution runs and matrix is initialised\n";

        explicitResutls = Matrix(numberOfSpacePoints, numberOfTimePoints);
        double actualValue = xMin;

        for (int i = 0; i < numberOfSpacePoints; ++i)
        {
            actualValue = (i * (*this).dx) + xMin;
            explicitResutls[i][0] = (1.0 / 2.0) * (*this).initializationFunction(setNumber, actualValue);

        }


        for (int i = 0; i < numberOfTimePoints; ++i)
        {
            explicitResutls[0][i] = 0;
            explicitResutls[numberOfSpacePoints - 1][i] = 2;

        }

//        (*this).isSetInitialised = true;
//        (*this).initializeSet(setNumber);
//        explicitResutls = Matrix((*this).getMatrix());

        for (int j = 0; j < numberOfTimePoints - 1; ++j) {
            for (int i = 1; i < numberOfSpacePoints; ++i) {
                explicitResutls[i][j + 1] = (explicitResutls[i][j] -
                                             CFL * (explicitResutls[i][j] - explicitResutls[i - 1][j]));
            }

        }
        double timeOfSolutionsEnd = MPI_Wtime();
        double solutionTime = timeOfSolutionsEnd - timeOfStart;
        std::cout << "Single Explicit upwind scheme solution is: " << solutionTime << std::endl;

        //GeneralScheme::solve(setNumber);
        //calculateNorms((*this).explicitResutls);

    }

    catch (std::exception &e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }
}

/*Matrix ExplicitUpwindScheme::getUpwindMatrix() {
    return explicitResutls;
}*/

std::string ExplicitUpwindScheme::getName() {
    return methodName;
}

std::vector<double> ExplicitUpwindScheme::getLastExplicitMatrixColumn() {
    return explicitResutls.getColumn(numberOfTimePoints - 1);
}


