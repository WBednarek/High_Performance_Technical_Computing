#include <mpi/mpi.h>
#include "ExplicitUpwindParallel.h"


ExplicitUpwindParallel::ExplicitUpwindParallel(double xMin,
                                               double xMax,
                                               double time,
                                               double numberOfSpacePoints,
                                               double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                                          numberOfSpacePoints, CFL),
                                                             methodName("ExplicitUpwindScheme") {

}


ExplicitUpwindParallel::~ExplicitUpwindParallel() {
}


void ExplicitUpwindParallel::solve(int setNumber) {
    try {

        std::cout << "Explicit upwind scheme solution runs and matrix is initialised\n";
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Comm_size(MPI_COMM_WORLD, &npes);

        double timeOfStart;
        double timeOfEnd;


        (*this).initializeSet(setNumber);
        timeOfStart = MPI_Wtime();
        explicitResutls = Matrix((*this).getMatrix());
        timeOfEnd = MPI_Wtime();
        double programExecutionTime = timeOfEnd - timeOfStart;
        std::cout << "Total time of initializing the same Matris is: " << programExecutionTime << std::endl;

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

/*Matrix ExplicitUpwindScheme::getUpwindMatrix() {
    return explicitResutls;
}*/

std::string ExplicitUpwindParallel::getName() {
    return methodName;
}

std::vector<double> ExplicitUpwindParallel::getLastExplicitMatrixColumn() {
    return explicitResutls.getColumn(numberOfTimePoints - 1);
}


