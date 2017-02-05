#include <mpi/mpi.h>
#include "ExplicitUpwindParallel.h"


ExplicitUpwindParallel::ExplicitUpwindParallel(double xMin,
                                               double xMax,
                                               double time,
                                               double numberOfSpacePoints,
                                               double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                                          numberOfSpacePoints, CFL),
                                                             methodName("ExplicitUpwindParallelScheme") {

}

ExplicitUpwindParallel::~ExplicitUpwindParallel() {
}


void ExplicitUpwindParallel::solve(int setNumber) {
    try {

        std::cout << methodName << " solution runs and matrix is initialised\n";
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

        double timeOfStart = MPI_Wtime();
        explicitResutls = Matrix(numberOfSpacePoints, numberOfTimePoints);
        lastNode = worldSize - 1;
        if (numberOfSpacePoints % worldSize != 0) {
            throw "error";
        }

        localLimit = numberOfSpacePoints / worldSize;


        int limitLow = myRank * localLimit;
        int limitHigh = (myRank + 1) * localLimit;

        double actualValue = xMin;


        for (int i = limitLow; i < limitHigh; ++i) {

            actualValue = (i * (*this).dx) + xMin;
            explicitResutls[i][0] = (1.0 / 2.0) * (*this).initializationFunction(1, actualValue);

        }

        if (myRank == 0) {
            for (int i = 0; i < numberOfTimePoints; ++i) {
                explicitResutls[0][i] = 0;
            }
        }

        if (myRank == lastNode) {

            for (int i = 0; i < numberOfTimePoints; ++i) {

                explicitResutls[numberOfSpacePoints - 1][i] = 0;
            }

        }


        //explicitResutls = Matrix((*this).getMatrix());


        for (int j = 0; j < numberOfTimePoints - 1; ++j) {

            if (myRank != 0) {
                double localTmp;
                MPI_Recv(&localTmp, 1, MPI_DOUBLE, myRank - 1, 1, MPI_COMM_WORLD, &status);
                explicitResutls[limitLow][j + 1] = (explicitResutls[limitLow][j] -
                                                    CFL * (explicitResutls[limitLow][j] - localTmp));
            }


            for (int i = limitLow + 1; i < limitHigh; ++i) {

                explicitResutls[i][j + 1] = (explicitResutls[i][j] -
                                             CFL * (explicitResutls[i][j] - explicitResutls[i - 1][j]));
            }

            if (lastNode != myRank) {
                MPI_Send(&explicitResutls[limitHigh - 1][j], 1, MPI_DOUBLE, myRank + 1, 1, MPI_COMM_WORLD);

            }

        }


        gatherResults.resize(numberOfSpacePoints);

        std::vector<double> tempForReduce(numberOfSpacePoints);

        for (int i = 0; i < numberOfSpacePoints; ++i) {
            tempForReduce[i] = explicitResutls[i][numberOfTimePoints - 1];
        }


        MPI_Reduce(&tempForReduce[0], &gatherResults[0], numberOfSpacePoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



        //GeneralScheme::solve(setNumber);
        //calculateNorms((*this).explicitResutls);
        double timeOfEnd = MPI_Wtime();
        double totalTime = timeOfEnd - timeOfStart;
        std::cout << "Explicit upwind parallel time is: " << totalTime << std::endl;



    }

    catch (std::exception &e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
        return;
    }
}


/*Matrix ExplicitUpwindScheme::getUpwindMatrix() {
    return explicitResutls;
}*/

std::string ExplicitUpwindParallel::getName() {
    return methodName;
}

std::vector<double> ExplicitUpwindParallel::getLastExplicitParallelMatrixColumn() {
    return explicitResutls.getColumn(numberOfTimePoints - 1);
}

const std::vector<double> &ExplicitUpwindParallel::getGatherResults() const {
    return gatherResults;
}





