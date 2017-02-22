
#include "ExplicitUpwindParallel.h"


ExplicitUpwindParallel::ExplicitUpwindParallel(double xMin,
                                               double xMax,
                                               double time,
                                               double numberOfSpacePoints,
                                               double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                                          numberOfSpacePoints, CFL),
                                                             methodName("ExplicitUpwindParallelScheme")
{

}

ExplicitUpwindParallel::~ExplicitUpwindParallel()
{
}


void ExplicitUpwindParallel::solve(int setNumber)
{
    try
    {

        std::cout << methodName << " solution runs and matrix is initialised\n";
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);

        double timeOfStart = MPI_Wtime();
        explicitParallelResults = Matrix(numberOfSpacePoints, numberOfTimePoints);
        lastProc = numOfProc - 1;
        if (numberOfSpacePoints % numOfProc != 0)
        {
            throw "error";
        }

        localLimit = numberOfSpacePoints / numOfProc;


        int limitLow = myRank * localLimit;
        int limitHigh = (myRank + 1) * localLimit;

        double actualValue = xMin;


        for (int i = limitLow; i < limitHigh; ++i)
        {

            actualValue = (i * (*this).dx) + xMin;
            explicitParallelResults[i][0] = (1.0 / 2.0) * (*this).initializationFunction(2, actualValue);

        }

        if (myRank == 0)
        {
            for (int i = 0; i < numberOfTimePoints; ++i)
            {
                explicitParallelResults[0][i] = 0;
            }
        }

        if (myRank == lastProc)
        {

            for (int i = 0; i < numberOfTimePoints; ++i)
            {

                explicitParallelResults[numberOfSpacePoints - 1][i] = 0;
            }

        }


        //explicitResutls = Matrix((*this).getMatrix());


        for (int j = 0; j < numberOfTimePoints - 1; ++j)
        {

            if (myRank != 0)
            {
                double localTmp;
                MPI_Recv(&localTmp, 1, MPI_DOUBLE, myRank - 1, 1, MPI_COMM_WORLD, &status);
                explicitParallelResults[limitLow][j + 1] = (explicitParallelResults[limitLow][j] -
                                                            CFL * (explicitParallelResults[limitLow][j] - localTmp));
            }


            for (int i = limitLow + 1; i < limitHigh; ++i)
            {

                explicitParallelResults[i][j + 1] = (explicitParallelResults[i][j] -
                                                     CFL * (explicitParallelResults[i][j] -
                                                            explicitParallelResults[i - 1][j]));
            }

            if (lastProc != myRank)
            {
                MPI_Send(&explicitParallelResults[limitHigh - 1][j], 1, MPI_DOUBLE, myRank + 1, 1, MPI_COMM_WORLD);

            }

        }


        gatherResults.resize(numberOfSpacePoints);

        std::vector<double> tempForReduce(numberOfSpacePoints);

        for (int i = 0; i < numberOfSpacePoints; ++i)
        {
            tempForReduce[i] = explicitParallelResults[i][numberOfTimePoints - 1];
        }


        MPI_Reduce(&tempForReduce[0], &gatherResults[0], numberOfSpacePoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



        //GeneralScheme::solve(setNumber);
        //calculateNorms((*this).explicitResutls);
        double timeOfEnd = MPI_Wtime();
        double totalTime = timeOfEnd - timeOfStart;
        std::cout << methodName << " time is: " << totalTime << std::endl;


    }

    catch (std::exception &e)
    {
        std::cout << "Standard exception: " << e.what() << std::endl;
        return;
    }
}


/*Matrix ExplicitUpwindScheme::getUpwindMatrix() {
    return explicitResutls;
}*/

std::string ExplicitUpwindParallel::getName()
{
    return methodName;
}

std::vector<double> ExplicitUpwindParallel::getLastExplicitParallelMatrixColumn()
{
    return explicitParallelResults.getColumn(numberOfTimePoints - 1);
}

const std::vector<double> &ExplicitUpwindParallel::getGatherResults() const
{
    return gatherResults;
}
