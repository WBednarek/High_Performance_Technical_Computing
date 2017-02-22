
#include "ImplicitParallel.h"


ImplicitParallel::ImplicitParallel(double xMin,
                                   double xMax,
                                   double time,
                                   double numberOfSpacePoints,
                                   double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                              numberOfSpacePoints, CFL),
                                                 methodName("ImplicitUpwindParallelScheme")
{

}

ImplicitParallel::~ImplicitParallel()
{
}


void ImplicitParallel::solve(int setNumber)
{
    try
    {

        std::cout << methodName << " solution runs and matrix is initialised\n";
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
        implicitResultsParallel = Matrix(numberOfSpacePoints, numberOfTimePoints);
        double timeOfStart = MPI_Wtime();

        lastNode = numOfProc - 1;
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
            implicitResultsParallel[i][0] = (1.0 / 2.0) * (*this).initializationFunction(setNumber, actualValue);

        }

        if (myRank == 0)
        {
            for (int i = 0; i < numberOfTimePoints; ++i)
            {
                implicitResultsParallel[0][i] = 0;
            }
        }

        if (myRank == lastNode)
        {

            for (int i = 0; i < numberOfTimePoints; ++i)
            {

                implicitResultsParallel[numberOfSpacePoints - 1][i] = 0;
            }

        }


        //explicitResutls = Matrix((*this).getMatrix());


        for (int j = 0; j < numberOfTimePoints - 1; ++j)
        {

            if (myRank != 0)
            {
                double localTmp;
                MPI_Recv(&localTmp, 1, MPI_DOUBLE, myRank - 1, 1, MPI_COMM_WORLD, &status);
                implicitResultsParallel[limitLow][j + 1] = (-implicitResultsParallel[limitLow][j] -
                                                            CFL * localTmp) / -(1 + CFL);
            }


            for (int i = limitLow + 1; i < limitHigh; ++i)
            {

                implicitResultsParallel[i][j + 1] = (-implicitResultsParallel[i][j] -
                                                     CFL * implicitResultsParallel[i - 1][j]) / -(1 + CFL);
            }

            if (lastNode != myRank)
            {
                MPI_Send(&implicitResultsParallel[limitHigh - 1][j], 1, MPI_DOUBLE, myRank + 1, 1, MPI_COMM_WORLD);

            }

        }


        gatherResults.resize(numberOfSpacePoints);

        std::vector<double> tempForReduce(numberOfSpacePoints);

        for (int i = 0; i < numberOfSpacePoints; ++i)
        {
            tempForReduce[i] = implicitResultsParallel[i][numberOfTimePoints - 1];
        }


        MPI_Reduce(&tempForReduce[0], &gatherResults[0], numberOfSpacePoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



        //GeneralScheme::solve(setNumber);
        //calculateNorms((*this).implicitResultsParallel);
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

std::string ImplicitParallel::getName()
{
    return methodName;
}

std::vector<double> ImplicitParallel::getLastExplicitParallelMatrixColumn()
{
    return implicitResultsParallel.getColumn(numberOfTimePoints - 1);
}

const std::vector<double> &ImplicitParallel::getImplicitParallelResults() const
{
    return gatherResults;
}





