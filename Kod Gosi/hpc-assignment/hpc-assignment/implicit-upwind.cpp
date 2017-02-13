#include <mpi.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>


unsigned int divide(int myRank, int numOfCores, unsigned int taskComplexity);

void initOfVector(double *vector, double Boundary_L, double acceleration, unsigned int selectedChunk, double dx);

double solveAnalytical(double a, double x, double t);

int main(int argc, char *argv[])
{
    int myRank;
    int numOfCores;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfCores);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


    double minimumGridValue = -50;
    double maksimumGridValue = 50;
    double acceleration = 1.75;
    unsigned int numberOfPoints = 1000;
    double upperTimeBondary = 5.0;
    double CourantNumber = 0.999;
    double dx = (maksimumGridValue - minimumGridValue) / numberOfPoints;
    double dt = (CourantNumber * dx) / acceleration;

    double startTime = MPI_Wtime();

    double valueOfBoundaryLeft = 0;
    unsigned int selectedChunk = divide(myRank, numOfCores, numberOfPoints);
    double actualTime = 0.0;
    double *previousValues = new double[selectedChunk];
    double *actualValues = new double[selectedChunk];
    initOfVector(previousValues, minimumGridValue + myRank * dx * floor(numberOfPoints / numOfCores), acceleration,
                 selectedChunk, dx);

    while (upperTimeBondary > actualTime)
    {
        if (myRank == 0)
        {
            valueOfBoundaryLeft = 0;
        } else
        {
            MPI_Status myActualStatus;
            double a;
            MPI_Recv(&a, 1, MPI_DOUBLE, myRank - 1, 0, MPI_COMM_WORLD, &myActualStatus);
            valueOfBoundaryLeft = a;
        }

        actualValues[0] = (previousValues[0] + CourantNumber * valueOfBoundaryLeft) / (1.0 + CourantNumber);
        for (unsigned int pointInSpace = 1; pointInSpace < selectedChunk; pointInSpace++)
        {
            actualValues[pointInSpace] =
                    (previousValues[pointInSpace] + CourantNumber * actualValues[pointInSpace - 1]) /
                    (1.0 + CourantNumber);
        }

        if (myRank != numOfCores - 1)
        {
            MPI_Send(&actualValues[selectedChunk - 1], 1, MPI_DOUBLE, myRank + 1, 0, MPI_COMM_WORLD);
        }

        for (unsigned int i = 0; i < selectedChunk; i++)
        {
            previousValues[i] = actualValues[i];
        }

        actualTime += dt;
    }


    if (myRank == 0)
    {
        double *implicitOutput = new double[numberOfPoints];
        double *it = implicitOutput;

        for (unsigned int i = 0; i < selectedChunk; i++)
        {
            implicitOutput[i] = previousValues[i];
        }

        for (unsigned int i = 1; i < numOfCores; i++)
        {
            MPI_Status myActualStatus;
            unsigned int taskComplexity = divide(i - 1, numOfCores, numberOfPoints);
            it += taskComplexity;
            MPI_Recv(it, divide(i, numOfCores, numberOfPoints), MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myActualStatus);
        }

        std::fstream storeResults;
        storeResults.open("implicitScheme.csv", std::fstream::trunc | std::fstream::out);
        for (unsigned int i = 0; i < numberOfPoints; i++)
        {
            storeResults << solveAnalytical(acceleration, minimumGridValue + i * dx, actualTime - dt)
                         << ","
                         << implicitOutput[i]
                         << std::endl;
        }

        double CalculateNorm = 0.0;
        for (unsigned int i = 0; i < numberOfPoints; i++)
        {
            CalculateNorm += fabs(
                    solveAnalytical(acceleration, minimumGridValue + i * dx, actualTime - dt) - implicitOutput[i]);
        }

        std::cout << "One CalculateNorm is: " << CalculateNorm << std::endl;
        storeResults << "One CalculateNorm is: " << CalculateNorm << std::endl;

        std::cout << "Whole program execution time:" << MPI_Wtime() - startTime << std::endl;
        storeResults << "Whole program execution time:" << MPI_Wtime() - startTime << std::endl;

        storeResults.close();

        delete[] implicitOutput;
    } else
    {

        MPI_Send(previousValues, selectedChunk, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    delete[] previousValues;
    delete[] actualValues;
    MPI_Finalize();
    return 0;
}

void initOfVector(double *vector, double Boundary_L, double acceleration, unsigned int selectedChunk, double dx)
{
    for (unsigned int i = 0; i < selectedChunk; i++)
    {
        vector[i] = solveAnalytical(acceleration, Boundary_L + i * dx, 0.0);
    }
}

double solveAnalytical(double a, double x, double t)
{
    return 0.5 * exp(-(x - a * t) * (x - a * t));
}

unsigned int divide(int myRank, int numOfCores, unsigned int taskComplexity)
{
    if (numOfCores == 1)
        return taskComplexity;
    else if (taskComplexity % numOfCores == 0)
        return taskComplexity / numOfCores;
    else
    {
        if (myRank == numOfCores - 1)
            return taskComplexity / numOfCores + taskComplexity % numOfCores;
        else
            return taskComplexity / numOfCores;
    }
}

