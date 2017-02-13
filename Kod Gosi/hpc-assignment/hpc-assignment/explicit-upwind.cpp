#include <mpi.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>

double solveAnalytical(double tmp, double x, double t);

unsigned int divide(int myRank, int numOfNodes, unsigned int taskComplexity);

void divide(double *vector, double Boundary_L, double acceleration, unsigned int selectedChunk, double dx);


std::string getCurrentPath()
{
    char charCurrentPath[1024];
    std::string path;
    if (getcwd(charCurrentPath, sizeof(charCurrentPath)) != NULL)
    {
        path = std::string(charCurrentPath);
    } else
        perror("getcwd() error");


    return path;

}


int main(int argc, char *argv[])
{
    int myRank;
    int numOfNodes;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfNodes);

    double minimumGridValue = -50;
    double maksimumGridValue = 50;
    double acceleration = 1.75;
    unsigned int numberOfPoints = 1000;
    double upperTimeBondary = 5.0;
    double CourantNumber = 0.9;
    double dx = (maksimumGridValue - minimumGridValue) / numberOfPoints;
    double dt = (CourantNumber * dx) / acceleration;

    double startTime = MPI_Wtime();

    double valueOfBoundaryLeft = 0;
    unsigned int selectedChunk = divide(myRank, numOfNodes, numberOfPoints);
    double actualTime = 0.0;
    double *previousValues = new double[selectedChunk];
    double *actualValues = new double[selectedChunk];
    divide(previousValues, minimumGridValue + myRank * dx * floor(numberOfPoints / numOfNodes), acceleration,
           selectedChunk, dx);

    while (upperTimeBondary > actualTime)
    {
        if (myRank == 0)
        {
            valueOfBoundaryLeft = 0;
        } else
        {
            MPI_Status myActualStatus;
            double tmp;
            MPI_Recv(&tmp, 1, MPI_DOUBLE, myRank - 1, 0, MPI_COMM_WORLD, &myActualStatus);
            valueOfBoundaryLeft = tmp;
        }

        if (myRank != numOfNodes - 1)
        {
            MPI_Send(&previousValues[selectedChunk - 1], 1, MPI_DOUBLE, myRank + 1, 0, MPI_COMM_WORLD);
        }

        actualValues[0] = previousValues[0] - CourantNumber * (previousValues[0] - valueOfBoundaryLeft);
        for (unsigned int pointInSpace = 1; pointInSpace < selectedChunk; pointInSpace++)
        {
            actualValues[pointInSpace] = previousValues[pointInSpace] - CourantNumber * (previousValues[pointInSpace] -
                                                                                         previousValues[pointInSpace -
                                                                                                        1]);
        }

        for (unsigned int i = 0; i < selectedChunk; i++)
        {
            previousValues[i] = actualValues[i];
        }

        actualTime += dt;
    }


    if (myRank == 0)
    {
        double *explicitOutput = new double[numberOfPoints];
        double *it = explicitOutput;

        for (unsigned int i = 0; i < selectedChunk; i++)
        {
            explicitOutput[i] = previousValues[i];
        }

        for (unsigned int i = 1; i < numOfNodes; i++)
        {
            MPI_Status myActualStatus;
            unsigned int taskComplexity = divide(i - 1, numOfNodes, numberOfPoints);
            it += taskComplexity;
            MPI_Recv(it, divide(i, numOfNodes, numberOfPoints), MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myActualStatus);
        }

        std::fstream storeResults;
        std::string tmpPath = getCurrentPath() + "/x";
        const char *currentPath = tmpPath.c_str();
        storeResults.open(currentPath, std::fstream::trunc | std::fstream::out);
        for (unsigned int i = 0; i < numberOfPoints; i++)
        {
            storeResults << solveAnalytical(acceleration, minimumGridValue + i * dx, actualTime - dt)
                         << ","
                         << explicitOutput[i]
                         << std::endl;
        }

        double CalculateNorm = 0.0;
        for (unsigned int i = 0; i < numberOfPoints; i++)
        {
            CalculateNorm += fabs(
                    solveAnalytical(acceleration, minimumGridValue + i * dx, actualTime - dt) - explicitOutput[i]);
        }

        std::cout << "One CalculateNorm is: " << CalculateNorm << std::endl;
        storeResults << "One CalculateNorm is: " << CalculateNorm << std::endl;

        std::cout << "Whole program execution time:" << MPI_Wtime() - startTime << std::endl;
        storeResults << "Whole program execution time: " << MPI_Wtime() - startTime << std::endl;

        storeResults.close();

        delete[] explicitOutput;
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

void divide(double *vector, double Boundary_L, double acceleration, unsigned int selectedChunk, double dx)
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

unsigned int divide(int myRank, int numOfNodes, unsigned int taskComplexity)
{
    if (numOfNodes == 1)
        return taskComplexity;
    else if (taskComplexity % numOfNodes == 0)
        return taskComplexity / numOfNodes;
    else
    {
        if (myRank == numOfNodes - 1)
            return taskComplexity / numOfNodes + taskComplexity % numOfNodes;
        else
            return taskComplexity / numOfNodes;
    }
}

