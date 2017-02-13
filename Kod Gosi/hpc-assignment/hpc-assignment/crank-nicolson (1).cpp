#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <unistd.h>

#define TAG_ONE 0
#define TAG_TWO 1
#define TAG_THREE 2
#define TAG_FOUR 3
#define TAG_FIVE 4
#define TAG_SIX 5
#define TAG_SEVEN 6
#define TAG_EIGHT 7

using namespace std;

int getProcArraySize(int myNode, int totalNodes, int totalNumberOfPoints);

double getInitialValue(double totalNumberOfPointsX, double deltaX, int j);

double calculateAnalyticalSolution(double totalNumberOfPointsX, double timeValue);

void solve(int numOfPoints, double deltaX, double CFL, double *initData, double *finalWave, int myNode, int totalNodes,
           double *l, double *u, double *d, double *y);

void LUPreparation(double *l, double *u, double *d, double c, int grid, int myNode, int totalNodes);


std::string getCurrentPath()
{
    char charCurrentPath[1024];
    std::string filePath;
    if (getcwd(charCurrentPath, sizeof(charCurrentPath)) != NULL)
    {
        filePath = std::string(charCurrentPath);
    } else
        perror("getcwd() error");


    return filePath;

}


int main(int argc, char *argv[])
{


    double maxGridValue = 50;
    double gridMinimumValue = -50;
    double acceleration = 1.75;
    double upperTimeBondary = 5.0;
    //int numberOfPoints = 1000;
    int numberOfPoints = 10000;
    double CourantNumber = 0.999;


    double *l, *u, *d, *y;

    double i;
    int j, procNumberOfPoints;
    int myRank, numOfNodes;
    double dt, dx, minimumXProcesor;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    double startTime = MPI_Wtime();

    dx = (maxGridValue - gridMinimumValue) / (numberOfPoints);
    dt = (CourantNumber * dx) / acceleration;

    procNumberOfPoints = getProcArraySize(myRank, numOfNodes, numberOfPoints);
    minimumXProcesor = gridMinimumValue + dx * procNumberOfPoints * myRank;

    l = new double[procNumberOfPoints];
    u = new double[procNumberOfPoints];
    d = new double[procNumberOfPoints];
    y = new double[procNumberOfPoints];

    double c = 0.25 * CourantNumber;
    LUPreparation(l, u, d, c, procNumberOfPoints, myRank, numOfNodes);

    double *initValue = new double[procNumberOfPoints];
    double *lastValue = new double[procNumberOfPoints];

    for (j = 0; j < procNumberOfPoints; j++)
    {

        initValue[j] = getInitialValue(minimumXProcesor, dx, j);
    }

    double *valueOfWave = new double[procNumberOfPoints];
    double currentTime = 0.0;

    while (upperTimeBondary > currentTime)
    {
        solve(procNumberOfPoints, dx, CourantNumber, initValue, valueOfWave, myRank, numOfNodes, l, u, d, y);
        for (int i = 0; i < procNumberOfPoints; i++)
        {
            initValue[i] = valueOfWave[i];
        }
        currentTime += dt;
    }

    if (myRank == 0)
    {
        double *finalSolution = new double[numberOfPoints];
        double *index = finalSolution;

        for (int i = 0; i < procNumberOfPoints; i++)
        {
            finalSolution[i] = initValue[i];
        }

        for (unsigned int i = 1; i < numOfNodes; i++)
        {
            MPI_Status myActualStatus;
            int itProcSize = getProcArraySize(i - 1, numOfNodes, numberOfPoints);
            index += itProcSize;
            MPI_Recv(index, getProcArraySize(i, numOfNodes, numberOfPoints), MPI_DOUBLE, i, 9, MPI_COMM_WORLD,
                     &myActualStatus);
        }

        std::cout << "Current path is: " << getCurrentPath() << std::endl;

        fstream fs;
        std::string tmpPath = getCurrentPath() + "/Res/x";
        const char *currentPath = tmpPath.c_str();
        fs.open(currentPath, fstream::out | fstream::trunc);
        double currentGridXValue = gridMinimumValue;

        for (int i = 0; i < numberOfPoints; i++)
        {
            fs << currentGridXValue << "," << finalSolution[i] << "," << calculateAnalyticalSolution(currentGridXValue,
                                                                                                     currentTime)
               << endl;
            currentGridXValue += dx;
        }

        double calculateNorm = 0.0;
        for (unsigned int i = 0; i < numberOfPoints; i++)
        {
            calculateNorm += fabs(
                    calculateAnalyticalSolution(gridMinimumValue + i * dx, currentTime - dt) - finalSolution[i]);
        }

        std::cout << "One norm is: " << calculateNorm << std::endl;
        fs << "Value of one norm" << calculateNorm << std::endl;

        std::cout << "Whole program execution time: " << MPI_Wtime() - startTime << std::endl;
        fs << "Whole program execution time: " << MPI_Wtime() - startTime << std::endl;

        fs.close();

        delete[] finalSolution;
    } else
    {
        MPI_Send(&initValue[0], procNumberOfPoints, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
    }

    delete[] l;
    delete[] u;
    delete[] d;
    delete[] y;
    delete[] initValue;
    delete[] lastValue;
    delete[] valueOfWave;

    MPI_Finalize();
    return 0;
}


void solve(int numOfPoints, double deltaX, double CFL, double *initData, double *finalWave, int myNode, int totalNodes,
           double *l, double *u, double *d, double *y)
{

    double tmp;
    MPI_Status status;
    double *q = finalWave;


    double lInitialBoundary = 0, rInitialBoundary = 0;

    long sRNode = (myNode + 1) % totalNodes;
    long rRNode = (myNode - 1 + totalNodes) % totalNodes;
    long sLNode = (myNode - 1 + totalNodes) % totalNodes;
    long rLNode = (myNode + 1) % totalNodes;

    tmp = initData[numOfPoints - 1];
    MPI_Sendrecv(&tmp, 1, MPI_DOUBLE, sRNode, TAG_EIGHT,
                 &lInitialBoundary, 1, MPI_DOUBLE, rRNode, TAG_EIGHT, MPI_COMM_WORLD, &status);

    tmp = initData[0];
    MPI_Sendrecv(&tmp, 1, MPI_DOUBLE, sLNode, TAG_SEVEN,
                 &rInitialBoundary, 1, MPI_DOUBLE, rLNode, TAG_SEVEN, MPI_COMM_WORLD, &status);

    for (int i = 1; i < numOfPoints - 1; i++)
    {
        q[i] = initData[i] - 0.25 * CFL * (initData[i + 1] - initData[i - 1]);
    }

    q[0] = initData[0] - 0.25 * CFL * (initData[1] - lInitialBoundary);
    q[numOfPoints - 1] = initData[numOfPoints - 1] - 0.25 * CFL * (rInitialBoundary - initData[numOfPoints - 2]);

    if (myNode == 0)
    {
        y[0] = q[0];
    } else
    {
        double LeftValueY;
        double leftValueL;

        MPI_Recv(&LeftValueY, 1, MPI_DOUBLE, myNode - 1, TAG_THREE, MPI_COMM_WORLD, &status);
        MPI_Recv(&leftValueL, 1, MPI_DOUBLE, myNode - 1, TAG_ONE, MPI_COMM_WORLD, &status);

        y[0] = q[0] - leftValueL * LeftValueY;
    }

    for (int i = 1; i < numOfPoints; i++)
        y[i] = q[i] - l[i - 1] * y[i - 1];


    if (myNode != totalNodes - 1)
    {
        tmp = y[numOfPoints - 1];
        MPI_Send(&tmp, 1, MPI_DOUBLE, myNode + 1, TAG_THREE, MPI_COMM_WORLD);
        tmp = l[numOfPoints - 1];
        MPI_Send(&tmp, 1, MPI_DOUBLE, myNode + 1, TAG_ONE, MPI_COMM_WORLD);
    }

    double *x = finalWave;

    if (myNode == totalNodes - 1)
    {
        x[numOfPoints - 1] = y[numOfPoints - 1] / d[numOfPoints - 1];
    } else
    {
        double valueOfXRight;
        MPI_Recv(&valueOfXRight, 1, MPI_DOUBLE, myNode + 1, TAG_SIX, MPI_COMM_WORLD, &status);

        x[numOfPoints - 1] = (y[numOfPoints - 1] - u[numOfPoints - 1] * valueOfXRight) / d[numOfPoints - 1];
    }

    for (int i = numOfPoints - 2; i >= 0; i--)
        x[i] = (y[i] - u[i] * x[i + 1]) / d[i];

    if (myNode != 0)
    {
        tmp = x[0];
        MPI_Send(&tmp, 1, MPI_DOUBLE, myNode - 1, TAG_SIX, MPI_COMM_WORLD);
    }

}

void LUPreparation(double *l, double *u, double *d, double c, int grid, int myNode, int totalNodes)
{
    int N = grid;

    double a = 1.0;
    double b = -c;
    double tmp;
    if (myNode == 0)
    {
        d[0] = 1.0;
    } else
    {

        MPI_Status status;
        MPI_Recv(&tmp, 1, MPI_DOUBLE, myNode - 1, TAG_FIVE, MPI_COMM_WORLD, &status);
        d[0] = a - tmp * c;
    }
    u[0] = c;
    for (int i = 0; i < N - 1; i++)
    {
        l[i] = b / d[i];
        d[i + 1] = a - l[i] * u[i];
        u[i + 1] = c;
    }
    l[N - 1] = b / d[N - 1];

    if (myNode != totalNodes - 1)
    {
        tmp = l[N - 1];
        MPI_Send(&tmp, 1, MPI_DOUBLE, myNode + 1, TAG_FIVE, MPI_COMM_WORLD);
    }
}

double getInitialValue(double totalNumberOfPointsX, double deltaX, int j)
{

    return calculateAnalyticalSolution(totalNumberOfPointsX + j * deltaX, 0.0);
}

double calculateAnalyticalSolution(double totalNumberOfPointsX, double timeValue)
{

    return 0.5 * exp(-(totalNumberOfPointsX - 1.75 * timeValue) * (totalNumberOfPointsX - 1.75 * timeValue));
}

int getProcArraySize(int myNode, int totalNodes, int totalNumberOfPoints)
{

    if (totalNodes == 1)
    {

        return totalNumberOfPoints;

    } else if (totalNumberOfPoints % totalNodes == 0)
    {

        return totalNumberOfPoints / totalNodes;

    } else
    {

        if (myNode == totalNodes - 1)
        {

            return (totalNumberOfPoints / totalNodes) + (totalNumberOfPoints % totalNodes);

        } else
        {

            return totalNumberOfPoints / totalNodes;
        }
    }

}

