

#include "CrankParallel.h"


CrankParallel::CrankParallel(double xMin,
                             double xMax,
                             double time,
                             double numberOfSpacePoints,
                             double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                        numberOfSpacePoints, CFL),
                                           methodName("CrankParallelScheme")
{


    /**
       * from CN:
       * C * f(n+1)(i+1) + A * f(n+1)(i) + B * f(n+1)(i-1) = ...
       */
    A = 1;
    C = CFL / 4;
    B = -CFL / 4;

    BR = 1;
    CR = CFL / 4;
    AR = -CFL / 4;

}

CrankParallel::~CrankParallel()
{
}


void CrankParallel::ThomasAlgorithm_P_LUDecomposition(int myCurrNode, int numberOfNodes, int N, double *b,
                                                      double *a, double *c, double *l, double *d)
{
    int i;
    int rowsOfLocal, offsetLocal;
    double S[2][2], T[2][2], s1TMPVal, s2TMPVal;
    MPI_Status myStatus;

    for (i = 0; i < N; i++)
        l[i] = d[i] = 0.0;

    S[0][0] = S[1][1] = 1.0;
    S[1][0] = S[0][1] = 0.0;
    rowsOfLocal = (int) floor(N / numberOfNodes);
    offsetLocal = myCurrNode * rowsOfLocal;
// Form local products of R_k matrices
    if (myCurrNode == 0)
    {
        s1TMPVal = a[offsetLocal] * S[0][0];
        S[1][0] = S[0][0];
        S[1][1] = S[0][1];
        S[0][1] = a[offsetLocal] * S[0][1];
        S[0][0] = s1TMPVal;
        for (i = 1; i < rowsOfLocal; i++)
        {
            s1TMPVal = a[i + offsetLocal] * S[0][0] -
                       b[i + offsetLocal - 1] * c[i + offsetLocal - 1] * S[1][0];
            s2TMPVal = a[i + offsetLocal] * S[0][1] -
                       b[i + offsetLocal - 1] * c[i + offsetLocal - 1] * S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1TMPVal;

            S[0][1] = s2TMPVal;
        }
    } else
    {
        for (i = 0; i < rowsOfLocal; i++)
        {
            s1TMPVal = a[i + offsetLocal] * S[0][0] -
                       b[i + offsetLocal - 1] * c[i + offsetLocal - 1] * S[1][0];
            s2TMPVal = a[i + offsetLocal] * S[0][1] -
                       b[i + offsetLocal - 1] * c[i + offsetLocal - 1] * S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1TMPVal;
            S[0][1] = s2TMPVal;
        }
    }

    for (i = 0; i <= log2(numberOfNodes); i++)
    {
        if (myCurrNode + pow(2, i) < numberOfNodes)
            MPI_Send(S, 4, MPI_DOUBLE, int(myCurrNode + pow(2, i)), 0,
                     MPI_COMM_WORLD);
        if (myCurrNode - pow(2, i) >= 0)
        {
            MPI_Recv(T, 4, MPI_DOUBLE, int(myCurrNode - pow(2, i)), 0,
                     MPI_COMM_WORLD, &myStatus);
            s1TMPVal = S[0][0] * T[0][0] + S[0][1] * T[1][0];
            S[0][1] = S[0][0] * T[0][1] + S[0][1] * T[1][1];
            S[0][0] = s1TMPVal;
            s1TMPVal = S[1][0] * T[0][0] + S[1][1] * T[1][0];
            S[1][1] = S[1][0] * T[0][1] + S[1][1] * T[1][1];
            S[1][0] = s1TMPVal;
        }
    }

    d[offsetLocal + rowsOfLocal - 1] = (S[0][0] + S[0][1]) /
                                       (S[1][0] + S[1][1]);
    if (myCurrNode == 0)
    {
        MPI_Send(&d[offsetLocal + rowsOfLocal - 1], 1, MPI_DOUBLE,
                 1, 0, MPI_COMM_WORLD);
    } else
    {
        MPI_Recv(&d[offsetLocal - 1], 1, MPI_DOUBLE, myCurrNode - 1, 0,
                 MPI_COMM_WORLD, &myStatus);
        if (myCurrNode != numberOfNodes - 1)

            MPI_Send(&d[offsetLocal + rowsOfLocal - 1], 1, MPI_DOUBLE,
                     myCurrNode + 1, 0, MPI_COMM_WORLD);
    }


    if (myCurrNode == 0)
    {
        l[0] = 0;
        d[0] = a[0];
        for (i = 1; i < rowsOfLocal - 1; i++)
        {
            l[offsetLocal + i] = b[offsetLocal + i - 1] /
                                 d[offsetLocal + i - 1];
            d[offsetLocal + i] = a[offsetLocal + i] -
                                 l[offsetLocal + i] * c[offsetLocal + i - 1];
        }
        l[offsetLocal + rowsOfLocal - 1] = b[offsetLocal + rowsOfLocal - 2] /
                                           d[offsetLocal + rowsOfLocal - 2];
    } else
    {
        for (i = 0; i < rowsOfLocal - 1; i++)
        {
            l[offsetLocal + i] = b[offsetLocal + i - 1] /
                                 d[offsetLocal + i - 1];
            d[offsetLocal + i] = a[offsetLocal + i] -
                                 l[offsetLocal + i] * c[offsetLocal + i - 1];
        }
        l[offsetLocal + rowsOfLocal - 1] = b[offsetLocal + rowsOfLocal - 2] /
                                           d[offsetLocal + rowsOfLocal - 2];
    }


    if (myCurrNode > 0)
        d[offsetLocal - 1] = 0;
    // Distribute d_k and l_k to all processes
    double *tmp = new double[N];
    for (i = 0; i < N; i++)
        tmp[i] = d[i];
    MPI_Allreduce(tmp, d, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (i = 0; i < N; i++)
        tmp[i] = l[i];
    MPI_Allreduce(tmp, l, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] tmp;

}

/**
 *
 * @param N
 * @param lStore
 * @param dStore
 * @param c
 * @param x output
 * @param q
 */
void CrankParallel::ThomasAlgorithm_P_solve(int N, double *lStore, double *dStore, double *c, double *q, double *x)
{
    int i;
    double *y = new double[N];

    for (i = 0; i < N; i++)
        y[i] = 0.0;

    /* Forward Substitution [L][y] = [q] */
    y[0] = q[0];
    for (i = 1; i < N; i++)
        y[i] = q[i] - lStore[i] * y[i - 1];
    /* Backward Substitution [U][x] = [y] */
    x[N - 1] = y[N - 1] / dStore[N - 1];
    for (i = N - 2; i >= 0; i--)
        x[i] = (y[i] - c[i] * x[i + 1]) / dStore[i];

    delete[] y;
    return;
}

/*
 * from CN:
 * ... = AR * f(n)(i+1) + BR * f(n)(i) + CR * f(n)(i-1)
 *
 * @param y input (n)
 * @param q output (n+1)
 */
void CrankParallel::CrankNicholson_RS(int N, double AR, double BR, double CR, double *y, double *q)
{
    for (int i = 1; i < N - 1; ++i)
    {
        q[i] = AR * y[i - 1] + BR * y[i] + CR * y[i + 1];
    }
}

void CrankParallel::solve(int setNumber)
{
    try
    {

        std::cout << methodName << " solution runs and matrix is initialised\n";
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
        crankResutls = Matrix(numberOfSpacePoints, numberOfTimePoints);
        double timeOfStart = MPI_Wtime();

        lastNode = numOfProc - 1;
        if (numberOfSpacePoints % numOfProc != 0)
        {
            throw "error";
        }

        //only for numOfProc = 2^n; n = 1, 2, ...


        if (myRank == 0)
        {
            double actualValue = xMin;

            for (int i = 0; i < numberOfSpacePoints; ++i)
            {

                actualValue = (i * (*this).dx) + xMin;
                crankResutls[i][0] = (1.0 / 2.0) * (*this).initializationFunction(setNumber, actualValue);

            }

            for (int i = 0; i < numberOfTimePoints; ++i)
            {
                crankResutls[0][i] = 0;
            }

            for (int i = 0; i < numberOfTimePoints; ++i)
            {

                crankResutls[numberOfSpacePoints - 1][i] = 0;
            }

        }

        double *a = new double[numberOfSpacePoints];
        double *b = new double[numberOfSpacePoints];
        double *c = new double[numberOfSpacePoints];
        double *d = new double[numberOfSpacePoints];
        double *l = new double[numberOfSpacePoints];
        double *tmpData = new double[numberOfSpacePoints];
        double *tmpPrev = new double[numberOfSpacePoints];
        double *tmpNext = new double[numberOfSpacePoints];


        b[0] = 0.0;
        c[numberOfSpacePoints - 1] = 0.0;

        /**
         * from CN:
         * C * f(n+1)(i+1) + A * f(n+1)(i) + B * f(n+1)(i-1) = ...
         */
        for (int i = 0; i < numberOfSpacePoints; ++i)
        {
            a[i] = A;
            if (i != 0) b[i] = B;
            if (i != numberOfSpacePoints - 1) c[i] = C;
            l[i] = 0.0;
            d[i] = 0.0;
        }

        ThomasAlgorithm_P_LUDecomposition(myRank, numOfProc, numberOfSpacePoints, b, a, c, l, d);

        for (int j = 0; j < numberOfTimePoints - 1; ++j)
        {
            for (int i = 0; i < numberOfSpacePoints; ++i)
                tmpPrev[i] = crankResutls[i][j];

            CrankNicholson_RS(numberOfSpacePoints, AR, BR, CR, tmpPrev, tmpData);
            ThomasAlgorithm_P_solve(numberOfSpacePoints, l, d, c, tmpData, tmpNext);

            for (int i = 0; i < numberOfSpacePoints; ++i)
                crankResutls[i][j + 1] = tmpNext[i];
        }

        delete (a);
        delete (b);
        delete (c);
        delete (d);
        delete (l);
        delete (tmpData);
        delete (tmpPrev);
        delete (tmpNext);

        //GeneralScheme::solve(setNumber);
        calculateNorms((*this).crankResutls);
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

std::vector<double> CrankParallel::getCrankResutls()
{
    return crankResutls.getColumn(numberOfTimePoints - 1);
}

const std::string &CrankParallel::getMethodName() const
{
    return methodName;
}
