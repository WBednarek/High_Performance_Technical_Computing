//
// Created by wik on 05.02.17.
//
#include <mpi.h>
#include <math.h>
#include "GeneralScheme.h"

#ifndef HIGH_PERFORMANCE_TECHNICAL_COMPUTING_CRANKPARALLEL_H
#define HIGH_PERFORMANCE_TECHNICAL_COMPUTING_CRANKPARALLEL_H


class CrankParallel : public GeneralScheme
{


    Matrix crankResutls;
    std::string methodName;
public:
    const std::string &getMethodName() const;

private:
    int myRank;
    int numOfProc;
    double lastNode;
    double localLimit;
    double A, AR;
    double B, BR;
    double C, CR;

public:

    CrankParallel(double xMin,
                  double xMax,
                  double time,
                  double numberOfSpacePoints,
                  double CFL);

    ~CrankParallel();

    void solve(int setNumber);

    void ThomasAlgorithm_P_LUDecomposition(int myCurrNode, int numberOfNodes, int N, double *b,
                                           double *a, double *c, double *l, double *d);

    void ThomasAlgorithm_P_solve(int N, double *l, double *d, double *c, double *x, double *q);

    void CrankNicholson_RS(int N, double AR, double BR, double CR, double *y, double *q);


    std::vector<double> getCrankResutls();


};


#endif //HIGH_PERFORMANCE_TECHNICAL_COMPUTING_CRANKPARALLEL_H
