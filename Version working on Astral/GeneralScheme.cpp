#include <mpi.h>
#include "GeneralScheme.h"

/**
Default constructor
*/

GeneralScheme::GeneralScheme()
{

}

GeneralScheme::~GeneralScheme()
{

}


GeneralScheme::GeneralScheme(double xMin, double xMax, double time, double numberOfSpacePoints, double CFL)
        : xMin(xMin), xMax(xMax), time(time), numberOfSpacePoints(numberOfSpacePoints), CFL(CFL),
          isSetInitialised(false), isAnaliticalSolutionSolved(false), name("GeneralScheme")
{

    (*this).dt = (*this).calculateDtValue();
    (*this).dx = (*this).calculateDxValue();
    (*this).numberOfTimePoints = std::ceil((time) / (((*this).CFL * (*this).dx) / u));

    matrixOfResults = Matrix(numberOfSpacePoints, numberOfTimePoints);
    (*this).calculateDtValue();

}

double GeneralScheme::calculateDtValue()
{
    return (*this).dt = ((*this).CFL * (*this).dx) / u;
}


double GeneralScheme::calculateDxValue()
{
    return (*this).dx = (std::abs((*this).xMin) + std::abs((*this).xMax) + 1) / (*this).numberOfSpacePoints;
}


//double GeneralScheme::getDx() {
//    return dx;
//}

double GeneralScheme::initializationFunction(int numberOfSet, double functionValue)
{
    switch (numberOfSet)
    {
        case 1:
            return (MathFunctions::sign(functionValue) + 1);
            break;
        case 2:
            return std::exp((-1.0) * std::pow(functionValue, 2));
            break;
        default:
            std::cout << "There is no such choice using first condition set" << std::endl;
            return (MathFunctions::sign(functionValue) + 1);
    }


}

double GeneralScheme::solutionFunctionAnalytical(int numberOfSet, double actualSpaceValue, double actualTimeValue)
{
    switch (numberOfSet)
    {
        case 1:
            return 0.5 * (MathFunctions::sign(actualSpaceValue - 1.75 * actualTimeValue) + 1);
            break;
        case 2:
            return 0.5 * std::exp((-1.0) * std::pow(actualSpaceValue - 1.75 * actualTimeValue, 2));
            break;
        default:
            return 0.5 * (MathFunctions::sign(actualSpaceValue - 1.75 * actualTimeValue) + 1);
            break;
    }
}


//Error and norms calculation
void GeneralScheme::calculateNorms(Matrix &toCalculateError)
{
    errorVector.resize(toCalculateError.getNumOfRows());
    for (int i = 0; i < numberOfSpacePoints; ++i)
    {
        errorVector[i] =
                toCalculateError[i][numberOfTimePoints - 1] - (*this).matrixOfResults[i][numberOfTimePoints - 1];
    }
    /*
    Second approach - quite complicated
    (*this).normInfiniteValue = std::max_element(errorVector.begin(), errorVector.end(), MathFunctions::compareTwoAbsElements);
    normInfiniteValue1 = std::distance(errorVector.begin(), normInfiniteValue);
    double normInf = errorVector.at(normInfiniteValue1);
    */
    normInfiniteValue = 0;
    for (unsigned int i = 0; i < errorVector.size(); ++i)
    {
        if (abs(normInfiniteValue) < abs(errorVector.at(i)))
        {
            normInfiniteValue = abs(errorVector.at(i));
        }

        (*this).normOneValue += std::abs(errorVector[i]);
        (*this).normTwoValue += std::pow(std::abs(errorVector[i]), 2);
    }


    normTwoValue = std::sqrt(normTwoValue);
    normOneValue = normOneValue / numberOfSpacePoints;
    normTwoValue = normTwoValue / numberOfSpacePoints;

    std::vector<double> norms;
    norms.push_back(normInfiniteValue);
    norms.push_back(normOneValue);
    norms.push_back(normTwoValue);
    norms.push_back(dx);

    toCalculateError.resizeMat(numberOfSpacePoints, numberOfTimePoints + 1);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if (myRank == 0)
    {
        for (unsigned int i = 0; i < norms.size() - 1; ++i)
        {
            //Adding norms results to last matrix colunm
            std::cout << "Norm " << i + 1 << ": " << norms.at(i) << std::endl;
            toCalculateError[i][numberOfTimePoints] = norms.at(i);
        }
    }


}


void GeneralScheme::put_timeValues()
{

    double actualValue = 0;
    for (int i = 0; i < numberOfTimePoints; ++i)
    {
        matrixOfResults[0][i] = actualValue;
        actualValue += (*this).dt;
    }

}

std::string GeneralScheme::getName()
{
    return name;
}

void GeneralScheme::initializeSet(int setNumber)
{


    try
    {

        double actualValue = xMin;
        for (int i = 0; i < numberOfSpacePoints; ++i)
        {
            matrixOfResults[i][0] = (1.0 / 2.0) * (*this).initializationFunction(setNumber, actualValue);
            actualValue += (*this).dx;
        }


        if (setNumber == 1)
        {
            for (int i = 0; i < numberOfTimePoints; ++i)
            {
                matrixOfResults[0][i] = 0;
                matrixOfResults[numberOfSpacePoints - 1][i] = 1;
            }
        } else
        {
            for (int i = 0; i < numberOfTimePoints; ++i)
            {
                matrixOfResults[0][i] = 0;
                matrixOfResults[numberOfSpacePoints - 1][i] = 0;
            }
        }

        (*this).isSetInitialised = true;

    }
    catch (std::exception &e)
    {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }


}


//Analytical scheme solving
void GeneralScheme::solve(int numberOfBoundaryConditionSet)
{
    (*this).initializeSet(numberOfBoundaryConditionSet);
    try
    {
        if ((*this).isSetInitialised == true)
        {
            std::cout << "Analytical solution runs and matrix is initialised\n";

            //Variables hold values below 0. Thanks to that negative values could be passed to sign function, it makes loop iteration easier.
            double actualSpaceValue = xMin;
            //Variable assinged to dt because time at 0 point is initialised in function initializeSet()
            double actualTimeValue = dt;
            for (int i = 1; i < numberOfSpacePoints; ++i)
            {
                for (int j = 1; j < numberOfTimePoints; ++j)
                {
                    matrixOfResults[i][j] = solutionFunctionAnalytical(numberOfBoundaryConditionSet, actualSpaceValue,
                                                                       actualTimeValue);
                    actualTimeValue += dt;

                }
                actualTimeValue = dt;
                actualSpaceValue += dx;
            }

            isAnaliticalSolutionSolved = true;
        } else
        {
            std::cout << "Matrix is not initialised\n";
        }

    }

    catch (std::exception &e)
    {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }


}


std::vector<double> GeneralScheme::getResults()
{
    return matrixOfResults.getColumn(numberOfTimePoints - 1);
}

const Matrix &GeneralScheme::getMatrix() const
{
    return matrixOfResults;
}






