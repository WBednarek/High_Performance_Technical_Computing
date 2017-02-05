#include <iostream>
#include <vector>
#include <stdio.h>
#include <iterator>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <map>
#include <iomanip>
#include <mpi/mpi.h>
#include "MathFunctions.h"
#include "GeneralScheme.h"
#include "ExplicitUpwindScheme.h"
#include "Display.h"
#include "ExplicitUpwindParallel.h"
#include "ImplicitParallel.h"
#include "ImplicitUpwindScheme.h"


using std::vector;
using std::cin;
using std::cout;
using std::endl;
//Own display namespace with user friendly displaying functions
using namespace display;


template<typename CharT>
class DecimalSeparator : public std::numpunct<CharT> {
public:
    DecimalSeparator(CharT Separator)
            : m_Separator(Separator) {}

protected:
    CharT do_decimal_point() const {
        return m_Separator;
    }

private:
    CharT m_Separator;
};


template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, "\n"));
    return out;
}

std::string getCurrentPath() {
    char charCurrentPath[1024];
    std::string path;
    if (getcwd(charCurrentPath, sizeof(charCurrentPath)) != NULL) {
        path = std::string(charCurrentPath);
    } else
        perror("getcwd() error");


    return path;

}


void runSchemes(int numberOfBoundaryConditionSet, vector<double> initialSettings, std::string typeOfExtension) {
    //Create Results folder first
    std::string path = getCurrentPath() + "/Results/";


    GeneralScheme general = GeneralScheme(initialSettings[0], initialSettings[1], initialSettings[2],
                                          initialSettings[3], initialSettings[4]);
    general.solve(numberOfBoundaryConditionSet);

    ExplicitUpwindScheme explicitScheme(initialSettings[0], initialSettings[1], initialSettings[2], initialSettings[3],
                                        initialSettings[4]);
    explicitScheme.solve(numberOfBoundaryConditionSet);

    ExplicitUpwindParallel parallel(initialSettings[0], initialSettings[1], initialSettings[2], initialSettings[3],
                                    initialSettings[4]);
    parallel.solve(numberOfBoundaryConditionSet);


    ImplicitUpwindScheme implicitUpwindScheme(initialSettings[0], initialSettings[1], initialSettings[2],
                                              initialSettings[3], initialSettings[4]);
    implicitUpwindScheme.solve(numberOfBoundaryConditionSet);

    ImplicitParallel implicitParallel(initialSettings[0], initialSettings[1], initialSettings[2], initialSettings[3],
                                      initialSettings[4]);
    implicitParallel.solve(numberOfBoundaryConditionSet);






    /*
     *


     osImplicitScheme.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));

     std::string implicitUpwindSchemeFileName =
            path + getInitialBoundaryConditionName(numberOfBoundaryConditionSet) + "_" +
            implicitUpwindScheme.getName() + "Results_t=" + streams[2].str() + "_points=" + streams[3].str() + "_CFL=" +
            streams[4].str() + typeOfExtension;
    const char *CharImplicitUpwindSchemeFileName = implicitUpwindSchemeFileName.c_str();

    osImplicitScheme.open(CharImplicitUpwindSchemeFileName);

    osImplicitScheme << implicitUpwindScheme.getLastImplicitMatrixColumn();

     osImplicitScheme.close();
    */


    //Generating open files
    std::ofstream osGeneralScheme;
    std::ofstream osExplicitScheme;
    std::ofstream osImplicitScheme;


    std::ofstream osParallel;
    std::ofstream osImplicitParallel;



    //Operation helps to plot charts in programs such as Exel. Setting type of decimal separator depending on current geographical location. In some countries comma in default separator in numbers in others dot
    osGeneralScheme.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));
    osExplicitScheme.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));
    osImplicitScheme.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));


    osParallel.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));
    osImplicitParallel.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));



    std::stringstream streams[initialSettings.size()];

    for (unsigned int i = 0; i < initialSettings.size(); ++i) {
        streams[i] << (int) initialSettings[i];
    }

    //Because I use c++98 I need to convert strings to chars
    std::string generalSchemeFileName =
            path + getInitialBoundaryConditionName(numberOfBoundaryConditionSet) + "_" + general.getName() +
            "Results_t=" + streams[2].str() + "_points=" + streams[3].str() + "_CFL=" + streams[4].str() +
            typeOfExtension;
    const char *CharGeneralSchemeFileName = generalSchemeFileName.c_str();

    std::string UpwindSchemeFileName =
            path + getInitialBoundaryConditionName(numberOfBoundaryConditionSet) + "_" + explicitScheme.getName() +
            "Results_t=" + streams[2].str() + "_points=" + streams[3].str() + "_CFL=" + streams[4].str() +
            typeOfExtension;
    const char *CharUpwindSchemeFileName = UpwindSchemeFileName.c_str();

    std::string ExplicitParallel =
            path + getInitialBoundaryConditionName(numberOfBoundaryConditionSet) + "_" + parallel.getName() +
            "Results_t=" + streams[2].str() + "_points=" + streams[3].str() + "_CFL=" + streams[4].str() +
            typeOfExtension;
    const char *CharExplicitParallelFileName = ExplicitParallel.c_str();

    std::string implicitParallelFileName =
            path + getInitialBoundaryConditionName(numberOfBoundaryConditionSet) + "_" + implicitParallel.getName() +
            "Results_t=" + streams[2].str() + "_points=" + streams[3].str() + "_CFL=" + streams[4].str() +
            typeOfExtension;
    const char *charImplicitParallelFileName = implicitParallelFileName.c_str();

    std::string implicitFileName =
            path + getInitialBoundaryConditionName(numberOfBoundaryConditionSet) + "_" +
            implicitUpwindScheme.getName() +
            "Results_t=" + streams[2].str() + "_points=" + streams[3].str() + "_CFL=" + streams[4].str() +
            typeOfExtension;
    const char *charImplicitFileName = implicitFileName.c_str();




    //std::string path = "C:/Users/Domowy/Desktop/Results/";
    //Open/create file with selected extension. It clold be for instance exel files extensions (.xls; .xlsx).
    //Single
    osGeneralScheme.open(CharGeneralSchemeFileName);
    osExplicitScheme.open(CharUpwindSchemeFileName);
    osImplicitScheme.open(charImplicitFileName);

    //Parallel
    osParallel.open(CharExplicitParallelFileName);
    osImplicitParallel.open(charImplicitParallelFileName);



    //Saving schemes calculated results

    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank == 0) {
        osGeneralScheme << general.getResults();

        //Explicit
        osExplicitScheme << explicitScheme.getLastExplicitMatrixColumn();
        osParallel << parallel.getGatherResults();

        //Implicit
        osImplicitParallel << implicitParallel.getImplicitParallelResults();
        osImplicitScheme << implicitUpwindScheme.getResults();



        //Closing all opened streams at the end
        osGeneralScheme.close();
        osExplicitScheme.close();
        osImplicitScheme.close();

        osImplicitParallel.close();
        osParallel.close();

    }




}


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    double timeOfStart;
    timeOfStart = MPI_Wtime();
    double timeOfEnd;


    //Number of boundary condition set. 1 for sign boundary set type ; 2 for exp boundary set type
//    vector<int> setNumber;
//    setNumber.push_back(1);
    int setNumber = 2;


    //Initial setings values are respectively: xMin, xMax, time, number of spacePoints, CFL value


    cout << "Actual parameters are: ";

    //Extension type of file which storing results of schemes computation. It could be Exel (.xls; .xlsx) file type for instance.
    std::string typeOfExtension = ".xls";

    //Running program using above settings for all initial boundary types

    vector<double> courantNumberSet;
    courantNumberSet.push_back(0.5);
    double courantNumber = 0.5;

    vector<double> pointsSet;
    pointsSet.push_back(10000);
    double numOfPoints = 10000;


    double simulationTime = 5;
/*
    for (int v = 0; v < setNumber.size(); ++v)
    {
        for (int j = 0; j < pointsSet.size(); ++j)
        {

            for (int i = 0; i < sizeof(timeSet) / sizeof(timeSet[0]); ++i)
            {

                for (int k = 0; k < courantNumberSet.size(); ++k)
                {
                    vector<double> initialSettings;
                    initialSettings.push_back(-50);
                    initialSettings.push_back(50);
                    initialSettings.push_back(simulationTime);
                    initialSettings.push_back(numOfPoints);
                    initialSettings.push_back(courantNumber);

                    runSchemes(setNumber[v], initialSettings, typeOfExtension);


                }

            }


        }


    }
*/

    vector<double> initialSettings;
    initialSettings.push_back(-50);
    initialSettings.push_back(50);
    initialSettings.push_back(simulationTime);
    initialSettings.push_back(numOfPoints);
    initialSettings.push_back(courantNumber);

    runSchemes(setNumber, initialSettings, typeOfExtension);
    timeOfEnd = MPI_Wtime();
    double programExecutionTime = timeOfEnd - timeOfStart;

    cout << "Total time of program execution is: " << programExecutionTime << endl;
    MPI_Finalize();
    return 0;

}


//system("pause");

	