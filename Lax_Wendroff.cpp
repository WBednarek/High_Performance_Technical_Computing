#include "Lax_Wendroff.h"


Lax_Wendroff::Lax_Wendroff() {
}

Lax_Wendroff::Lax_Wendroff(double xMin,
                           double xMax,
                           double time,
                           double numberOfSpacePoints,
                           double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time, numberOfSpacePoints, CFL),
                                         methodName("Lax_Wendroff") {
}


Lax_Wendroff::~Lax_Wendroff() {
}

void Lax_Wendroff::solve(int setNumber) {
    try {

        std::cout << "Lax_Wendroff upwind scheme solution runs and matrix is initialised\n";

        (*this).initializeSet(setNumber);
        Lax_WendroffResutls = Matrix((*this).getMatrix());


        double T1 = (CFL * (CFL + 1)) / 2;
        double T2 = 1 - (CFL * CFL);
        double T3 = (CFL * (CFL - 1)) / 2;


        for (int j = 0; j < numberOfTimePoints - 1; ++j) {
            for (int i = 1; i < numberOfSpacePoints - 1; ++i) {

                Lax_WendroffResutls[i][j + 1] = T1 * Lax_WendroffResutls[i - 1][j] + T2 * Lax_WendroffResutls[i][j] +
                                                T3 * Lax_WendroffResutls[i + 1][j];
            }

        }
        GeneralScheme::solve(setNumber);
        calculateNorms((*this).Lax_WendroffResutls);

    }

    catch (std::exception &e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }
}

Matrix Lax_Wendroff::getLax_WendroffdMatrix() {
    return Lax_WendroffResutls;
}

std::string Lax_Wendroff::getName() {
    return methodName;
}
