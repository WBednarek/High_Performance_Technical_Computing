#include "Richtmyer_multi_step.h"


Richtmyer_multi_step::Richtmyer_multi_step() {
}

Richtmyer_multi_step::Richtmyer_multi_step(double xMin,
                                           double xMax,
                                           double time,
                                           double numberOfSpacePoints,
                                           double CFL) : GeneralScheme::GeneralScheme(xMin, xMax, time,
                                                                                      numberOfSpacePoints, CFL),
                                                         methodName("RichtmyerMultiStepScheme") {
}


Richtmyer_multi_step::~Richtmyer_multi_step() {
}

void Richtmyer_multi_step::solve(int setNumber) {
    try {

        std::cout << "Richtmyer_multi_step scheme solution runs and matrix is initialised\n";

        //Preparing initial data with choosen setNumber: 1 - sign set; 2 - exp set
        (*this).initializeSet(setNumber);
        RichtmyerResutls = Matrix((*this).getMatrix());
        halfStepRichtmyer = Matrix((*this).getMatrix());


        //Calculating Richtmyer_multi_step scheme cooeficients before loop for code claryti and performance profit
        double coef1 = 0.5 * (1 - (CFL / 2));
        double coef2 = 0.5 * (1 + (CFL / 2));
        double coef3 = CFL / 2;



        //Main time loop iterating for each time point
        for (int j = 0; j < numberOfTimePoints - 1; ++j) {
            //First it is needed to calculate half step according to Richtmyer equation. Those results will be used in final computation.
            for (int i = 1; i < numberOfSpacePoints - 1; ++i) {

                halfStepRichtmyer[i][j] = coef1 * halfStepRichtmyer[i + 1][j] + coef2 * halfStepRichtmyer[i - 1][j];
            }

            //In this loop previosly calculated half step is used
            for (int i = 1; i < numberOfSpacePoints - 1; ++i) {
                RichtmyerResutls[i][j + 1] =
                        RichtmyerResutls[i][j] - coef3 * (halfStepRichtmyer[i + 1][j] - halfStepRichtmyer[i - 1][j]);
            }

        }
        GeneralScheme::solve(setNumber);
        calculateNorms((*this).RichtmyerResutls);

    }

    catch (std::exception &e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }

}

Matrix Richtmyer_multi_step::getRichtmyer_multi_stepdMatrix() {
    return RichtmyerResutls;
}

std::string Richtmyer_multi_step::getName() {
    return methodName;
}