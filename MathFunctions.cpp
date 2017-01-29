#include "MathFunctions.h"


MathFunctions::MathFunctions() {
}


MathFunctions::~MathFunctions() {
}

/**
@brief Static function sign

Fucntion gives an output of sign function
1. -1 When input walue < 0
2.  0 when input value = 0
3.  1 When input value > 0

@param Input value to sign function
@return result of sign function
*/
int MathFunctions::sign(double x) {
    if (x < 0) {
        return -1;
    } else if (x == 0) {
        return 0;
    } else {
        return 1;
    }

}

bool MathFunctions::compareTwoAbsElements(double first, double second) {
    return (std::abs(first) < std::abs(second));
}
