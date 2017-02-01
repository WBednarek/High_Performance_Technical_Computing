#include <iterator>
#include "Display.h"


/*
void  display::displayVector(std::vector<double> vector)
{
	for (auto v : vector)
	{
		std::cout << "\n" << v;
	}
	std::cout << std::endl;

}
*/
std::string display::getInitialBoundaryConditionName(int &numberOfBoundaryCondition) {
    switch (numberOfBoundaryCondition) {
        case 1:
            return "Sign";
            break;
        case 2:
            return "Exp";
            break;
        default:
            return "Sign";
            break;
    }

}

//Overloaded operator << for easy loading vector to file


