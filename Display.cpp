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
            break;
    }

}

//Overloaded operator << for easy loading vector to file
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, "\n"));
    return out;
}