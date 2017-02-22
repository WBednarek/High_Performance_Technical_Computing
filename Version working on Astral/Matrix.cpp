#include "Matrix.h"

/**
Matrix class allows to  inspired by Dr Peter Sherar's Matrix class
*/

typedef std::vector<std::vector<double> > mat;

/**
Default Constructor
*/
Matrix::Matrix() : mat()
{}

/**
Constructor
*/
Matrix::Matrix(int numOfRows, int numOfColumns) : mat()
{


    (*this).resize(numOfRows);

    for (int i = 0; i < numOfRows; ++i)
    {
        (*this)[i].resize(numOfColumns);
    }

}

//Copy constructor
Matrix::Matrix(const Matrix &m) : std::vector<std::vector<double> >()
{
    // set the size of the rows
    (*this).resize(m.size());
    // set the size of the columns
    std::size_t i;
    for (i = 0; i < m.size(); i++) (*this)[i].resize(m[0].size());

    // copy the elements
    for (int i = 0; i < m.getNumOfRows(); i++)
        for (int j = 0; j < m.getNumOfColumns(); j++)
            (*this)[i][j] = m[i][j];
}

int Matrix::getNumOfRows() const
{
    return (*this).size();
}

int Matrix::getNumOfColumns() const
{
    return (*this)[0].size();
}


std::vector<double> &Matrix::getRow(int rowNumber)
{
    static std::vector<double> tmp;

    for (int i = 0; i < (*this).getNumOfColumns(); ++i)
    {
        tmp.push_back((*this)[rowNumber][i]);
    }

    return tmp;
}


std::vector<double> &Matrix::getColumn(int columnNumber)
{


    for (int i = 0; i < (*this).getNumOfRows(); ++i)
    {
        column.push_back((*this)[i][columnNumber]);
    }

    return column;
}


std::ostream &operator<<(std::ostream &os, Matrix &mat)
{

    for (int i = 0; i < mat.getNumOfRows(); ++i)
    {
        for (int j = 0; j < mat.getNumOfColumns(); ++j)
        {
            os << mat[i][j] << " ";
        }
        os << std::endl;
    }

    return os;
}


std::ofstream &operator<<(std::ofstream &ofs, const Matrix &m)
{
    //put matrix rownumber in first line (even if it is zero)
    //ofs << "dt" << std::endl;
    //put matrix columnnumber in second line (even if it is zero)
    //ofs << m.getNumOfColumns() << std::endl;
    //put data in third line (if size==zero nothing will be put)
    for (int i = 0; i < m.getNumOfRows(); i++)
    {
        for (int j = 0; j < m.getNumOfColumns(); j++) ofs << m[i][j] << "\t";
        ofs << std::endl;
    }
    return ofs;
}


Matrix &Matrix::operator=(const Matrix &m)
{
    (*this).resize(m.size());
    std::size_t i;
    std::size_t j;
    for (i = 0; i < m.size(); i++) (*this)[i].resize(m[0].size());

    for (i = 0; i < m.size(); i++)
        for (j = 0; j < m[0].size(); j++)
            (*this)[i][j] = m[i][j];
    return *this;
}


void Matrix::resizeMat(int numOfRows, int numOfColumns)
{


    (*this).resize(numOfRows);

    for (int i = 0; i < numOfRows; ++i)
    {
        (*this)[i].resize(numOfColumns);
    }

}

