#include <vector>
#include <ostream>
#include <iostream>  
#include <fstream>  
#include <memory>


class Matrix : private std::vector<std::vector<double> > {


public:
    using std::vector<std::vector<double> >::operator[];

    Matrix();

    //std::shared_ptr < std::vector<double> > tmpVec;
    Matrix(int numOfRows, int numOfColumns);

    Matrix(const Matrix &m);

    /**
    @brief Method returns number of rows in Matrix

    @return Number of rows in Matrix
    */
    int getNumOfRows() const;

    /**
    @brief Method returns number of columns in Matrix

    @return Number of columns in Matrix
    */
    int getNumOfColumns() const;

    /**
    @brief Method returns selected row as vector

    @return Returns selected row as vector
    */
    std::vector<double> &getRow(int rowNumber);

    /**
    @brief Method returns selected column as vector

    @return Selected column as vector
    */
    std::vector<double> &getColumn(int columnnumber);

    friend std::ostream &operator<<(std::ostream &os, Matrix &mat);

    friend std::ofstream &operator<<(std::ofstream &ofs,
                                     const Matrix &m);

    Matrix &operator=(const Matrix &m);

    /**
    @brief Method returns selected column as vector

    @param New rows quantity for Matrix resize

    @param New columns quantity for Matrix resize
    @return Selected column as vector
    */
    void resizeMat(int numOfRows, int numOfColumns);

};