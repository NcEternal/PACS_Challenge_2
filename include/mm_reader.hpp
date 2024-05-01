#ifndef HH_MM_READER_HH
#define HH_MM_READER_HH

#include "Matrix.hpp"

/*Reads a "MatrixMarket matrix coordinate real general" file and returns the matrix in it*/
template<algebra::Ordering OrderType = algebra::Ordering::RowMajor>
algebra::Matrix<double, OrderType> mm_reader(const std::string& filename);

extern template algebra::Matrix<double, algebra::Ordering::RowMajor> mm_reader(const std::string& filename);
extern template algebra::Matrix<double, algebra::Ordering::ColumnMajor> mm_reader(const std::string& filename);

#endif