# Challenge 2 A sparse matrix #

Solution to the $2^{nd}$ PACS challenge. Implementation of a sparse matrix class.

The makefile produces an executable called `test`, which collects multiple tests of the class functionalities.


## The Matrix Class ##

The class `Matrix` is implemented in the `algebra` namespace. It is a template class with parameters
`ValueType`, which must be a numeric (integral/floating point/complex<T>), and `OrderType`, which can be
either `Ordering::RowMajor` or `Ordering::ColumnMajor`. <br>
The class implements the following methods:

* A default constructor;

* A constructor that takes as arguments the number of rows and columns of the matrix;

* `ValueType& operator() (std::size_t i, std::size_t j)`: the random access operator, which takes the indices of the element as input arguments.
Out of bounds access will throw a `std::runtime_error`, as will attempting to access elements equal to 0 while the matrix is compressed;

* `ValueType operator() (std::size_t i, std::size_t j) const`: the constant version of the random access operator, which takes the indices of the element as input 
arguments. Out of bounds access will return `std::numeric_limits<ValueType>::quiet_NaN()`;

* `void compress()`: puts the matrix in its compressed state;

* `void uncompress()`: puts the matrix in its uncompressed state;

* `void resize(std::size_t new_rows, std::size_t new_cols)`: uncompresses and resizes the matrix to match the new dimensions;

* `template <Norm NormType> double norm() const`: a template method to calculate the norm of the matrix. The default template parameter is Norm::Frobenius. 
Other supported norms are Norm::One and Norm::Infinity;

* `bool is_compressed() const`: returns the current state of the matrix;

* `std::size_t rows() const`: returns the current number of rows of the matrix;

* `std::size_t cols() const`: returns the current number of columns of the matrix;  

* `get_row(std:.size_t i) const`: returns the i-th row as a `std::map<std::size_t, ValueType>`, with the key in this map representing the element's column;

* `get_col(std:.size_t j) const`: returns the j-th column as a `std::map<std::size_t, ValueType>`, with the key in this map representing the element's row.

The class also implements the following operators:

* `operator*`: Matrix-Vector multiplication (with a `std::vector`) or Matrix-Matrix multiplication (with another `Matrix`). The 2 objects being multiplied 
must have compatible sizes and share `ValueType`. `ValueType` is considered shared even if one of the 2 objects uses `std::complex<ValueType>`
(example: `Matrix<double, Ordering::RowMajor> * std::vector<std::complex<double>>`). <br> 
The result will be a `std::vector<T>` with the appropriate T in the case of Matrix-Vector multplication and a `Matrix<T, O>` with the appropriate `ValueType` T 
and the same `OrderType` O as the matrix to the left of the operator in the case of Matrix-Matrix multiplication;

* `operator<<`: the streaming operator to output the matrix to a stream object;

* `operator==`: a rudimentary comparison operator. Can only compare matrices with the exact same `ValueType` and `OrderType` and only if they are in the same
compression state.

To help compilation, the instatiation of `Matrix<double, Ordering::RowMajor>`, `Matrix<double, Ordering::ColumnMajor>`, 
`Matrix<std::complex<double>, Ordering::RowMajor>` and `Matrix<std::complex<double>, Ordering::ColumnMajor>` is placed in the file `common_matrix_types.cpp`.

## Matrix Market Reader ##

The code comes with a `template<Ordering O> mm_reader(std::string filename)` function, which reads a file in the *MatrixMarket matrix coordinate real general* 
format and returns a `Matrix<double, O>`. By default, the template argument is `Ordering::RowMajor`