#include "Matrix.hpp"

namespace algebra {

	template class Matrix<double, Ordering::RowMajor>;
	template class Matrix<double, Ordering::ColumnMajor>;
	template class Matrix<std::complex<double>, Ordering::RowMajor>;
	template class Matrix<std::complex<double>, Ordering::ColumnMajor>;

}