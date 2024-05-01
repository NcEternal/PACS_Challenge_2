#ifndef HH_MATRIX_HH
#define HH_MATRIX_HH

#include "Matrix_Traits.hpp"
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>

namespace algebra {

template <Number ValueType, Ordering OrderType = Ordering::RowMajor>
class Matrix {
public:
	/*Traits*/
	using T = MatrixTraits;

private:

	/*Matrix for Uncompressed Format*/
	T::matrix_type<ValueType, OrderType> m_mat;

	/*Vectors for Compressed Format*/
	std::vector<ValueType> m_values;
	std::vector<std::size_t> m_inner_idxs;
	std::vector<std::size_t> m_outer_idxs;

	/*State Variables*/
	bool compressed = false;
	std::size_t m_rows;
	std::size_t m_cols;


	/*Avoid repeating code for cases when the vector and matrix value types differ only because one of them is complex and the other isn't*/
	template<Number Res, Number N>
	std::vector<Res> vector_multiplication_helper(const std::vector<N>& v) const;

	/*Avoid repeating code for cases when the 2 matrices value types differ only because one of them is complex and the other isn't*/
	template<Number Res, Number N, Ordering O>
	Matrix<Res, OrderType> matrix_multiplication_helper(const Matrix<N, O>& rhs) const;

public:
	/*Constructors*/
	Matrix() = default;
	Matrix(std::size_t r, std::size_t c) : m_rows(r), m_cols(c) {}

	/*Access*/
	inline ValueType& operator() (std::size_t i, std::size_t j);		//NEEDS CHECKING
	inline ValueType operator() (std::size_t i, std::size_t j) const;

	/*Compression and Decompression*/
	void compress();
	void uncompress();

	/*Resize*/ 
	void resize(std::size_t new_r, std::size_t new_c);	// Uncompresses the Matrix in order to erase elements that go out of the new bounds

	/*Norm*/
	template<Norm NormType = Norm::Frobenius>
	double norm() const;

	/*Getters*/
	bool is_compressed() const { return compressed; }
	std::size_t rows() const { return m_rows; }
	std::size_t cols() const { return m_cols; }
	T::vector_type<ValueType> get_row(std::size_t i) const;	// Returns row as a sparse vector (i.e. map<size_t, ValueType>)
	T::vector_type<ValueType> get_col(std::size_t j) const; // Returns column as a sparse vector (i.e. map<size_t, ValueType>)



	/*Matrix-Vector Multiplication*/
	template <Number U, Ordering O>																					// Same Type
	friend std::vector<U> operator*(const Matrix<U, O>& m, const std::vector<U>& v);								

	template <Number U, Ordering O>																					// Same Base Type
	friend std::vector<std::complex<U>> operator*(const Matrix<U, O>& m, const std::vector<std::complex<U>>& v);	// But Vector is Complex

	template <Number U, Ordering O>																					// Same Base Type
	friend std::vector<std::complex<U>> operator*(const Matrix<std::complex<U>, O>& m, const std::vector<U>& v);	// But Matrix is Complex

	/*Matrix-Matrix Multiplication*/
	template <Number U, Ordering O1, Ordering O2>																	// Same Type
	friend Matrix<U, O1> operator*(const Matrix<U, O1>& lhs, const Matrix<U, O2>& rhs);

	template <Number U, Ordering O1, Ordering O2>																	// Same Base Type
	friend Matrix<std::complex<U>, O1> operator*(const Matrix<U, O1>& lhs, const Matrix<std::complex<U>, O2>& rhs);	// But Right-hand-side is Complex

	template <Number U, Ordering O1, Ordering O2>																	// Same Base Type
	friend Matrix<std::complex<U>, O1> operator*(const Matrix<std::complex<U>, O1>& lhs, const Matrix<U, O2>& rhs);	// But Matrix is Complex

	/*Printing*/
	template <Number U, Ordering O>
	friend std::ostream& operator<<(std::ostream& os, const Matrix<U, O>& m);

	/*Rudimentary Comparison Operator*/
	template <Number U, Ordering O>
	friend bool operator==(const Matrix<U, O>& lhs, const Matrix<U, O>& rhs);

};

	/*Helper Function for Matrix-Matrix Multiplication*/
	namespace matrix_matrix_multiplication_helper {
		template <Number ResType, Number N1, Number N2>
		ResType row_times_column(const MatrixTraits::vector_type<N1>& row, const MatrixTraits::vector_type<N2>& column);
	}







/*-------------------------------------*/
/*** ------- IMPLEMENTATIONS ------- ***/
/*-------------------------------------*/



/*** NON-CONST ACCESS ***/	//NEEDS CHECKING
template <Number ValueType, Ordering OrderType>
ValueType& Matrix<ValueType, OrderType>::operator() (std::size_t i, std::size_t j) {
	
	/*Check Index Bounds*/
	if (i >= m_rows || j >= m_cols) {
		std::cerr << "ERROR: Index out of bounds" << std::endl;
		return Matrix(1, 1)(0, 0);
	}
	
	/*Uncompressed Case*/
	if (!compressed) {
		T::matrix_index idxs = { i, j };
		return m_mat[idxs];
	}

	/*Deal with both RowMajor and ColumnMajor matrices with the same code by swapping the values of the indeces*/
	if constexpr (OrderType == Ordering::ColumnMajor) {
		std::swap(i, j);
	}
	/*Access for compressed matrix. Each row/column (depending on OrderType) within m_outer_idxs is sorted, so we can use std::lower_bound for efficiency*/
	auto start = m_outer_idxs.cbegin();
	auto search_start = start + m_inner_idxs[i], search_end = start + m_inner_idxs[i + 1];
	long int el_n = std::lower_bound(search_start, search_end, j) - start;
	if (el_n != search_end - start && m_outer_idxs[el_n] == j)
		return m_values[el_n];

	/*Error if trying to access a 0 element in compressed state*/
	std::cerr << "ERROR: Cannot change value of a 0 element in compressed state" << std::endl;
	return Matrix(1, 1)(0, 0);

}



/*** CONST ACCESS ***/
template <Number ValueType, Ordering OrderType>
ValueType Matrix<ValueType, OrderType>::operator() (std::size_t i, std::size_t j) const {
	
	/*Check Index Bounds*/
	if (i >= m_rows || j >= m_cols) {
		std::cerr << "ERROR: Index out of bounds" << std::endl;
		return std::numeric_limits<ValueType>::quiet_NaN();
	}

	/*Uncompressed Case*/
	if (!compressed) {
		T::matrix_index idxs = { i, j };
		auto el = m_mat.find(idxs);
		if (el == m_mat.cend())
			return static_cast<ValueType>(0.);
		return el->second;
	}

	/*Deal with both RowMajor and ColumnMajor matrices with the same code by swapping the values of the indeces*/
	if constexpr (OrderType == Ordering::ColumnMajor) {
		std::swap(i, j);
	}
	/*Access for compressed matrix. Each row/column (depending on OrderType) within m_outer_idxs is sorted, so we can use std::lower_bound for efficiency*/
	auto start = m_outer_idxs.cbegin();
	auto search_start = start + m_inner_idxs[i], search_end = start + m_inner_idxs[i + 1];
	long int el_n = std::lower_bound(search_start, search_end, j) - start;
	if (el_n != search_end - start && m_outer_idxs[el_n] == j)
		return m_values[el_n];
	return static_cast<ValueType>(0.);

}



/*** GET ROW ***/
template <Number ValueType, Ordering OrderType>
MatrixTraits::vector_type<ValueType> Matrix<ValueType, OrderType>::get_row(std::size_t i) const {
	
	T::vector_type<ValueType> row;

	/*RowMajor case*/
	if constexpr (OrderType == Ordering::RowMajor) {
		if (!compressed) {
			T::matrix_index start = { i, 0 }, end = { i, m_cols };
			auto start_it = m_mat.lower_bound(start), end_it = m_mat.upper_bound(end);
			for (auto it = start_it; it != end_it; ++it)
				row[it->first[1]] = it->second;
		}
		else {
			for (std::size_t el_n = m_inner_idxs[i]; el_n < m_inner_idxs[i + 1]; ++el_n) {
				row[m_outer_idxs[el_n]] = m_values[el_n];
			}
		}
	}
	/*ColumnMajor case*/
	else {
		for (std::size_t j = 0; j < m_cols; ++j) {
			ValueType a_ij = this->operator()(i, j);	//Calls the const version, always within bounds
			if (a_ij != static_cast<ValueType>(0.)) {
				row[j] = a_ij;
			}
		}
	}

	return row;

}



/*** GET COLUMN ***/
template <Number ValueType, Ordering OrderType>
MatrixTraits::vector_type<ValueType> Matrix<ValueType, OrderType>::get_col(std::size_t j) const {
	
	T::vector_type<ValueType> column;

	/*ColumnMajor case*/
	if constexpr (OrderType == Ordering::ColumnMajor) {
		if (!compressed) {
			T::matrix_index start = { 0, j }, end = { m_rows, j };
			auto start_it = m_mat.lower_bound(start), end_it = m_mat.upper_bound(end);
			for (auto it = start_it; it != end_it; ++it)
				column[it->first[0]] = it->second;
		}
		else {
			for (std::size_t el_n = m_inner_idxs[j]; el_n < m_inner_idxs[j + 1]; ++el_n) {
				column[m_outer_idxs[el_n]] = m_values[el_n];
			}
		}
	}
	/*RowMajor case*/
	else {
		for (std::size_t i = 0; i < m_rows; ++i) {
			ValueType a_ij = this->operator()(i, j);	//Calls the const version, always within bounds
			if (a_ij != static_cast<ValueType>(0.)) {
				column[i] = a_ij;
			}
		}
	}

	return column;

}



/*** COMPRESSION ***/
template <Number ValueType, Ordering OrderType>
void Matrix<ValueType, OrderType>::compress() {
	
	/*If matrix is already compressed we don't need to do anything*/
	if (compressed)
		return;

	/*Get size of the major dimension and position of the outer index*/
	std::size_t major_size = m_rows;
	std::size_t outer_idx = 1;
	if constexpr (OrderType == Ordering::ColumnMajor) {
		major_size = m_cols;
		outer_idx = 0;
	}

	/*Lambdas to get lower and upper bounds in both orderings*/
	auto low_bound = [this](std::size_t major_idx) {if constexpr (OrderType == Ordering::RowMajor) return m_mat.lower_bound({ major_idx, 0 }); else return m_mat.lower_bound({ 0, major_idx }); };
	auto upp_bound = [this](std::size_t major_idx) {if constexpr (OrderType == Ordering::RowMajor) return m_mat.upper_bound({ major_idx, m_cols }); else return m_mat.upper_bound({ m_rows, major_idx }); };

	/*Resize inner index vector*/
	m_inner_idxs.resize(major_size + 1);

	std::size_t nnz = 0;

	for (std::size_t i = 0; i < major_size; ++i) {

		m_inner_idxs[i] = nnz;

		/* Get Row/Column */
		auto low = low_bound(i), high = upp_bound(i);

		/* Loop over row/column to get values and outer indices */
		for (auto j = low; j != high; ++j) {
			
			if (j->second == static_cast<ValueType>(0.)) 	//Ignore 0s that were added in uncompressed mode
				continue;
			
			m_values.push_back(j->second);
			
			m_outer_idxs.push_back(j->first[outer_idx]);
			
			++nnz;
		
		}

	}

	m_inner_idxs[major_size] = nnz;

	/*Change State and clear*/
	compressed = true;
	m_mat.clear();

}



/*** DECOMPRESSION ***/
template <Number ValueType, Ordering OrderType>
void Matrix<ValueType, OrderType>::uncompress() {

	/*If matrix is already uncompressed we don't need to do anything*/
	if (!compressed)
		return;

	/*Get positions of inner and outer indeces*/
	std::size_t major_idx = 0, minor_idx = 1;
	if constexpr (OrderType == Ordering::ColumnMajor) {
		major_idx = 1;
		minor_idx = 0;
	}

	for (std::size_t i = 0; i < m_inner_idxs.size() - 1; ++i) {
		for (std::size_t j = m_inner_idxs[i]; j < m_inner_idxs[i + 1]; ++j) {
			
			if (m_values[j] == static_cast<ValueType>(0.))			//Ignore elements that were changed to 0
				continue;
			
			T::matrix_index idx;
			idx[major_idx] = i;
			idx[minor_idx] = m_outer_idxs[j];

			m_mat[idx] = m_values[j];
		}
	}

	/*Change State and Clear*/
	compressed = false;
	m_values.clear();
	m_inner_idxs.clear();
	m_outer_idxs.clear();

}



/*** RESIZE ***/
template <Number ValueType, Ordering OrderType>
void Matrix<ValueType, OrderType>::resize(std::size_t new_r, std::size_t new_c) {
	
	/*Uncompress to handle elements that go outside the new bounds*/
	uncompress();

	/*If the matrix only gets larger than there's no need to do anything other than update the bounds*/
	if (new_r >= m_rows && new_c >= m_cols) {
		m_rows = new_r;
		m_cols = new_c;
		return;
	}

	/*If the matrix gets less rows than it had, delete elements out of bounds.
	For RowMajor matrices we can just delete all elements starting from the first element of the first out of bounds row
	For ColumnMajor matrices we need to check all elements one by one*/
	if (new_r < m_rows) {
		if constexpr (OrderType == Ordering::RowMajor) {
			T::matrix_index erase_start = { new_r, 0 };
			auto er_it = m_mat.lower_bound(erase_start);
			m_mat.erase(er_it, m_mat.end());
		}
		else {
			auto it = m_mat.begin();
			while (it != m_mat.end()) {
				if (it->first[0] >= new_r)
					it = m_mat.erase(it);
				else
					++it;
			}
		}
	}

	m_rows = new_r;

	/*Same as what happens in the row case, but with ColumnMajor and RowMajor inverted*/
	if (new_c < m_cols) {
		if constexpr (OrderType == Ordering::ColumnMajor) {
			T::matrix_index erase_start = { 0, new_c };
			auto er_it = m_mat.lower_bound(erase_start);
			m_mat.erase(er_it, m_mat.end());
		}
		else {
			auto it = m_mat.begin();
			while (it != m_mat.end()) {
				if (it->first[1] >= new_c)
					it = m_mat.erase(it);
				else
					++it;
			}
		}
	}

	m_cols = new_c;

}



/*** NORM ***/
template <Number ValueType, Ordering OrderType>
template<Norm NormType>
double Matrix<ValueType, OrderType>::norm() const {

	double res = 0.;

	/*Frobenius Norm*/
	if constexpr (NormType == Norm::Frobenius) {
		

		if (!compressed) {
			for (auto it = m_mat.cbegin(); it != m_mat.cend(); ++it) {
				double abs_value = std::abs(it->second);
				res += abs_value * abs_value;
			}
		}
		else {
			for (std::size_t el = 0; el < m_values.size(); ++el){
				double abs_value = std::abs(m_values[el]);
				res += abs_value * abs_value; 
			}
		}

		res = std::sqrt(res);

	}
	/*One and Infinity Norms*/
	else {

		/*Variables to differentiate each case*/
		std::size_t idx_of_interest = 0;
		std::size_t size_of_interest = m_rows;
		if constexpr (NormType == Norm::One) {
			idx_of_interest = 1;
			size_of_interest = m_cols;
		}

		/*In these cases the way the matrix is stored makes things easier to write with a single for loop.*/
		constexpr bool special_case = (NormType == Norm::One && OrderType == Ordering::RowMajor) || (NormType == Norm::Infinity && OrderType == Ordering::ColumnMajor);

		/*Generator for column/row index*/
		std::size_t idx = 0;
		auto generator = [&idx, this](const std::size_t el_n) {if constexpr (special_case) return m_outer_idxs[el_n]; else if (el_n >= m_inner_idxs[idx + 1]) ++idx;  return idx; };

		/*Store partial results. We could avoid using it when special_case == false, by using 2 loops, with the outer one looping over rows/columns 
		while the inner one computes the sum per row/column and then having res = std::max(res, sum) at the end of each inner loop, but it would introduce more branches*/
		std::vector<double> temp(size_of_interest, 0.);

		if (!compressed) {
			for (auto it = m_mat.cbegin(); it != m_mat.cend(); ++it)
				temp[it->first[idx_of_interest]] += std::abs(it->second);
		}
		else {
			for (std::size_t el_n = 0; el_n < m_values.size(); ++el_n)
				temp[generator(el_n)] += std::abs(m_values[el_n]);
		}

		res = *std::max_element(temp.cbegin(), temp.cend());

	}

	return res;

}



/*** MATRIX VECTOR MULTIPLICATION ***/
/*Same Type*/
template <Number ValueType, Ordering OrderType>
std::vector<ValueType> operator*(const Matrix<ValueType, OrderType>& m, const std::vector<ValueType>& v) {
	return m.template vector_multiplication_helper<ValueType>(v);
}

/*Same Base Type, but vector is Complex*/
template <Number ValueType, Ordering OrderType>
std::vector<std::complex<ValueType>> operator*(const Matrix<ValueType, OrderType>& m, const std::vector<std::complex<ValueType>>& v) {
	return m.template vector_multiplication_helper<std::complex<ValueType>>(v);
}

/*Same Base Type, but matrix is Complex*/
template <Number ValueType, Ordering OrderType>
std::vector<std::complex<ValueType>> operator*(const Matrix<std::complex<ValueType>, OrderType>& m, const std::vector<ValueType>& v) {
	return m.template vector_multiplication_helper<std::complex<ValueType>>(v);
}

/*Helper to handle all 3 cases*/
template <Number ValueType, Ordering OrderType>
template <Number Res, Number N>
std::vector<Res> Matrix<ValueType, OrderType>::vector_multiplication_helper(const std::vector<N>& v) const {

	/*Size Compatibility check*/
	if (m_cols != v.size()) {
		std::cerr << "ERROR: Incompatible sizes" << std::endl;
		return std::vector<Res>();
	}

	std::vector<Res> result(m_rows, 0.);

	/*Uncompressed case*/
	if (!compressed) {
		for (auto i = m_mat.cbegin(); i != m_mat.cend(); ++i) {
			result[i->first[0]] += i->second * v[i->first[1]];
		}
		return result;
	}

	/*Compressed Case, Row Major*/
	if constexpr (OrderType == Ordering::RowMajor) {
		for (std::size_t i = 0; i < m_rows; ++i)
			for (std::size_t el_n = m_inner_idxs[i]; el_n < m_inner_idxs[i + 1]; ++el_n)
				result[i] += m_values[el_n] * v[m_outer_idxs[el_n]];
		return result;
	}
	/*Compressed Case Column Major*/
	else {
		for (std::size_t j = 0; j < m_cols; ++j)
			for (std::size_t el_n = m_inner_idxs[j]; el_n < m_inner_idxs[j + 1]; ++el_n)
				result[m_outer_idxs[el_n]] += m_values[el_n] * v[j];
		return result;
	}

}



/*** MATRIX MATRIX MULTIPLICATION ***/
/*Same Type*/
template <Number ValueType, Ordering O1, Ordering O2>																	
Matrix<ValueType, O1> operator*(const Matrix<ValueType, O1>& lhs, const Matrix<ValueType, O2>& rhs) {
	return lhs.template matrix_multiplication_helper<ValueType>(rhs);
}

/*Same Base Type, but RHS is Complex*/
template <Number ValueType, Ordering O1, Ordering O2>																	
Matrix<std::complex<ValueType>, O1> operator*(const Matrix<ValueType, O1>& lhs, const Matrix<std::complex<ValueType>, O2>& rhs) {
	return lhs.template matrix_multiplication_helper<std::complex<ValueType>>(rhs);
}

/*Same Base Type, but LHS is Complex*/
template <Number ValueType, Ordering O1, Ordering O2>
Matrix<std::complex<ValueType>, O1> operator*(const Matrix<std::complex<ValueType>, O1>& lhs, const Matrix<ValueType, O2>& rhs) {
	return lhs.template matrix_multiplication_helper<std::complex<ValueType>>(rhs);
}

/*Helper to handle all 3 cases*/
template<Number ValueType, Ordering OrderType>
template<Number Res, Number N, Ordering O>
Matrix<Res, OrderType> Matrix<ValueType, OrderType>::matrix_multiplication_helper(const Matrix<N, O>& rhs) const {
	
	if (m_cols != rhs.rows()) {
		std::cerr << "ERROR: Incompatible sizes" << std::endl;
		return Matrix<Res, OrderType>();
	}

	Matrix<Res, OrderType> result(m_rows, rhs.cols());

	/*Try to have the more efficient calculation in the inner loop / avoid repeating the inner loop when lhs has 1 row or rhs has 1 column*/
	if ((OrderType == Ordering::RowMajor && m_rows != 1) || rhs.cols() == 1) {
		for (std::size_t j = 0; j < rhs.cols(); ++j) {
			auto col = rhs.get_col(j);
			for (std::size_t i = 0; i < m_rows; ++i) {
				auto row = get_row(i);
				Res row_times_col = matrix_matrix_multiplication_helper::row_times_column<Res>(row, col);
				if (row_times_col != static_cast<Res>(0.))
					result(i, j) = row_times_col;
			}
		}
	}
	else {
		for (std::size_t i = 0; i < m_rows; ++i) {
			auto row = get_row(i);
			for (std::size_t j = 0; j < rhs.cols(); ++j) {
				auto col = rhs.get_col(j);
				Res row_times_col = matrix_matrix_multiplication_helper::row_times_column<Res>(row, col);
				if (row_times_col != static_cast<Res>(0.))
					result(i, j) = row_times_col;
			}
		}
	}

	return result;

}

	/*Row-Column Multiplication*/
	namespace matrix_matrix_multiplication_helper {
		template <Number ResType, Number N1, Number N2>
		ResType row_times_column(const MatrixTraits::vector_type<N1>& row, const MatrixTraits::vector_type<N2>& column) {

			auto col_it = column.cbegin(), col_end = column.cend();
			ResType res(0.);

			for (auto row_it = row.cbegin(); row_it != row.cend() && col_it != col_end; ++row_it) {
				while (col_it != col_end && col_it->first < row_it->first)
					++col_it;
				if (col_it != col_end && col_it->first == row_it->first)
					res += row_it->second * col_it->second;
			}

			return res;

		}
	}



/*** PRINTER ***/
template <Number ValueType, Ordering OrderType>
std::ostream& operator<<(std::ostream& os, const Matrix<ValueType, OrderType>& m) {
	
	/*Uncompressed Printing*/
	if (!m.is_compressed()) {
		os << "i\tj\tval\n";
		for (auto i = m.m_mat.cbegin(); i != m.m_mat.cend(); ++i)
			os << i->first[0] << "\t" << i->first[1] << "\t" << i->second << "\n";
	}
	/*Compressed Printing, Distinguishes Ordering*/
	else {

		os << "Values:\t\t\t\t";
		for (const auto& i : m.m_values)
			os << i << "\t";

		if constexpr (OrderType == Ordering::RowMajor) {
			os << "\nColumn Indeces:\t\t\t";
			for (const auto& i : m.m_outer_idxs)
				os << i << "\t";
			os << "\nFirst Element of Each Row:\t";
			for (const auto& i : m.m_inner_idxs)
				os << i << "\t";
		}
		else {
			os << "\nRow Indeces:\t\t\t";
			for (const auto& i : m.m_outer_idxs)
				os << i << "\t";
			os << "\nFirst Element of Each Column:\t";
			for (const auto& i : m.m_inner_idxs)
				os << i << "\t";
		}

		os << "\n";
	}

	os << std::endl;
	return os;
}



template <Number U, Ordering O>
bool operator==(const Matrix<U, O>& lhs, const Matrix<U, O>& rhs) {
	if (lhs.is_compressed() != rhs.is_compressed()) {
		std::cerr << "ERROR: Cannot compare matrices in different compression state\n";
		return false;
	}
	if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
		return false;
	if (lhs.is_compressed())
		return lhs.m_inner_idxs == rhs.m_inner_idxs && lhs.m_outer_idxs == rhs.m_outer_idxs && lhs.m_values == rhs.m_values;
	return lhs.m_mat == rhs.m_mat;
}



extern template class Matrix<double, Ordering::RowMajor>;
extern template class Matrix<double, Ordering::ColumnMajor>;
extern template class Matrix<std::complex<double>, Ordering::RowMajor>;
extern template class Matrix<std::complex<double>, Ordering::ColumnMajor>;

}
#endif