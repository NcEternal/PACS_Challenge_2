#ifndef HH_MATRIX_TRAITS_HH
#define HH_MATRIX_TRAITS_HH

#include <array>
#include <map>
#include <complex>
#include <type_traits>
#include <concepts>

namespace algebra {

/*Norm types*/
enum class Norm { One, Infinity, Frobenius };


/*Ordering types*/
enum class Ordering { RowMajor, ColumnMajor };


/*Definition of the concept of a number*/
template <typename T>
concept Real_Number = std::integral<T> || std::floating_point<T>;

template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

template <typename T>
concept Complex_Number = is_complex<T>::value;

template <typename T>
concept Number = Real_Number<T> || Complex_Number<T>;


/*Forward Declaration so that MatrixTraits can use ordering_comp and ordering_comp can use MatrixTraits::matrix_index*/
template <Ordering O>
struct ordering_comp;

/*Matrix Traits*/
struct MatrixTraits {

	using matrix_index = std::array<std::size_t, 2>;
	
	template <Number N, Ordering O = Ordering::RowMajor>
	using matrix_type = std::map<matrix_index, N, ordering_comp<O>>;

	template <Number N>
	using vector_type = std::map<std::size_t, N>;	// Sparse Vector used to get rows and columns in a matrix

};

/*Comparison operator for the different Ordering types*/
template <Ordering O = Ordering::RowMajor>
struct ordering_comp {
	constexpr bool operator() (const MatrixTraits::matrix_index& a, const MatrixTraits::matrix_index& b) const {
		if constexpr (O == Ordering::RowMajor)
			return a < b;
		else
			return (a[1] < b[1]) || (a[1] == b[1] && a[0] < b[0]);
	};
};

}

#endif
