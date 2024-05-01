#include "mm_reader.hpp"
#include <string>
#include <fstream>
#include <sstream>

template <algebra::Ordering OrderType>
algebra::Matrix<double, OrderType> mm_reader(const std::string& filename) {

	using namespace algebra;

	std::ifstream infile(filename);

	/*Check if file exists*/
	if (!infile) {
		std::cerr << "ERROR: Matrix file " << filename << " not found" << std::endl;
		infile.close();
		return Matrix<double, OrderType>(0, 0);
	}

	using namespace std::string_literals;

	std::string line;
	std::getline(infile, line);

	/*Check if file is in any Matrix Market Format*/
	if (line.find("%%MatrixMarket"s) == std::string::npos) {
		std::cerr << "ERROR: File " << filename << " does not follow the MatrixMarket format" << std::endl;
		infile.close();
		return Matrix<double, OrderType>(0, 0);
	}

	/*Check if file is not in the only Matrix Market format currently supported*/
	if (line.compare("%%MatrixMarket matrix coordinate real general"s)) {
		std::cerr << "ERROR: Matrix Type not supported" << std::endl;
		infile.close();
		return Matrix<double, OrderType>(0, 0);
	}

	/*Ignore all commented lines*/
	std::getline(infile, line);
	while (line.at(0) == '%')
		std::getline(infile, line);

	/*Read Matrix sizes*/
	std::size_t rows, cols;

	std::stringstream line_s(line);

	line_s >> rows >> cols;

	/*Read Matrix elements*/
	Matrix<double, OrderType> result(rows, cols);

	while (std::getline(infile, line)) {

		std::stringstream element(line);

		std::size_t i, j;
		double val;

		element >> i >> j >> val;

		result(i - 1, j - 1) = val;

	}

	infile.close();

	return result;

}

template algebra::Matrix<double, algebra::Ordering::RowMajor> mm_reader(const std::string& filename);
template algebra::Matrix<double, algebra::Ordering::ColumnMajor> mm_reader(const std::string& filename);