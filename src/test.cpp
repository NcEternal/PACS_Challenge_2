#include "mm_reader.hpp"
#include "chrono.hpp"
#include <string>

int main(int argc, char** argv) {

	using namespace algebra;

	Matrix<double, Ordering::RowMajor> m1(4, 4);
	Matrix<double, Ordering::ColumnMajor> m2(4, 4);

	m1(0, 0) = m2(0, 0) = 5;
	m1(0, 1) = m2(0, 1) = 1;
	m1(1, 0) = m2(1, 0) = 2;
	m1(1, 1) = m2(1, 1) = 0;
	m1(2, 1) = 4;
	m1(0, 3) = m2(2, 2) = 16;
	m1(3, 1) = m2(3, 3) = 14;
	m2(3, 2) = 56;
	m2(0, 3) = 12;
	m2(3, 0) = 17;

	std::cout << "**** ROW ORDERED ****\n\n";
	std::cout << "Uncompressed:\n";
	std::cout << m1;

	std::cout << "Compression:\n";
	m1.compress();

	std::cout << m1;

	std::cout << "Decompression:\n";
	m1.uncompress();

	std::cout << m1;

	std::cout << "Decompression:\n";
	m1.uncompress();

	std::cout << m1;

	std::cout << "Random Access To 0 Element Decompressed:\n";
	std::cout << "(1, 1): " << m1(1, 1) << "\n";

	std::cout << "Random Access To Non 0 Element Decompressed:\n";
	std::cout << "(1, 0): " << m1(1, 0) << "\n";

	std::cout << "Random Access To 0 Element Compressed:\n";
	m1.compress();
	std::cout << "(1, 1): " << static_cast<const Matrix<double, Ordering::RowMajor>>(m1)(1, 1) << "\n";		//operator<< would use the non-const version without the cast and thus
																											//would print an error since m(1, 1) = 0 and 0 elements can't be changed 
																											//while the matrix is compressed
	std::cout << "Random Access To Non 0 Element Compressed:\n";
	std::cout << "(0, 0): " << m1(0, 0) << "\n\n";

	std::cout << "Change Value in Compressed State:\n";
	m1(0, 0) = 10;
	std::cout << m1;

	std::cout << "Resize Matrix (From (4, 4) to (3, 2)):\n";
	m1.resize(3, 2);
	std::cout << m1;
	m1.compress();

	/*--------------*/

	std::cout << "\n\n**** COLUMN ORDERED ****\n\n";
	std::cout << "Uncompressed:\n";
	std::cout << m2;

	std::cout << "Compression:\n";
	m2.compress();

	std::cout << m2;

	std::cout << "Decompression:\n";
	m2.uncompress();

	std::cout << m2;

	std::cout << "Decompression:\n";
	m2.uncompress();

	std::cout << m2;

	std::cout << "Random Access To 0 Element Decompressed:\n";
	std::cout << "(1, 1): " << m2(1, 1) << "\n";

	std::cout << "Random Access To Non 0 Element Decompressed:\n";
	std::cout << "(0, 3): " << m2(0, 3) << "\n";

	std::cout << "Random Access To 0 Element Compressed:\n";
	m2.compress();
	std::cout << "(1, 1): " << static_cast<const Matrix<double, Ordering::ColumnMajor>>(m2)(1, 1) << "\n";	//operator<< would use the non-const version without the cast and thus
																											//would print an error since m(1, 1) = 0 and 0 elements can't be changed 
																											//while the matrix is compressed
	
	std::cout << "Random Access To Non 0 Element Compressed:\n";
	std::cout << "(3, 2): " << m2(3, 2) << "\n\n";

	std::cout << "Change Value in Compressed State:\n";
	m2(0, 0) = 10;
	std::cout << m2;
	
	std::cout << "Resize Matrix (From (4, 4) to (2, 2)):\n";
	m2.resize(2, 2);
	std::cout << m2;
	m2.compress();

	/*-----------*/

	std::cout << "\n\n**** Matrix-Vector multiplication ****\n\n";

	std::vector<double> v(2, 1.);
	auto res1 = m1 * v;
	auto res2 = m2 * v;

	m1.uncompress();
	m2.uncompress();
	
	auto res3 = m1 * v;
	auto res4 = m2 * v;

	std::cout << "Vector:\t";
	for (auto i : v)
		std::cout << i << "\t";
	std::cout << "\n\n";

	std::cout << "M(RowMajor) * V (Compressed):\t";
	for (auto i : res1)
			std::cout << i << "\t";
	std::cout << "\n\n";

	std::cout << "M(ColMajor) * V (Compressed):\t";
	for (auto i : res2)
		std::cout << i << "\t";
	std::cout << "\n\n";
	
	std::cout << "M(RowMajor) * V (Uncompressed):\t";
	for (auto i : res3)
		std::cout << i << "\t";
	std::cout << "\n\n";

	std::cout << "M(RowMajor) * V (Uncompressed):\t";
	for (auto i : res4)
		std::cout << i << "\t";
	std::cout << "\n\n";

	std::vector<std::complex<double>> vi(2);

	vi[0] = std::complex<double>(1, 1);
	vi[1] = std::complex<double>(1, 2);

	std::cout << "Complex Vector:\t";
	for (auto i : vi)
		std::cout << i << "\t";
	std::cout << "\n\n";

	auto resi = m1 * vi;
	std::cout << "M1 * Complex V:\t";
	for (auto i : resi)
		std::cout << i << "\t";
	std::cout << "\n\n";

	Matrix<std::complex<double>, Ordering::RowMajor> mi(2, 2);

	mi(0, 0) = std::complex<double>(1, 1);
	mi(1, 1) = std::complex<double>(2, 3);
	mi(0, 1) = std::complex<double>(1.5, 2.1);
	mi(1, 0) = std::complex<double>(0.5, 0.4);
	
	std::cout << "Complex Matrix:\n" << mi;

	std::cout << "Complex M * Complex V:\t";
	auto resii = mi * vi;
	for (auto i : resii)
		std::cout << i << "\t";
	std::cout << "\n\n";

	std::cout << "Complex M * V:\t";
	auto res5i = mi * v;
	for (auto i : res5i)
		std::cout << i << "\t";
	std::cout << "\n\n";

	/*-----------*/

	std::cout << "\n\n**** Matrix-Matrix multiplication ****\n\n";

	std::cout << "M(RowMajor) * M(ColumnMajor) (Both Uncompressed):\n" << m1 * m2;
	m2.compress();
	std::cout << "M(RowMajor) * M(ColumnMajor) (Row Uncompressed, Column Compressed):\n" << m1 * m2;
	m1.compress();
	std::cout << "M(RowMajor) * M(ColumnMajor) (Both Compressed):\n" << m1 * m2;
	m2.uncompress();
	std::cout << "M(RowMajor) * M(ColumnMajor) (Row Compressed, Column Uncompressed):\n" << m1 * m2;
	std::cout << "M(ColumnMajor) * M(RowMajor) (Both Uncompressed, RHS Complex):\n" << m2 * mi;
	m2.compress();
	std::cout << "M(ColumnMajor) * M(RowMajor) (Column Compressed, Row Uncompressed, RHS Complex):\n" << m2 * mi;
	mi.compress();
	std::cout << "M(ColumnMajor) * M(RowMajor) (Both Compressed, RHS Complex):\n" << m2 * mi;
	m2.uncompress();
	std::cout << "M(ColumnMajor) * M(RowMajor) (Column Uncompressed, Row Compressed, RHS Complex):\n" << m2 * mi;
	mi.uncompress();
	std::cout << "M(RowMajor) * M(ColumnMajor) (Both Uncompressed, LHS Complex):\n" << mi * m2;
	m2.compress();
	std::cout << "M(RowMajor) * M(ColumnMajor) (Row Uncompressed, Column Compressed, LHS Complex):\n" << mi * m2;
	mi.compress();
	std::cout << "M(RowMajor) * M(ColumnMajor) (Both Compressed, LHS Complex):\n" << mi * m2;
	m2.uncompress();
	std::cout << "M(RowMajor) * M(ColumnMajor) (Row Compressed, Column Uncompressed, LHS Complex):\n" << mi * m2;
	mi.uncompress();
	m1.uncompress();
	std::cout << "M(RowMajor) * M(RowMajor) (Both Uncompressed, LHS Complex):\n" << m1 * mi;
	mi.compress();
	std::cout << "M(RowMajor) * M(RowMajor) (LHS Uncompressed, RHS Compressed, LHS Complex):\n" << m1 * mi;
	m1.compress();
	std::cout << "M(RowMajor) * M(RowMajor) (Both Compressed, RHS Complex):\n" << m1 * mi;
	mi.uncompress();
	std::cout << "M(RowMajor) * M(RowMajor) (LHS Compressed, RHS Uncompressed, LHS Complex):\n" << m1 * mi;

	std::cout << "Size Compatibility Check: Mi(2, 2) with M1(3, 2): " << mi * m1;

	m1.uncompress();
	m2.uncompress();
	mi.uncompress();


	/*-------------------*/

	std::cout << "\n**** Norms ****\n\n";

	m1.compress();
	m2.compress();
	mi.compress();

	std::cout << "RowMajor Compressed:\nOne: " << m1.template norm<Norm::One>() << "\nInfinity: " << m1.template norm<Norm::Infinity>() << "\nFrobenius: " << m1.norm();
	std::cout << "\n\nColumnMajor Compressed:\nOne: " << m2.template norm<Norm::One>() << "\nInfinity: " << m2.template norm<Norm::Infinity>() << "\nFrobenius: " << m2.norm();
	std::cout << "\n\nComplex RowMajor Compressed:\nOne: " << mi.template norm<Norm::One>() << "\nInfinity: " << mi.template norm<Norm::Infinity>() << "\nFrobenius: " << mi.norm();

	m1.uncompress();
	m2.uncompress();
	mi.uncompress();

	std::cout << "\n\nRowMajor Uncompressed:\nOne: " << m1.template norm<Norm::One>() << "\nInfinity: " << m1.template norm<Norm::Infinity>() << "\nFrobenius: " << m1.norm();
	std::cout << "\n\nColumnMajor Uncompressed:\nOne: " << m2.template norm<Norm::One>() << "\nInfinity: " << m2.template norm<Norm::Infinity>() << "\nFrobenius: " << m2.norm();
	std::cout << "\n\nComplex RowMajor Uncompressed:\nOne: " << mi.template norm<Norm::One>() << "\nInfinity: " << mi.template norm<Norm::Infinity>() << "\nFrobenius: " << mi.norm();


	/*-------------------*/

	std::cout << "\n\n\n**** Timing ****\n\n";

	std::string filename("lnsp_131.mtx");	
	auto test_un = mm_reader(filename);
	auto test_un_col = mm_reader<Ordering::ColumnMajor>(filename);
	auto test_co = test_un;
	auto test_co_col = test_un_col;
	test_co.compress();	//For some reason the operation following a Matrix::compress gets slowed down, so I keep a compressed matrix ready to better compare 
						//the matrix-vector multiplication
	test_co_col.compress();
	std::vector<double> v_test(test_un.cols(), 1);
	std::vector<std::complex<double>> vi_test(test_un.cols(), std::complex<double>(1, 1));
	Timings::Chrono watch;

	std::cout << "Matrix Sizes:\nRows: " << test_un.rows() << "\nColumns: " << test_un.cols();

	watch.start();
	auto wow = test_un * v_test;
	watch.stop();

	std::cout << "\n\nFirst Test: Matrix * Vector (Uncompressed)\n" << watch;

	watch.start();
	auto wow2 = test_co * v_test;
	watch.stop();

	std::cout << "\n\nSecond Test: Matrix * Vector (Compressed)\n" << watch;

	watch.start();
	auto wow3 = test_un * vi_test;
	watch.stop();

	std::cout << "\n\nThird Test: Matrix * Complex Vector (Uncompressed)\n" << watch;

	watch.start();
	auto wow4 = test_co * vi_test;
	watch.stop();

	std::cout << "\n\nFourth Test: Matrix * Complex Vector (Compressed)\n" << watch;

	std::cout << "\n\n" << std::boolalpha << "Check Result Equality: " << (wow == wow2) << " " << (wow3 == wow4) << "\n";

	watch.start();
	auto wow5 = test_un * test_co;
	watch.stop();

	std::cout << "\n\nFifth Test: Matrix * Matrix (Uncompressed, Compressed, Both RowMajor)\n" << watch;

	watch.start();
	auto wow6 = test_co * test_un;
	watch.stop();

	std::cout << "\n\nSixth Test: Matrix * Matrix (Compressed, Uncompressed, Both RowMajor)\n" << watch;

	watch.start();
	auto wow7 = test_un * test_un;
	watch.stop();

	std::cout << "\n\nSeventh Test: Matrix * Matrix (Both Uncompressed, Both RowMajor)\n" << watch;

	watch.start();
	auto wow8 = test_co * test_co;
	watch.stop();

	std::cout << "\n\nEighth Test: Matrix * Matrix (Both Compressed, Both RowMajor)\n" << watch;

	watch.start();
	auto wow9 = test_un_col * test_co_col;
	watch.stop();

	std::cout << "\n\nNinth Test: Matrix * Matrix (Uncompressed, Compressed, Both ColumnMajor)\n" << watch;

	watch.start();
	auto wow10 = test_co_col * test_un_col;
	watch.stop();

	std::cout << "\n\nTenth Test: Matrix * Matrix (Compressed, Uncompressed, Both ColumnMajor)\n" << watch;

	watch.start();
	auto wow11 = test_un_col * test_un_col;
	watch.stop();

	std::cout << "\n\nEleventh Test: Matrix * Matrix (Both Uncompressed, Both ColumnMajor)\n" << watch;

	watch.start();
	auto wow12 = test_co_col * test_co_col;
	watch.stop();

	std::cout << "\n\nTwelfth Test: Matrix * Matrix (Both Compressed, Both ColumnMajor)\n" << watch;

	watch.start();
	auto wow13 = test_un * test_co_col;
	watch.stop();

	std::cout << "\n\nThirteenth Test: Matrix * Matrix (Uncompressed RowMajor, Compressed ColumnMajor)\n" << watch;

	watch.start();
	auto wow14 = test_co * test_un_col;
	watch.stop();

	std::cout << "\n\nFourteenth Test: Matrix * Matrix (Compressed RowMajor, Uncompressed ColumnMajor)\n" << watch;

	watch.start();
	auto wow15 = test_un * test_un_col;
	watch.stop();

	std::cout << "\n\nFifteenth Test: Matrix * Matrix (Both Uncompressed, LHS RowMajor, RHS ColumnMajor)\n" << watch;

	watch.start();
	auto wow16 = test_co * test_co_col;
	watch.stop();

	std::cout << "\n\nSixteenth Test: Matrix * Matrix (Both Compressed, LHS RowMajor, RHS ColumnMajor)\n" << watch;

	watch.start();
	auto wow17 = test_un_col * test_co;
	watch.stop();

	std::cout << "\n\nSeventeenth Test: Matrix * Matrix (Uncompressed ColumnMajor, Compressed RowMajor)\n" << watch;

	watch.start();
	auto wow18 = test_co_col * test_un;
	watch.stop();

	std::cout << "\n\nEighteenth Test: Matrix * Matrix (Compressed ColumnMajor, Uncompressed RowMajor)\n" << watch;

	watch.start();
	auto wow19 = test_un_col * test_un;
	watch.stop();

	std::cout << "\n\nNineteenth Test: Matrix * Matrix (Both Uncompressed, LHS ColumnMajor, RHS RowMajor)\n" << watch;

	watch.start();
	auto wow20 = test_co_col * test_co;
	watch.stop();

	std::cout << "\n\nTwentieth Test: Matrix * Matrix (Both Compressed, LHS ColumnMajor, RHS RowMajor)\n" << watch;

	std::cout << "\n\n" << std::boolalpha << "Check Result Equality: " << ((wow5 == wow6) && (wow6 == wow7) && (wow7 == wow8) && (wow13 == wow14) && (wow14 == wow15) && (wow15 == wow16) && (wow5 == wow16)) << " "
																	   << ((wow9 == wow10) && (wow10 == wow11) && (wow11 == wow12) && (wow17 == wow18) && (wow18 == wow19) && (wow19 == wow20) && (wow9 == wow20)) << "\n";

	std::cout << "\n\n\n**** Norm Timings ****\n";

	std::cout << "\nNorm One Uncompressed:\n";

	watch.start();
	double norm1_un = test_un.template norm<Norm::One>();
	watch.stop();

	std::cout  << watch;
	std::cout << "\nNorm Infinity Uncompressed:\n";

	watch.start();
	double norm2_un = test_un.template norm<Norm::Infinity>();
	watch.stop();

	std::cout << watch;
	std::cout << "\nNorm Frobenius Uncompressed:\n";

	watch.start();
	double norm3_un = test_un.template norm<Norm::Frobenius>();
	watch.stop();

	std::cout << watch;
	std::cout << "\nNorm One Compressed:\n";

	watch.start();
	double norm1_co = test_co.template norm<Norm::One>();
	watch.stop();

	std::cout << watch;
	std::cout << "\nNorm Infinity Compressed:\n";

	watch.start();
	double norm2_co = test_co.template norm<Norm::Infinity>();
	watch.stop();

	std::cout << watch;
	std::cout << "\nNorm Frobenius Compressed:\n";

	watch.start();
	double norm3_co = test_co.template norm<Norm::Frobenius>();
	watch.stop();

	std::cout << watch;

	std::cout << "\n\n" << std::boolalpha << "Check Equalities: " << (norm1_un == norm1_co) << " " << (norm2_un == norm2_co) << " " << (norm3_un == norm3_co) << "\n";


	return 0;

}