#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <stdexcept>

#include "matrix.h"
#include "matrix_io.h"

/*--------------------------- I/O Functions ---------------------------*/
template<typename _Tp> Matrix<_Tp> read_matrix() {
	int i, j;
	int m, n;
	std::vector< std::vector <_Tp> > V;

	std::cout << "Enter m : ";
	std::cin >> m;
	std::cout << "Enter n : ";
	std::cin >> n;

	if (m < 1 || n < 1) {
		throw std::runtime_error("Must specify unsigned non-zero value");
	}

	std::cout << "\n";
	for (i = 0; i < m; ++i) {
		std::vector <_Tp> x;
		for (j = 0; j < n; ++j) {
			_Tp v;
			std::cout << "Enter element[" << i << "][" << j << "] : ";
			std::cin >> v;
			x.push_back(v);
		}
		V.push_back(x);
	}
	Matrix<_Tp> A(static_cast<size_t>(m),
		static_cast<size_t>(n), V);

	return A;
}

template<typename _Tp> void print_matrix(const Matrix<_Tp> &A) throw() {
	size_t i, j, k;
	size_t w1, w2;

	if ((A.getM() == 0) || (A.getN() == 0)) {
		std::cout << "\n" << (unsigned char)LUCNR;
		std::cout << "     ";
		std::cout << (unsigned char)RUCNR;
		std::cout << "\n" << (unsigned char)VRT;
		std::cout << "Empty";
		std::cout << (unsigned char)VRT << "\n";
		std::cout << (unsigned char)LLCNR;
		std::cout << "     ";
		std::cout << (unsigned char)RLCNR;
	}
	else {
		w1 = 0;
		for (i = 0; i < A.getM(); ++i) {
			for (j = 0; j < A.getN(); ++j) {
				size_t l;
				std::ostringstream tmpstream;
				tmpstream << std::right << std::fixed << std::setprecision(2)
					<< std::noshowpos << A(i, j);
				l = tmpstream.str().length();
				w1 = (w1 > l ? w1 : l);
			}
		}
		w2 = ((w1 + 1)*A.getN());

		std::cout << "\n" << (unsigned char)LUCNR;
		for (k = 0; k < w2; ++k) {
			std::cout << " ";
		}
		std::cout << (unsigned char)RUCNR;

		for (i = 0; i < A.getM(); ++i) {
			std::cout << "\n";
			for (j = 0; j < A.getN(); ++j) {
				if (j == 0) {
					std::cout << (unsigned char)VRT;
				}
				std::cout << std::right << std::setfill(' ') << std::setw(w1)
					<< std::fixed << std::setprecision(2)
					<< std::noshowpos << A(i, j) << " ";
				if (j == ((A.getN()) - 1)) {
					std::cout << (unsigned char)VRT;
				}
			}
		}

		std::cout << "\n" << (unsigned char)LLCNR;
		for (k = 0; k < w2; ++k) {
			std::cout << " ";
		}
		std::cout << (unsigned char)RLCNR;
	}
}
/*---------------------------------------------------------------------*/

/* Specializations */
template Matrix<double> read_matrix();
template Matrix<float> read_matrix();
template void print_matrix(const Matrix<double>&) throw();
template void print_matrix(const Matrix<float>&) throw();
/* End of Specializations */