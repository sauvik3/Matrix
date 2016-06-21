#include <vector>
#include <functional>
#include <string>
#include <iostream>
#include <stdexcept>

#include "matrix.h"
#include "matrix_exception.h"
#include "matrix_io.h"

/*---------------------------- Prototypes -----------------------------*/
template<typename _Tp> void add_matrix();
template<typename _Tp> void multiply_matrix();
template<typename _Tp> void transpose_matrix();
template<typename _Tp> void adjoint_matrix();
template<typename _Tp, typename _Dt> void inverse_matrix();
template<typename _Tp> void determinant_matrix();
/*---------------------------------------------------------------------*/


/*---------------------------- Main Routine ---------------------------*/
int main() {
	std::vector< std::string > menu = {
		"Add Matrices",
		"Multiply Matrices",
		"Matrix Transpose",
		"Matrix Adjoint",
		"Matrix Inverse",
		"Determinant"
	};

	std::vector< std::function<void()> > menuAction =
		std::vector<std::function<void()> >(0);
	menuAction.push_back(add_matrix<double>);
	menuAction.push_back(multiply_matrix<double>);
	menuAction.push_back(transpose_matrix<double>);
	menuAction.push_back(adjoint_matrix<double>);
	menuAction.push_back(inverse_matrix<double, float>);
	menuAction.push_back(determinant_matrix<float>);

	size_t choice;
	size_t menuCount = menu.size();

	/* Display Menu Endlessly */
	for (;;) {
		size_t i;
		std::cout << "\nMatrix Menu";
		for (i = 1; i <= menuCount; ++i) {
			std::cout << "\n" << i << ". " << menu.at(i-1);
		}
		std::cout << "\n" << i << ". " << "Exit";
		std::cout << "\nEnter choice : ";
		std::cin >> choice;

		if (choice == i) {
			return 0;
		}
		else {
			if (choice >0 && choice <= menuCount) {
				if (menuAction.at(choice - 1) != NULL) {
					try {
						menuAction.at(choice - 1)();
					}
					catch (MatrixException &e) {
						std::cout << "\nError : " << e.what() << std::endl;
					}
					catch (std::runtime_error &e) {
						std::cout << "\nError : " << e.what() << std::endl;
					}
				}
			}
			else {
				std::cout << "\nInvalid Input !!!\n";
			}
		}
	}
}
/*---------------------------- I/O Wrappers ---------------------------*/
template<typename _Tp> void add_matrix() {
	Matrix<_Tp> A, B, C;

	/* Read details for Matrix A */
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	/* Read details for Matrix B */
	std::cout << "\nEnter details for Matrix B :-\n";
	B = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);
	std::cout << "\nMatrix B :-";
	print_matrix<_Tp>(B);

	/* Perform Addition */
	C = A + B;

	/* Print out result */
	std::cout << "\nSum :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp> void multiply_matrix() {
	Matrix<_Tp> A, B, C;

	/* Read details for Matrix A */
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	/* Read details for Matrix B */
	std::cout << "\nEnter details for Matrix B :-\n";
	B = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);
	std::cout << "\nMatrix B :-";
	print_matrix<_Tp>(B);

	/* Perform Multiplication */
	C = A * B;

	/* Print out result */
	std::cout << "\nProduct :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp> void transpose_matrix() {
	Matrix<_Tp> A, C;

	/* Read details for Matrix A */
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);
	
	/* Perform Transpose Operation */
	C = transpose<_Tp>(A);

	/* Print out result */
	std::cout << "\nTranspose :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp> void adjoint_matrix() {
	Matrix<_Tp> A, C;

	/* Read details for Matrix A */
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);

	/* Perform Adjoint Operation */
	C = adjoint<_Tp>(A);

	/* Print out result */
	std::cout << "\nAdjoint :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp, typename _Dt> void inverse_matrix() {
	Matrix<_Tp> A;
	Matrix<_Dt> C;

	/* Read details for Matrix A */
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);

	/* Perform Inverse Operation */
	C = inverse<_Tp, _Dt>(A);

	/* Print out result */
	std::cout << "\nInverse :-";
	print_matrix<_Dt>(C);
}

template<typename _Tp> void determinant_matrix() {
	Matrix<_Tp> A;
	_Tp det;

	/* Read details for Matrix A */
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);

	/* Find Determinant of Matrix */
	det = determinant<_Tp>(A);

	/* Print out result */
	std::cout << "\nDeterminant :-\n";
	std::cout << det << "\n";
}
/*---------------------------------------------------------------------*/