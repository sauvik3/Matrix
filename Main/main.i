template<typename _Tp> class __declspec(dllimport) Matrix;

template<typename _Tp> __declspec(dllimport) _Tp __stdcall determinant(const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllimport) Matrix<_Tp> __stdcall transpose(const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllimport) Matrix<_Tp> __stdcall cofactor(const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllimport) Matrix<_Tp> __stdcall adjoint(const Matrix<_Tp>&);
template<typename _Tp, typename _Dt> __declspec(dllimport) Matrix<_Dt> __stdcall inverse(const Matrix<_Tp>&);

template<typename _Tp, typename _Dt> __declspec(dllimport) Matrix<_Dt> __stdcall operator*(const Matrix<_Tp>&, const _Dt&);
template<typename _Tp, typename _Dt> __declspec(dllimport) Matrix<_Dt> __stdcall operator/(const Matrix<_Tp>&, const _Dt&);

template<typename _Tp> __declspec(dllimport) bool __stdcall operator==(const Matrix<_Tp>&, const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllimport) bool __stdcall operator!=(const Matrix<_Tp>&, const Matrix<_Tp>&);
template<typename _Tp> bool operator>(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator<(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator>=(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator<=(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;

template<typename _Tp>
class __declspec(dllimport) Matrix {
	std::size_t m;
	std::size_t n;
	std::vector< std::vector <_Tp> > data;
public:
	std::size_t __stdcall getM() const throw();
	std::size_t __stdcall getN() const throw();
	std::vector< std::vector< _Tp > > __stdcall getData() const throw();

	__stdcall Matrix();
	__stdcall Matrix(const std::size_t &, const std::size_t &);
	__stdcall Matrix(const std::size_t &, const std::size_t &, std::vector< std::vector <_Tp> > &);
	__stdcall Matrix(const Matrix<_Tp>&);

	__stdcall ~Matrix<_Tp>() throw() {};

	_Tp& __stdcall operator ()(const std::size_t &, const std::size_t &);
	_Tp __stdcall operator ()(const std::size_t &, const std::size_t &) const;
	Matrix<_Tp>& __stdcall operator=(const Matrix<_Tp>&) throw();
	Matrix<_Tp>& __stdcall operator+=(const Matrix<_Tp>&);
	Matrix<_Tp>& __stdcall operator-=(const Matrix<_Tp>&);
	Matrix<_Tp>& __stdcall operator*=(const Matrix<_Tp>&);
	template<typename _Scalar> Matrix<_Tp>& __stdcall operator*=(const _Scalar&);
	template<typename _Scalar> Matrix<_Tp>& __stdcall operator/=(const _Scalar&);

	Matrix<_Tp> __stdcall operator+(const Matrix<_Tp>&) const;
	Matrix<_Tp> __stdcall operator-(const Matrix<_Tp>&) const;
	Matrix<_Tp> __stdcall operator*(const Matrix<_Tp>&) const;
};

#line 75 ".\\include\\matrix.h"
#line 8 ".\\src\\main.cpp"
#line 1 ".\\include\\matrix_exception.h"








class __declspec(dllimport) MatrixException : public std::exception {
public:
	enum class MatrixError : size_t {
		MATRIX_NO_ERROR,
		MATRIX_NOT_INITIALIZED = 600,
		MATRIX_INVALID_INDEX,
		MATRIX_INCOMPATIBLE_ORDER,
		MATRIX_NOT_SQUARE,
	};

	explicit MatrixException(MatrixError reason) : errorCode(reason) {};
	virtual ~MatrixException() throw() {};
	virtual const char* what() const throw ();

private:
	MatrixError errorCode;
};


#line 29 ".\\include\\matrix_exception.h"
#line 9 ".\\src\\main.cpp"
#line 1 ".\\include\\matrix_io.h"












template<typename _Tp> Matrix<_Tp> read_matrix();
template<typename _Tp> void print_matrix(const Matrix<_Tp>&) throw();


#line 18 ".\\include\\matrix_io.h"
#line 10 ".\\src\\main.cpp"


template<typename _Tp> void add_matrix();
template<typename _Tp> void multiply_matrix();
template<typename _Tp> void transpose_matrix();
template<typename _Tp> void adjoint_matrix();
template<typename _Tp, typename _Dt> void inverse_matrix();
template<typename _Tp> void determinant_matrix();




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
				if (menuAction.at(choice - 1) != 0) {
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

template<typename _Tp> void add_matrix() {
	Matrix<_Tp> A, B, C;

	
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	
	std::cout << "\nEnter details for Matrix B :-\n";
	B = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);
	std::cout << "\nMatrix B :-";
	print_matrix<_Tp>(B);

	
	C = A + B;

	
	std::cout << "\nSum :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp> void multiply_matrix() {
	Matrix<_Tp> A, B, C;

	
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	
	std::cout << "\nEnter details for Matrix B :-\n";
	B = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);
	std::cout << "\nMatrix B :-";
	print_matrix<_Tp>(B);

	
	C = A * B;

	
	std::cout << "\nProduct :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp> void transpose_matrix() {
	Matrix<_Tp> A, C;

	
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);
	
	
	C = transpose<_Tp>(A);

	
	std::cout << "\nTranspose :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp> void adjoint_matrix() {
	Matrix<_Tp> A, C;

	
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);

	
	C = adjoint<_Tp>(A);

	
	std::cout << "\nAdjoint :-";
	print_matrix<_Tp>(C);
}

template<typename _Tp, typename _Dt> void inverse_matrix() {
	Matrix<_Tp> A;
	Matrix<_Dt> C;

	
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);

	
	C = inverse<_Tp, _Dt>(A);

	
	std::cout << "\nInverse :-";
	print_matrix<_Dt>(C);
}

template<typename _Tp> void determinant_matrix() {
	Matrix<_Tp> A;
	_Tp det;

	
	std::cout << "\nEnter details for Matrix A :-\n";
	A = read_matrix<_Tp>();

	std::cout << "\nMatrix A :-";
	print_matrix<_Tp>(A);

	
	det = determinant<_Tp>(A);

	
	std::cout << "\nDeterminant :-\n";
	std::cout << det << "\n";
}