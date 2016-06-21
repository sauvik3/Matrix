#line 20 ".\\include\\matrix.h"





template<typename _Tp> class __declspec(dllexport) Matrix;

template<typename _Tp> __declspec(dllexport) _Tp __stdcall determinant(const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllexport) Matrix<_Tp> __stdcall transpose(const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllexport) Matrix<_Tp> __stdcall cofactor(const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllexport) Matrix<_Tp> __stdcall adjoint(const Matrix<_Tp>&);
template<typename _Tp, typename _Dt> __declspec(dllexport) Matrix<_Dt> __stdcall inverse(const Matrix<_Tp>&);

template<typename _Tp, typename _Dt> __declspec(dllexport) Matrix<_Dt> __stdcall operator*(const Matrix<_Tp>&, const _Dt&);
template<typename _Tp, typename _Dt> __declspec(dllexport) Matrix<_Dt> __stdcall operator/(const Matrix<_Tp>&, const _Dt&);

template<typename _Tp> __declspec(dllexport) bool __stdcall operator==(const Matrix<_Tp>&, const Matrix<_Tp>&);
template<typename _Tp> __declspec(dllexport) bool __stdcall operator!=(const Matrix<_Tp>&, const Matrix<_Tp>&);
template<typename _Tp> bool operator>(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator<(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator>=(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator<=(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;

template<typename _Tp>
class __declspec(dllexport) Matrix {
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
#line 8 ".\\src\\matrix.cpp"
#line 1 ".\\include\\matrix_exception.h"







class __declspec(dllexport) MatrixException : public std::exception {
public:
	enum class MatrixError : std::size_t {
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


#line 28 ".\\include\\matrix_exception.h"
#line 9 ".\\src\\matrix.cpp"




template<typename _Tp> std::size_t __stdcall Matrix<_Tp> ::getM() const throw() {
	return this->m;
}

template<typename _Tp> std::size_t __stdcall Matrix<_Tp> ::getN() const throw() {
	return this->n;
}

template<typename _Tp> std::vector< std::vector< _Tp > > __stdcall Matrix<_Tp> ::getData() const throw() {
	return this->data;
}




template<typename _Tp> Matrix<_Tp> ::Matrix() : m(0), n(0), data(0, std::vector<_Tp>(0)) {}

template<typename _Tp> Matrix<_Tp> ::Matrix(const std::size_t &_m, const std::size_t &_n) :
	m(_m),
	n(_n),
	data(_m, std::vector<_Tp>(_n)) {}

template<typename _Tp> Matrix<_Tp> ::Matrix(const std::size_t &_m, const std::size_t &_n, std::vector< std::vector <_Tp> > &_data) :
	m(_m),
	n(_n),
	data(_data) {}

template<typename _Tp> Matrix<_Tp> ::Matrix(const Matrix<_Tp>& A) :
	m(A.getM()),
	n(A.getN())
{
	data = A.data;
}




template<typename _Tp> _Tp& __stdcall Matrix<_Tp> :: operator ()(const std::size_t &i, const std::size_t &j) {
	if (this->getM() == 0 || this->getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	if (this->getM() < i || this->getN() < j) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INVALID_INDEX);
	}
	return this->data.at(i).at(j);
}

template<typename _Tp> _Tp __stdcall Matrix<_Tp> :: operator ()(const std::size_t &i, const std::size_t &j) const {
	if (this->getM() == 0 || this->getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	if (this->getM() < i || this->getN() < j) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INVALID_INDEX);
	}
	return this->data.at(i).at(j);
}

template<typename _Tp> Matrix<_Tp>& __stdcall Matrix<_Tp> :: operator=(const Matrix<_Tp>& A) throw() {
	this->m = A.getM();
	this->n = A.getN();
	this->data = A.getData();
	return *this;
}

template<typename _Tp> Matrix<_Tp>& __stdcall Matrix<_Tp> :: operator+=(const Matrix<_Tp>& A) {
	std::size_t i;

	if ((A.getM() != this->getM()) || (A.getN() != this->getN())) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INCOMPATIBLE_ORDER);
	}

	for (i = 0; i < this->m; ++i) {
		std::size_t j;
		for (j = 0; j < this->n; ++j) {
			(*this)(i, j) += A(i, j);
		}
	}

	return *this;
}

template<typename _Tp> Matrix<_Tp>& __stdcall Matrix<_Tp> :: operator-=(const Matrix<_Tp>& A) {
	std::size_t i;

	if ((A.getM() != this->getM()) || (A.getN() != this->getN())) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INCOMPATIBLE_ORDER);
	}

	for (i = 0; i < this->m; ++i) {
		std::size_t j;
		for (j = 0; j < this->n; ++j) {
			(*this)(i, j) -= A(i, j);
		}
	}

	return *this;
}

template<typename _Tp> Matrix<_Tp>& __stdcall Matrix<_Tp> :: operator*=(const Matrix<_Tp>& A) {
	Matrix<_Tp> C(this->getM(), A.getN());
	std::size_t i;

	if (A.getM() != this->getN()) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INCOMPATIBLE_ORDER);
	}

	for (i = 0; i < C.getM(); ++i) {
		std::size_t j;
		for (j = 0; j < C.getN(); ++j) {
			std::size_t k;
			_Tp sum = 0;
			for (k = 0; k < A.getM(); ++k) {
				sum += A(k, j) * (*this)(i, k);
			}
			C(i, j) = sum;
		}
	}

	this->m = C.getM();
	this->n = C.getN();
	this->data.clear();
	this->data.resize(this->m, std::vector<_Tp>(this->n, 0));
	this->data = C.data;

	return *this;
}

template<typename _Tp> template<typename _Scalar> Matrix<_Tp>& __stdcall Matrix<_Tp> :: operator*=(const _Scalar& k) {
	std::size_t i;

	for (i = 0; i < this->getM(); ++i) {
		std::size_t j;
		for (j = 0; j < this->getN(); ++j) {
			using type = typename std::conditional<sizeof(_Scalar) >= sizeof(_Tp), _Scalar, _Tp>::type;
			type t = static_cast<type>((*this)(i, j));
			type u = static_cast<type>(k);
			(*this)(i, j) = static_cast<_Tp>(t*u);
		}
	}

	return *this;
}

template<typename _Tp> template<typename _Scalar> Matrix<_Tp>& __stdcall Matrix<_Tp> :: operator/=(const _Scalar& k) {
	std::size_t i;

	for (i = 0; i < this->getM(); ++i) {
		std::size_t j;
		for (j = 0; j < this->getN(); ++j) {
			using type = typename std::conditional<sizeof(_Scalar) >= sizeof(_Tp), _Scalar, _Tp>::type;
			type t = static_cast<type>((*this)(i, j));
			type u = static_cast<type>(k);
			(*this)(i, j) = static_cast<_Tp>(t / u);
		}
	}

	return *this;
}

template<typename _Tp> Matrix<_Tp> __stdcall Matrix<_Tp> :: operator+(const Matrix<_Tp>& A) const {
	Matrix<_Tp> C(*this);
	C += A;

	return C;
}

template<typename _Tp> Matrix<_Tp> __stdcall Matrix<_Tp> :: operator-(const Matrix<_Tp>& A) const {
	Matrix<_Tp> C(*this);
	C -= A;

	return C;
}

template<typename _Tp> Matrix<_Tp> __stdcall Matrix<_Tp> :: operator*(const Matrix<_Tp>& A) const {
	Matrix<_Tp> C(*this);
	C *= A;

	return C;
}

template<typename _Tp, typename _Dt> __declspec(dllexport) Matrix<_Dt> __stdcall operator*(const Matrix<_Tp>& A, const _Dt& k) {
	Matrix<_Dt> C(A.getM(), A.getN());
	std::size_t i;

	for (i = 0; i < A.getM(); ++i) {
		std::size_t j;
		for (j = 0; j < A.getN(); ++j) {
			C(i, j) = static_cast<_Dt>(A(i, j)) * k;
		}
	}

	return C;
}

template<typename _Tp, typename _Dt> __declspec(dllexport) Matrix<_Dt> __stdcall operator/(const Matrix<_Tp>& A, const _Dt& k) {
	Matrix<_Dt> C(A.getM(), A.getN());
	std::size_t i;

	for (i = 0; i < A.getM(); ++i) {
		std::size_t j;
		for (j = 0; j < A.getN(); ++j) {
			C(i, j) = static_cast<_Dt>(A(i, j)) / k;
		}
	}

	return C;
}

template<typename _Tp> __declspec(dllexport) bool __stdcall operator==(const Matrix<_Tp>& A, const Matrix<_Tp>& B) {
	bool res;

	if ((A.getM() == 0 || B.getM() == 0) || (A.getN() == 0 || B.getN() == 0)) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	res = true;
	if ((A.getM() == B.getM()) && (A.getN() == B.getN())) {
		std::size_t i;

		for (i = 0; i < A.getM(); ++i) {
			std::size_t j;
			for (j = 0; j < A.getN(); ++j) {
				res = (A(i, j) == B(i, j));
				if (res == false) {
					break;
				}
			}
			if (res == false) {
				break;
			}
		}
	}
	else {
		res = false;
	}

	return res;
}

template<typename _Tp> __declspec(dllexport) bool __stdcall operator!=(const Matrix<_Tp>& A, const Matrix<_Tp>& B) {
	bool v;
	v = !(A == B);

	return v;
}




template<typename _Tp> __declspec(dllexport) Matrix<_Tp> __stdcall transpose(const Matrix<_Tp>& A) {
	Matrix<_Tp> C(A.getN(), A.getM());
	std::size_t i;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	for (i = 0; i < A.getN(); ++i) {
		std::size_t j;
		for (j = 0; j < A.getM(); ++j) {
			C(i, j) = A(j, i);
		}
	}

	return C;
}

template<typename _Tp> __declspec(dllexport) _Tp __stdcall determinant(const Matrix<_Tp>& A) {
	_Tp det;
	std::size_t j, j1;
	std::size_t n;
	int sign;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}
	else if (A.getM() != A.getN()) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_SQUARE);
	}

	n = A.getM();

	if (n < 1) { 
		return 0;
	}
	else if (n == 1) {
		det = A(0, 0);
	}
	else if (n == 2) {
		det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
	}
	else {
		det = 0;
		sign = 1;
		for (j1 = 0; j1 < n; j1++) {
			Matrix<_Tp> C(n - 1, n - 1);
			std::size_t i;

			for (i = 1; i < n; ++i) {
				std::size_t j2 = 0;
				for (j = 0; j < n; j++) {
					if (j == j1) {
						continue;
					}
					C(i - 1, j2) = A(i, j);
					++j2;
				}
			}

			
			_Tp ptr = determinant<_Tp>(C);

			det += sign * A(0, j1) * ptr;
			sign *= -1;
		}
	}

	return det;
}

template<typename _Tp> __declspec(dllexport) Matrix<_Tp> __stdcall cofactor(const Matrix<_Tp>& A) {
	Matrix<_Tp> B(A.getM(), A.getN());
	std::size_t n;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}
	else if (A.getM() != A.getN()) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_SQUARE);
	}
	n = A.getM();

	if (n == 1) {
		B(0, 0) = 1;
	}
	else {
		Matrix<_Tp> C(n - 1, n - 1);
		std::size_t i, j, ii, jj;
		int sign;

		sign = 1;
		for (j = 0; j < n; ++j) {
			for (i = 0; i < n; ++i) {
				_Tp det;
				int i1 = 0;

				for (ii = 0; ii < n; ++ii) {
					int j1 = 0;
					if (ii == i) {
						continue;
					}
					for (jj = 0; jj < n; ++jj) {
						if (jj == j) {
							continue;
						}
						C(i1, j1) = A(ii, jj);
						++j1;
					}
					++i1;
				}

				
				det = determinant<_Tp>(C);

				
				B(i, j) = sign * det;
				sign *= -1;
			}
		}
	}

	return B;
}

template<typename _Tp> __declspec(dllexport) Matrix<_Tp> __stdcall adjoint(const Matrix<_Tp>& A) {
	Matrix<_Tp> B;
	Matrix<_Tp> C;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}
	else if (A.getM() != A.getN()) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_SQUARE);
	}

	B = cofactor(A);
	C = transpose(B);

	return C;
}

template<typename _Tp, typename _Dt> __declspec(dllexport) Matrix<_Dt> __stdcall inverse(const Matrix<_Tp>& A) {
	Matrix<_Tp> B;
	Matrix<_Dt> C;
	_Tp det;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}
	else if (A.getM() != A.getN()) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_SQUARE);
	}

	det = determinant<_Tp>(A);
	B = adjoint(A);
	if (det != 0) {
		C = B / static_cast<_Dt>(det);
	}

	return C;
}



template class __declspec(dllexport) Matrix<float>;
template __declspec(dllexport) float __stdcall determinant(const Matrix<float>&);
template __declspec(dllexport) Matrix<float> __stdcall transpose(const Matrix<float>&);
template __declspec(dllexport) Matrix<float> __stdcall cofactor(const Matrix<float>&);
template __declspec(dllexport) Matrix<float> __stdcall adjoint(const Matrix<float>&);
template __declspec(dllexport) Matrix<float> __stdcall inverse(const Matrix<float>&);
template __declspec(dllexport) Matrix<double> __stdcall operator*(const Matrix<float>&, const double&);
template __declspec(dllexport) Matrix<double> __stdcall operator/(const Matrix<float>&, const double&);
template __declspec(dllexport) bool __stdcall operator==(const Matrix<float>&, const Matrix<float>&);
template __declspec(dllexport) bool __stdcall operator!=(const Matrix<float>&, const Matrix<float>&);

template class __declspec(dllexport) Matrix<double>;
template __declspec(dllexport) double __stdcall determinant(const Matrix<double>&);
template __declspec(dllexport) Matrix<double> __stdcall transpose(const Matrix<double>&);
template __declspec(dllexport) Matrix<double> __stdcall cofactor(const Matrix<double>&);
template __declspec(dllexport) Matrix<double> __stdcall adjoint(const Matrix<double>&);
template __declspec(dllexport) Matrix<double> __stdcall inverse(const Matrix<double>&);
template __declspec(dllexport) bool __stdcall operator==(const Matrix<double>&, const Matrix<double>&);
template __declspec(dllexport) bool __stdcall operator!=(const Matrix<double>&, const Matrix<double>&);

template __declspec(dllexport) Matrix<double> __stdcall inverse(const Matrix<float>&);
template __declspec(dllexport) Matrix<float> __stdcall inverse(const Matrix<double>&);


#line 454 ".\\src\\matrix.cpp"
