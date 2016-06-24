/*-------------------------- Include Directives -------------------------*/
#include <vector>

#include "matrix.h"
#include "matrix_exception.h"
/*-----------------------------------------------------------------------*/


/*------------------------------- Getters -------------------------------*/
template<typename _Tp> std::size_t MATRIX_CALL Matrix<_Tp> ::getM() const throw() {
	return this->m;
}

template<typename _Tp> std::size_t MATRIX_CALL Matrix<_Tp> ::getN() const throw() {
	return this->n;
}

template<typename _Tp> std::vector< std::vector< _Tp > > MATRIX_CALL Matrix<_Tp> ::getData() const throw() {
	return this->data;
}
/*-----------------------------------------------------------------------*/


/*------------------------------ Constructors ---------------------------*/
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
/*-----------------------------------------------------------------------*/


/*------------------------------- Operators -----------------------------*/
template<typename _Tp> _Tp& MATRIX_CALL Matrix<_Tp> :: operator ()(const std::size_t &i, const std::size_t &j) {
	if (this->getM() == 0 || this->getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	if (this->getM() < i || this->getN() < j) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INVALID_INDEX);
	}
	return this->data.at(i).at(j);
}

template<typename _Tp> _Tp MATRIX_CALL Matrix<_Tp> :: operator ()(const std::size_t &i, const std::size_t &j) const {
	if (this->getM() == 0 || this->getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	if (this->getM() < i || this->getN() < j) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INVALID_INDEX);
	}
	return this->data.at(i).at(j);
}

template<typename _Tp> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator=(const Matrix<_Tp>& A) throw() {
	this->m = A.getM();
	this->n = A.getN();
	this->data = A.getData();
	return *this;
}

template<typename _Tp> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator+=(const Matrix<_Tp>& A) {
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

template<typename _Tp> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator-=(const Matrix<_Tp>& A) {
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

template<typename _Tp> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator*=(const Matrix<_Tp>& A) {
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

template<typename _Tp> template<typename _Scalar> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator*=(const _Scalar& k) {
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

template<typename _Tp> template<typename _Scalar> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator/=(const _Scalar& k) {
	std::size_t i;

	if (k == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_DIVIDE_BY_ZERO);
	}

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

template<typename _Tp> Matrix<_Tp> MATRIX_CALL Matrix<_Tp> :: operator+(const Matrix<_Tp>& A) const {
	Matrix<_Tp> C(*this);
	C += A;

	return C;
}

template<typename _Tp> Matrix<_Tp> MATRIX_CALL Matrix<_Tp> :: operator-(const Matrix<_Tp>& A) const {
	Matrix<_Tp> C(*this);
	C -= A;

	return C;
}

template<typename _Tp> Matrix<_Tp> MATRIX_CALL Matrix<_Tp> :: operator*(const Matrix<_Tp>& A) const {
	Matrix<_Tp> C(*this);
	C *= A;

	return C;
}

template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL operator*(const Matrix<_Tp>& A, const _Dt& k) {
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

template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL operator/(const Matrix<_Tp>& A, const _Dt& k) {
	Matrix<_Dt> C(A.getM(), A.getN());
	std::size_t i;

	if (k == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_DIVIDE_BY_ZERO);
	}

	for (i = 0; i < A.getM(); ++i) {
		std::size_t j;
		for (j = 0; j < A.getN(); ++j) {
			C(i, j) = static_cast<_Dt>(A(i, j)) / k;
		}
	}

	return C;
}

template<typename _Tp> MATRIX_API bool MATRIX_CALL operator==(const Matrix<_Tp>& A, const Matrix<_Tp>& B) {
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

template<typename _Tp> MATRIX_API bool MATRIX_CALL operator!=(const Matrix<_Tp>& A, const Matrix<_Tp>& B) {
	bool v;
	v = !(A == B);

	return v;
}
/*-----------------------------------------------------------------------*/


/*--------------------------- Matrix Functions --------------------------*/
template<typename _Tp> MATRIX_API Matrix<_Tp> MATRIX_CALL transpose(const Matrix<_Tp>& A) {
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

template<typename _Tp> MATRIX_API _Tp MATRIX_CALL determinant(const Matrix<_Tp>& A) {
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

	if (n < 1) { /* Error */
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

			/* Recursively calculate determinant of sub-matrix */
			_Tp ptr = determinant<_Tp>(C);

			det += sign * A(0, j1) * ptr;
			sign *= -1;
		}
	}

	return det;
}

template<typename _Tp> MATRIX_API Matrix<_Tp> MATRIX_CALL cofactor(const Matrix<_Tp>& A) {
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

				/* Calculate the determinant */
				det = determinant<_Tp>(C);

				/* Fill in the elements of the cofactor */
				B(i, j) = sign * det;
				sign *= -1;
			}
		}
	}

	return B;
}

template<typename _Tp> MATRIX_API Matrix<_Tp> MATRIX_CALL adjoint(const Matrix<_Tp>& A) {
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

template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL inverse(const Matrix<_Tp>& A) {
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
	} else {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INVERTIBLE);
	}

	return C;
}
/*-----------------------------------------------------------------------*/

/* Specializations */
template class MATRIX_API Matrix<float>;
template MATRIX_API float MATRIX_CALL determinant(const Matrix<float>&);
template MATRIX_API Matrix<float> MATRIX_CALL transpose(const Matrix<float>&);
template MATRIX_API Matrix<float> MATRIX_CALL cofactor(const Matrix<float>&);
template MATRIX_API Matrix<float> MATRIX_CALL adjoint(const Matrix<float>&);
template MATRIX_API Matrix<float> MATRIX_CALL inverse(const Matrix<float>&);
template MATRIX_API Matrix<double> MATRIX_CALL operator*(const Matrix<float>&, const double&);
template MATRIX_API Matrix<double> MATRIX_CALL operator/(const Matrix<float>&, const double&);
template MATRIX_API bool MATRIX_CALL operator==(const Matrix<float>&, const Matrix<float>&);
template MATRIX_API bool MATRIX_CALL operator!=(const Matrix<float>&, const Matrix<float>&);

template class MATRIX_API Matrix<double>;
template MATRIX_API double MATRIX_CALL determinant(const Matrix<double>&);
template MATRIX_API Matrix<double> MATRIX_CALL transpose(const Matrix<double>&);
template MATRIX_API Matrix<double> MATRIX_CALL cofactor(const Matrix<double>&);
template MATRIX_API Matrix<double> MATRIX_CALL adjoint(const Matrix<double>&);
template MATRIX_API Matrix<double> MATRIX_CALL inverse(const Matrix<double>&);
template MATRIX_API bool MATRIX_CALL operator==(const Matrix<double>&, const Matrix<double>&);
template MATRIX_API bool MATRIX_CALL operator!=(const Matrix<double>&, const Matrix<double>&);

template MATRIX_API Matrix<double> MATRIX_CALL inverse(const Matrix<float>&);
template MATRIX_API Matrix<float> MATRIX_CALL inverse(const Matrix<double>&);
/* End of Specializations */