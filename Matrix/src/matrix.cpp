/*-------------------------- Include Directives -------------------------*/
#include <vector>
#include <thread>

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
	register std::size_t i;
	std::vector <std::thread> th;

	if ((A.getM() != this->getM()) || (A.getN() != this->getN())) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INCOMPATIBLE_ORDER);
	}

	for (i = 0; i < this->m; ++i) {
		// Capture everything by value.
		th.push_back(std::thread([=]() {
			register std::size_t j;

			for (j = 0; j < this->n; ++j) {
				(*this)(i, j) += A(i, j);
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

	return *this;
}

template<typename _Tp> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator-=(const Matrix<_Tp>& A) {
	register std::size_t i;
	std::vector <std::thread> th;

	if ((A.getM() != this->getM()) || (A.getN() != this->getN())) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INCOMPATIBLE_ORDER);
	}

	for (i = 0; i < this->m; ++i) {
		// Capture everything by value.
		th.push_back(std::thread([=]() {
			register std::size_t j;

			for (j = 0; j < this->n; ++j) {
				(*this)(i, j) -= A(i, j);
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

	return *this;
}

template<typename _Tp> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator*=(const Matrix<_Tp>& A) {
	Matrix<_Tp> C(this->getM(), A.getN());
	register std::size_t i;
	std::vector <std::thread> th;

	if (A.getM() != this->getN()) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_INCOMPATIBLE_ORDER);
	}

	for (i = 0; i < C.getM(); ++i) {
		// Capture everything by value.
		th.push_back(std::thread([=, &C]() {
			register std::size_t j;

			for (j = 0; j < C.getN(); ++j) {
				std::size_t k;
				_Tp sum = 0;
				for (k = 0; k < A.getM(); ++k) {
					sum += A(k, j) * (*this)(i, k);
				}
				C(i, j) = sum;
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

	/*for (i = 0; i < C.getM(); ++i) {
		std::size_t j;
		for (j = 0; j < C.getN(); ++j) {
			std::size_t k;
			_Tp sum = 0;
			for (k = 0; k < A.getM(); ++k) {
				sum += A(k, j) * (*this)(i, k);
			}
			C(i, j) = sum;
		}
	}*/

	this->m = C.getM();
	this->n = C.getN();
	this->data.clear();
	this->data.resize(this->m, std::vector<_Tp>(this->n, 0));
	this->data = C.data;

	return *this;
}

template<typename _Tp> template<typename _Scalar> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator*=(const _Scalar& k) {
	register std::size_t i;
	std::vector <std::thread> th;

	if (k == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_DIVIDE_BY_ZERO);
	}

	for (i = 0; i < this->getM(); ++i) {
		th.push_back(std::thread([=]() {
			register std::size_t j;

			for (j = 0; j < this->getN(); ++j) {
				using type = typename std::conditional<sizeof(_Scalar) >= sizeof(_Tp), _Scalar, _Tp>::type;
				type t = static_cast<type>((*this)(i, j));
				type u = static_cast<type>(k);
				(*this)(i, j) = static_cast<_Tp>(t * u);
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

	return *this;
}

template<typename _Tp> template<typename _Scalar> Matrix<_Tp>& MATRIX_CALL Matrix<_Tp> :: operator/=(const _Scalar& k) {
	register std::size_t i;
	std::vector <std::thread> th;

	if (k == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_DIVIDE_BY_ZERO);
	}

	for (i = 0; i < this->getM(); ++i) {
		th.push_back(std::thread([=]() {
			register std::size_t j;

			for (j = 0; j < this->getN(); ++j) {
				using type = typename std::conditional<sizeof(_Scalar) >= sizeof(_Tp), _Scalar, _Tp>::type;
				type t = static_cast<type>((*this)(i, j));
				type u = static_cast<type>(k);
				(*this)(i, j) = static_cast<_Tp>(t / u);
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

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

template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL operator* (const Matrix<_Tp>& A, const _Dt& k) {
	Matrix<_Dt> C(A.getM(), A.getN());
	register std::size_t i;
	std::vector <std::thread> th;

	for (i = 0; i < A.getM(); ++i) {
		th.push_back(std::thread([=, &C]() {
			register std::size_t j;

			for (j = 0; j < A.getN(); ++j) {
				C(i, j) = static_cast<_Dt>(A(i, j)) * k;
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

	return C;
}

template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL operator/ (const Matrix<_Tp>& A, const _Dt& k) {
	Matrix<_Dt> C(A.getM(), A.getN());
	register std::size_t i;
	std::vector <std::thread> th;

	for (i = 0; i < A.getM(); ++i) {
		th.push_back(std::thread([=, &C]() {
			register std::size_t j;

			for (j = 0; j < A.getN(); ++j) {
				C(i, j) = static_cast<_Dt>(A(i, j)) / k;
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

	return C;
}

template<typename _Tp> MATRIX_API bool MATRIX_CALL operator==(const Matrix<_Tp>& A, const Matrix<_Tp>& B) {
	volatile bool res;

	if ((A.getM() == 0 || B.getM() == 0) || (A.getN() == 0 || B.getN() == 0)) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	res = true;
	if ((A.getM() == B.getM()) && (A.getN() == B.getN())) {
		register std::size_t i;
		std::vector <std::thread> th;

		for (i = 0; i < A.getM(); ++i) {
			th.push_back(std::thread([=, &res]() {
				register std::size_t j;

				for (j = 0; j < A.getN(); ++j) {
					res = ((A(i, j) == B(i, j))? res: false); // Race eliminated...
					if (res == false) {
						break;
					}
				}
			}));
			if (res == false) {
				break;
			}
		}
		for(std::thread &t : th) {
			t.join();
		};
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
	register std::size_t i;
	std::vector <std::thread> th;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	for (i = 0; i < A.getN(); ++i) {
		th.push_back(std::thread([=, &C]() {
			register std::size_t j;

			for (j = 0; j < A.getM(); ++j) {
				C(i, j) = A(j, i);
			}
		}));
	}
	for(std::thread &t : th) {
		t.join();
	};

	return C;
}


template<typename _Tp> static Matrix<_Tp> MATRIX_CALL pivot(Matrix<_Tp>& A) {
	Matrix<_Tp> C(A);
	std::size_t i;
	std::size_t n;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}

	n = A.getM();
	for (i = 0; i < n; ++i) {
	}

	return C;
}


template<typename _Tp> MATRIX_API _Tp MATRIX_CALL determinant(const Matrix<_Tp>& A) {
	std::size_t i, j, k;
	std::size_t n;
	Matrix<_Tp> L(A.getM(), A.getN());
	Matrix<_Tp> U(A.getM(), A.getN());
	_Tp det;

	if (A.getM() == 0 || A.getN() == 0) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_INITIALIZED);
	}
	else if (A.getM() != A.getN()) {
		throw MatrixException(MatrixException::MatrixError::MATRIX_NOT_SQUARE);
	}

	n = A.getM();
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (j < i) {
				L(j, i) = 0;
			} else {
				L(j, i) = A(j, i);
				for (k = 0; k < i; k++) {
					L(j, i) = L(j, i) - L(j, k) * U(k, i);
				}
			}
		}
		for (j = 0; j < n; ++j) {
			if (j < i) {
				U(i, j) = 0;
			}
			else {
				if (j == i) {
					U(i, j) = 1;
				}
				else {
					U(i, j) = A(i, j) / L(i, i);
					for (k = 0; k < i; ++k) {
						U(i, j) = U(i, j) - ((L(i, k) * U(k, j)) / L(i, i));
					}
				}
			}
		}
	}

	det = 1;
	for (i = 0; i < n; ++i) {
		det *= L(i, i) *  U(i, i);
	}

	return det;
}


/*template<typename _Tp> MATRIX_API _Tp MATRIX_CALL determinant(const Matrix<_Tp>& A) {
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

	if (n < 1) { // Error
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

			// Recursively calculate determinant of sub-matrix
			_Tp ptr = determinant<_Tp>(C);

			det += sign * A(0, j1) * ptr;
			sign *= -1;
		}
	}

	return det;
}*/

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
#define SPEC_DECL_A(_A) \
	template class MATRIX_API Matrix<_A>; \
	template MATRIX_API _A MATRIX_CALL determinant(const Matrix<_A>&); \
	template MATRIX_API Matrix<_A> MATRIX_CALL transpose(const Matrix<_A>&); \
	template MATRIX_API Matrix<_A> MATRIX_CALL cofactor(const Matrix<_A>&); \
	template MATRIX_API Matrix<_A> MATRIX_CALL adjoint(const Matrix<_A>&); \
	template MATRIX_API Matrix<_A> MATRIX_CALL inverse(const Matrix<_A>&); \
	template MATRIX_API Matrix<_A> MATRIX_CALL operator*(const Matrix<_A>&, const _A&); \
	template MATRIX_API Matrix<_A> MATRIX_CALL operator/(const Matrix<_A>&, const _A&); \
	template Matrix<_A>& MATRIX_CALL Matrix<_A> :: operator*=(const _A&); \
	template Matrix<_A>& MATRIX_CALL Matrix<_A> :: operator/=(const _A&); \
	template MATRIX_API bool MATRIX_CALL operator==(const Matrix<_A>&, const Matrix<_A>&); \
	template MATRIX_API bool MATRIX_CALL operator!=(const Matrix<_A>&, const Matrix<_A>&);


#define MATRIX_SPECIALIZATION(_A, _B) \
	static_assert (sizeof(_A) <= sizeof(_B), "Expects (_A <= _B) ..."); \
	SPEC_DECL_A(_A) \
	SPEC_DECL_A(_B) \
	template Matrix<_B>& MATRIX_CALL Matrix<_B> :: operator*=(const _A&); \
	template Matrix<_B>& MATRIX_CALL Matrix<_B> :: operator/=(const _A&); \
	template MATRIX_API Matrix<_B> MATRIX_CALL inverse(const Matrix<_A>&); \
	template MATRIX_API Matrix<_A> MATRIX_CALL inverse(const Matrix<_B>&);


MATRIX_SPECIALIZATION(float, double)

/* End of Specializations */