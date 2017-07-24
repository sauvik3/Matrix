#ifndef __MATRIX_H_INCLUDED__
#define __MATRIX_H_INCLUDED__

#ifdef _WIN32
	#ifdef MATRIX_USE_DLL
		#ifdef DLL_EXPORTS
			#define MATRIX_API __declspec(dllexport)
		#else
			#define MATRIX_API __declspec(dllimport)
		#endif
	#else
		#define MATRIX_API
	#endif

	#define MATRIX_CALL __stdcall
#else
	#define MATRIX_API
	#define MATRIX_CALL
#endif

#include <vector>

/*--------------------------- Declarations ----------------------------*/
// Forward declaration
template<typename _Tp> class MATRIX_API Matrix;

template<typename _Tp> MATRIX_API _Tp MATRIX_CALL determinant(const Matrix<_Tp>&);
template<typename _Tp> MATRIX_API Matrix<_Tp> MATRIX_CALL transpose(const Matrix<_Tp>&);
template<typename _Tp> MATRIX_API Matrix<_Tp> MATRIX_CALL cofactor(const Matrix<_Tp>&);
template<typename _Tp> MATRIX_API Matrix<_Tp> MATRIX_CALL adjoint(const Matrix<_Tp>&);
template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL inverse(const Matrix<_Tp>&);

template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL operator*(const Matrix<_Tp>&, const _Dt&);
template<typename _Tp, typename _Dt> MATRIX_API Matrix<_Dt> MATRIX_CALL operator/(const Matrix<_Tp>&, const _Dt&);

template<typename _Tp> MATRIX_API bool MATRIX_CALL operator==(const Matrix<_Tp>&, const Matrix<_Tp>&);
template<typename _Tp> MATRIX_API bool MATRIX_CALL operator!=(const Matrix<_Tp>&, const Matrix<_Tp>&);
template<typename _Tp> bool operator>(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator<(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator>=(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;
template<typename _Tp> bool operator<=(const Matrix<_Tp>&, const Matrix<_Tp>&) = delete;

template<typename _Tp>
class MATRIX_API Matrix {
	std::size_t m;
	std::size_t n;
	std::vector< std::vector <_Tp> > data;
public:
	std::size_t MATRIX_CALL getM() const throw();
	std::size_t MATRIX_CALL getN() const throw();
	std::vector< std::vector< _Tp > > MATRIX_CALL getData() const throw();

	MATRIX_CALL Matrix();
	MATRIX_CALL Matrix(const std::size_t &, const std::size_t &);
	MATRIX_CALL Matrix(const std::size_t &, const std::size_t &, std::vector< std::vector <_Tp> > &);
	MATRIX_CALL Matrix(const Matrix<_Tp>&);

	MATRIX_CALL ~Matrix<_Tp>() throw() {};

	_Tp& MATRIX_CALL operator ()(const std::size_t &, const std::size_t &);
	_Tp MATRIX_CALL operator ()(const std::size_t &, const std::size_t &) const;
	Matrix<_Tp>& MATRIX_CALL operator=(const Matrix<_Tp>&) throw();
	Matrix<_Tp>& MATRIX_CALL operator+=(const Matrix<_Tp>&);
	Matrix<_Tp>& MATRIX_CALL operator-=(const Matrix<_Tp>&);
	Matrix<_Tp>& MATRIX_CALL operator*=(const Matrix<_Tp>&);
	template<typename _Scalar> Matrix<_Tp>& MATRIX_CALL operator*=(const _Scalar&);
	template<typename _Scalar> Matrix<_Tp>& MATRIX_CALL operator/=(const _Scalar&);

	Matrix<_Tp> MATRIX_CALL operator+(const Matrix<_Tp>&) const;
	Matrix<_Tp> MATRIX_CALL operator-(const Matrix<_Tp>&) const;
	Matrix<_Tp> MATRIX_CALL operator*(const Matrix<_Tp>&) const;
};

#endif