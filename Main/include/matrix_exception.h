#ifndef __MATRIX_EXCEPTION_H_INCLUDED__
#define __MATRIX_EXCEPTION_H_INCLUDED__

#include <exception>

#include "matrix.h"

/*--------------------------- Declarations ----------------------------*/
class MATRIX_API MatrixException : public std::exception {
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
/*---------------------------------------------------------------------*/

#endif