#ifndef __MATRIX_EXCEPTION_CPP_INCLUDED__
#define __MATRIX_EXCEPTION_CPP_INCLUDED__

#include <string>

#include "matrix_exception.h"
#include "matrix_error_messages.h"

const char * MatrixException::what() const throw()
{
	static std::string errMsg("");

	switch (this->errorCode)
	{
	case MatrixError::MATRIX_NO_ERROR:	errMsg = STR_MATRIX_NO_ERROR;
		break;
	case MatrixError::MATRIX_NOT_INITIALIZED:	errMsg = STR_MATRIX_NOT_INITIALIZED;
		break;
	case MatrixError::MATRIX_INVALID_INDEX:	errMsg = STR_MATRIX_INVALID_INDEX;
		break;
	case MatrixError::MATRIX_INCOMPATIBLE_ORDER:	errMsg = STR_MATRIX_INCOMPATIBLE_ORDER;
		break;
	case MatrixError::MATRIX_NOT_SQUARE:	errMsg = STR_MATRIX_NOT_SQUARE;
	default:
		break;
	}

	return errMsg.c_str();
}

#endif