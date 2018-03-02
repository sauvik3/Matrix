# Matrix : A Generic Matrix Library in C++  [![Build Status](https://travis-ci.org/sauvik3/Matrix.svg?branch=master)](https://travis-ci.org/sauvik3/Matrix)&nbsp; [![Build status](https://ci.appveyor.com/api/projects/status/d2jnxclg20vd0o50?svg=true)](https://ci.appveyor.com/project/sauvik3/matrix)&nbsp; [![Coverage Status](https://coveralls.io/repos/github/sauvik3/Matrix/badge.svg?branch=master)](https://coveralls.io/github/sauvik3/Matrix?branch=master)
Description
-------------------------------------------------------

Builds as a C++ Library Module for performing following
basic Matrix Operations:

1. Matrix Addition
2. Matrix Subtraction
3. Matrix Multiplication
4. Scalar Matrix Multiplication
5. Scalar Matrix Division
6. Matrix Determinant
7. Matrix Transpose
8. Matrix Cofactor
9. Matrix Adjoint
10. Matrix Inverse

Functionality is implemented as generic templatized class Matrix.
Operations are implemented using intitutive Operator Overloading
and generic constructs.

Usage
-------------------------------------------------------

To use the library in your program, you must build it and include
respective header files (matrix.h, matrix_exception.h, matrix_error_messages.h)
and link against the static/export library (matrix.lib (,matrix.dll)).
