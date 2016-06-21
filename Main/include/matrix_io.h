#ifndef __MATRIX_IO_H_INCLUDED__
#define __MATRIX_IO_H_INCLUDED__

/* Use ASCII Box-Drawing Characters */

#define LLCNR 0xC0   /* ? */
#define LUCNR 0xDA   /* ? */
#define RLCNR 0xD9   /* ? */
#define RUCNR 0xBF   /* ? */
#define VRT   0xB3   /* ? */

/*---------------------------- Prototypes -----------------------------*/
template<typename _Tp> Matrix<_Tp> read_matrix();
template<typename _Tp> void print_matrix(const Matrix<_Tp>&) throw();
/*---------------------------------------------------------------------*/

#endif