#ifndef FAST_MATRIX_UPDATE_HPP
#define FAST_MATRIX_UPDATE_HPP

#include "Config.hpp"


namespace DMFT
{
    namespace util
    {
    void swapRows(MatrixT *A, int from, int to);
    
    void MInc(MatrixT *A, const RowVectorT &R, const VectorT &Q, const RealT Sp, const int index, bool swap = true);
    
    void MDec(MatrixT* A, const unsigned int index);
    }
}
#endif
