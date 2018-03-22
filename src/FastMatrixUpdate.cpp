#include "FastMatrixUpdate.hpp"


namespace DMFT
{
    namespace util
    {
    void swapRows(MatrixT *A, int from, int to)
    {
        auto P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic >(A->rows());
        P.setIdentity();
        if(to < from)
        {
            for(int i = to; i < from; i++)
            {
                P.applyTranspositionOnTheRight(i, i+1);
            }
        }
        else
        {
            for(int i = to; i > from; i--)
            {
                P.applyTranspositionOnTheRight(i-1, i);
            }
        }
        *A = (P*(*A)*(P.transpose())).eval();
    }

    void MInc(MatrixT *A, const RowVectorT &R, const VectorT &Q, const RealT Sp, const int index, bool swap)
    {
        const unsigned int n = A->rows();
        MatrixT P = *A;
        (*A).conservativeResize(n+1, n+1);
        (*A)(n,n)                = Sp;
        A->topRightCorner(n,1)   = -(P*Q)*Sp;
        A->bottomLeftCorner(1,n) = -(R*P)*Sp;
        A->topLeftCorner(n,n)    = P + A->topRightCorner(n,1) * A->bottomLeftCorner(1,n)/Sp;
        if(swap)
            swapRows(A, n, index);
    } 

    void MDec(MatrixT* A, const unsigned int index)
    {
        const unsigned int n = A->rows();
        if(n > 1)
        {
            swapRows(A, index, n-1);
            A->topLeftCorner(n-1,n-1) = (A->topLeftCorner(n-1,n-1) - A->topRightCorner(n-1,1)*A->bottomLeftCorner(1,n-1)  / (*A)(n-1,n-1)).eval();
        }
        A->conservativeResize(n-1,n-1);
    } 
    }
}
