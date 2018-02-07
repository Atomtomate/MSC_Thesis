#ifndef
#define IGFMATRIX_HPP_

#include "Config.hpp"
#include "GreensFct.hpp"

template<unsigned int NFLAVORS>
class IGFMatrix
{

    void swapRows(unsigned int from, unsigned int to, unsigned int f)
    {
        auto P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic >(M[f].rows());
        P.setIdentity();
        if(to < from)
        {
            for(int i = to; i < from; i++)
                P.applyTranspositionOnTheRight(i, i+1);
        }
        else
        {
            for(int i = to; i > from; i--)
                P.applyTranspositionOnTheRight(i-1, i);
        }
        M[f] = (P*M[f]*(P.transpose())).eval();
    }

    void MInc(const RowVectorT &R, const VectorT &Q, const RealT Sp, const int index, unsigned int f)
    {
        const unsigned int n = M[f].rows();
        MatrixT P = M[f];
        M[f].conservativeResize(n+1, n+1);
        M[f](n,n)                = Sp;
        M[f].topRightCorner(n,1)   = -(P*Q)*Sp;
        M[f].bottomLeftCorner(1,n) = -(R*P)*Sp;
        M[f].topLeftCorner(n,n)    = P + M[f].topRightCorner(n,1) * M[f].bottomLeftCorner(1,n)/Sp;
        swapRows(n, index, f);
    } 
    

    void MDec(const unsigned int index, unsigned int f)
    {
        const unsigned int n = M[f].rows();
        // swap rows to obtain the correct position
        if(n > 1)
        {
            swapRows(index, n-1, f);
            M[f].topLeftCorner(n-1,n-1) = (M[f].topLeftCorner(n-1,n-1) - M[f].topRightCorner(n-1,1)*M[f].bottomLeftCorner(1,n-1)  / M[f](n-1,n-1)).eval();
        }
        M[f].conservativeResize(n-1,n-1);
    } 

    private:
        std::array<MatrixT, NFLAVORS> M;
};

#endif
