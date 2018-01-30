#include "StrongCoupling.hpp"
#include <iostream>


namespace  DMFT
{
        StrongCoupling::StrongCoupling(GreensFct &hybr, GreensFct &gImp, const Config& config, const unsigned int burninSteps):
        hyb(hyb), gImp(gImp), conf(config), burninSteps(burninSteps), segments(config)
        {
            totN 		= 0;
            steps 		= 0;
            totalSign	= 0;
            lastSign	= 1;
            r_time.split(6, 0);
            r_spin.split(6, 1);
            r_insert.split(6, 2);
            r_accept.split(6, 3);
            r_shift.split(6, 4);
            r_timep.split(6, 5);
            r_time.split(config.local.size() , config.local.rank());           // choose sub−stream no. rank out of size streams
            r_spin.split(config.local.size() , config.local.rank());           // choose sub−stream no. rank out of size streams
            r_insert.split(config.local.size() , config.local.rank());         // choose sub−stream no. rank out of size streams
            r_accept.split(config.local.size() , config.local.rank());         // choose sub−stream no. rank out of size streams
            r_shift.split(config.local.size() , config.local.rank());         // choose sub−stream no. rank out of size streams
            r_timep.split(config.local.size() , config.local.rank());           // choose sub−stream no. rank out of size streams
        }

    RealT StrongCoupling::tryInc(const RealT t, const RealT tp, const unsigned char f, const RealT fac, const RealT zetap)
    {
        const unsigned int n = M[f].rows();
        RealT Spi = hybCall(t, tp, f);
        if(n > 0)
        {
            MatrixT     P = M[f];
            RowVectorT  R(n);
            VectorT     Q(n);
            for(int i = 0; i < n; i ++)
            {
                auto sec = segmentCache[f].at(i);
                R(i) = hybCall(sec.first, tp, f);                                     //M_ij = Delta(tp_i - t_j) 
                Q(i) = hybCall(t, sec.second,  f);
            }
            Spi = Spi - R*P*Q;
            VLOG(5) << "Acceptance rate for insertion: " << fac*Spi;
            if(fac*Spi <= 0) LOG(WARNING) << "negative acceptance rate for insertion. A = " << fac*Spi;
            if(zetap < 0 || fac*Spi > zetap)
            {
                int index = segments.insertSegment(t, tp, f);
                if(index < 0)
                {
                    LOG(WARNING) << "trying to insert invalid segment, this is a bug.";
                }
                else
                {
                    MInc(&(M[f]), R, Q, 1./Spi, index);
                    segmentCache[f].emplace_back(t,tp);
                    VLOG(5) << "Inserted segment of flavor " << f << " at order " << n;
                }
            }
        }
        else
        {
            VLOG(5) << "Acceptance rate for insertion: " << fac*Spi;
            if(zetap < 0 || fac*Spi > zetap)
            {
                M[f].conservativeResize(1,1);
                M[f](1,1) = 1./Spi;
                segmentCache[f].emplace_back(t,tp);
                VLOG(5) << "Inserted segment of flavor " << f << " at order 0";
            }
        }
        return Spi;
    }

    //TODO: insert at correct row/column?
    RealT StrongCoupling::MInc(MatrixT *A, const RowVectorT &R, const VectorT &Q, const RealT Sp, const int index)
    {
        const unsigned int n = A->rows();
        MatrixT P = *A;
        (*A).conservativeResize(n+1, n+1);
        (*A)(n,n)                = Sp;
        A->topRightCorner(n,1)   = -(P*Q)*Sp;
        A->bottomLeftCorner(1,n) = -(R*P)*Sp;
        A->topLeftCorner(n,n)   += (A->topRightCorner(n,1))*(A->bottomLeftCorner(1,n))/Sp;
        return 1./Sp; 
    }

    RealT StrongCoupling::tryDec(const unsigned int row, const unsigned char f, const RealT fac, const RealT zetap)
    {
        RealT Spi = M[f](row,row);
        VLOG(5) << "Acceptance rate for removal: " << fac*Spi;
        if(fac*Spi <= 0) LOG(WARNING) << "negative acceptance rate for removal. A = " << fac*Spi;
        if(fac*Spi > zetap)
        {
            MDec(&(M[f]), row);
        segmentCache[f].erase(segmentCache[f].begin() + row);
        }
        return fac*Spi;
    }

    void StrongCoupling::MDec(MatrixT* A, const unsigned int row)
    {
        const unsigned int n = A->rows();
        if(n > 1)
        {
            (*A).topLeftCorner(n-1,n-1) = (A->topLeftCorner(n-1,n-1) - A->topRightCorner(n-1,1)\
                                * A->bottomLeftCorner(1,n-1) / (*A)(n-1,n-1)).eval();	//M = P-Q*R/S (8.45)
        }
        A->conservativeResize(n-1,n-1);
        VLOG(5) << "Deleted segment at order " << n+1;
    }

    bool StrongCoupling::tryInsAntiSeg(RealT t_n, RealT tp_n, const int f_n, const RealT fac, const RealT zetap)
    {
        int index = 0;
        std::pair<RealT,RealT> old = segmentCache[f_n][0];
        for(auto el : segmentCache[f_n])
        {
            if(segments.inSegment(t_n, el.first, el.second))
            {
                auto old = el;
                break;
            }
            index += 1;
        }
        RealT tf = old.second;
        if(t_n < old.first)
        {
            t_n += conf.beta;
            tp_n = std::fmod(tp_n, conf.beta);
            tf = std::fmod(tf, conf.beta);
        }
        if(tf < tp_n)
            tf += conf.beta;
        VLOG(5) << "Trying to insert [" << old.first << "] == [" << t_n << "] -- [" << tp_n << "] == [" << old.second << "]";
        RealT fac1    = tryDec(index, f_n, fac, -1.);
        RealT fac2    = tryInc(old.first, t_n, f_n, fac1, -1.);
        RealT fac3    = tryInc(tp_n, old.second, f_n, fac2, -1);
        if(fac3 < zetap)   // reset
        {
            tryRemAntiSeg(index, f_n, 0., -1.);
        }
        else
        {
            return true;
        }
    }

    bool StrongCoupling::tryRemAntiSeg(int index, const int f_n, const RealT fac, const RealT zetap)
    {
        LOG(ERROR) << "not implemented yet.";
        return false;
    }

    StrongCoupling::~StrongCoupling()
    {
            //ptr->free_sprng();                                                    //only needed for default interface
    }

    void StrongCoupling::doMeasurement(void)
    {
        //(M.transpose() * M2).diagonal(); // m[i,:].dot(m[:,i])
    }

    RealT StrongCoupling::acceptanceR(const RealT U,const RealT beta) const
    {
            return 0.0;
    }

    int StrongCoupling::update(const unsigned long int iterations)
    {
        for(unsigned long it=0l;it < iterations;it++)
        {
            //TODO: if totalSign > threshold: break;
            steps 	+= 1;
            const RealT f_n     = static_cast<int>(u(r_spin) + 0.5);    	    	// auxiliary spin variable (-1 or 1)
            const RealT zeta 	= u(r_insert);
            const RealT zetap 	= u(r_accept);
            const int mSize     = M[f_n].rows();
            VLOG(2) << "Updating Configuration, current size: " << mSize << ", zetap: " << zetap;
            if(zeta < 0.25)                     // insert segment
            {
                RealT t_n  = conf.beta*u(r_time);
                const RealT ml      = segments.maxl(t_n, f_n);
                RealT tp_n = ml*u(r_time);
                if(ml > 0)
                {
                    // force minimum size of segment ( tp_n \in [0, max_len) ) 
                    if(tp_n == t_n) tp_n += std::numeric_limits<RealT>::epsilon(); 
                    const RealT d_ov = segments.overlap(t_n, tp_n, f_n);
                    const RealT Sp = std::exp(conf.mu*(tp_n - t_n) - conf.U * d_ov)*conf.beta*ml/(mSize + 1);
                    tryInc(t_n, tp_n, f_n, Sp, zetap);
                }
            }
            else if(zeta >= 0.25 < 0.5)         // insert anti segment
            {
                RealT t_n       = conf.beta*u(r_time);
                const RealT ml  = -segments.maxl(t_n, f_n);
                RealT tp_n      = ml*u(r_time);
                if(ml > 0)
                {  
                    // force minimum size of segment ( tp_n \in [0, max_len) ) 
                    if(tp_n == t_n) tp_n += std::numeric_limits<RealT>::epsilon(); 
                    const RealT d_ov= segments.overlap(t_n, tp_n, f_n);
                    const RealT Sp = std::exp(conf.mu*(t_n - tp_n) + conf.U * d_ov)*conf.beta*ml/(mSize + 1);
                    auto res = tryInsAntiSeg(t_n, tp_n, f_n, Sp, zetap);
                    //TODO: finish
                }
            }
            else if(zeta >= 0.5 < 0.75)         // remove segment
            {
                if(M[f_n].rows() > 0)
                {
                    RealT n_n       = (int)(M[f_n].rows()*u(r_time));
                    auto seg        = segments.getTimeOrdered(n_n, f_n);
                    RealT t_n       = seg.first;
                    RealT tp_n      = seg.second;
                    const RealT d_ov= segments.overlap(t_n, tp_n, f_n);
                    const RealT ml  = segments.getTimeOrdered(n_n+1, f_n).first - t_n; 
                    if(ml <= 0) LOG(WARNING) << "Negative maximum length for laready inserted segment. This is a bug.";
                    const RealT Sp = std::exp(conf.mu*(t_n-tp_n) + conf.U*d_ov)*M[f_n].rows()/(ml*conf.beta);
                    tryDec(n_n, f_n, Sp, zetap);
                }
            }
            else                                // remove anti segment
            {
                if(M[f_n].rows() > 0)
                {

                    //TODO: finish
                }
            }
        }
        return 0;	
    }

} //end namespace DMFT
