#include "StrongCoupling.hpp"
#include <iostream>

#define _CONFIG_USE_LEGENDRE_P 0
#define _CONFIG_HYB_DEBUG 0
#define _CONFIG_HYB_USE_MMEAS 0
namespace  DMFT
{
    StrongCoupling::StrongCoupling(GreensFct* const hyb, GreensFct* const gImp, Config const * const config, const unsigned int burninSteps):
    hyb(hyb), gImp(gImp), conf(config), burninSteps(burninSteps), segments(config), gImpLPoly(gImp, config), expOrd{config}
    {
        //TODO: check for correct initialization of totN and totalSign
        steps 		= 0;
        for(int f =0; f < _CONFIG_spins;f++)
        {
            lastSign	= 1;
            totalSign = 0;
            gl_c[f].fill(0);
            mfBins[f].fill(ComplexT(0.,0.));
            runningDet[f] = 0;
        }
        r_time.split(6, 0);
        r_spin.split(6, 1);
        r_insert.split(6, 2);
        r_accept.split(6, 3);
        r_shift.split(6, 4);
        r_timep.split(6, 5);
        r_time.split(config->local.size() , config->local.rank());           // choose sub−stream no. rank out of size streams
        r_spin.split(config->local.size() , config->local.rank());           // choose sub−stream no. rank out of size streams
        r_insert.split(config->local.size() , config->local.rank());         // choose sub−stream no. rank out of size streams
        r_accept.split(config->local.size() , config->local.rank());         // choose sub−stream no. rank out of size streams
        r_shift.split(config->local.size() , config->local.rank());         // choose sub−stream no. rank out of size streams
        r_timep.split(config->local.size() , config->local.rank());           // choose sub−stream no. rank out of size streams
    }

    StrongCoupling::~StrongCoupling()
    {
    }

    RealT StrongCoupling::acceptanceR(const RealT U,const RealT beta) const
    {
            return 0.0;
    }

    void StrongCoupling::swapRows(MatrixT *A, int from, int to)
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

    void StrongCoupling::MInc(MatrixT *A, const RowVectorT &R, const VectorT &Q, const RealT Sp, const int index)
    {
        const unsigned int n = A->rows();
        MatrixT P = *A;
        (*A).conservativeResize(n+1, n+1);
        (*A)(n,n)                = Sp;
        A->topRightCorner(n,1)   = -(P*Q)*Sp;
        A->bottomLeftCorner(1,n) = -(R*P)*Sp;
        A->topLeftCorner(n,n)    = P + A->topRightCorner(n,1) * A->bottomLeftCorner(1,n)/Sp;
        swapRows(A, n, index);
    } 

    void StrongCoupling::MDec(MatrixT* A, const unsigned int index)
    {
        const unsigned int n = A->rows();
        // swap rows to obtain the correct position
        if(n > 1)
        {
            swapRows(A, index, n-1);
            A->topLeftCorner(n-1,n-1) = (A->topLeftCorner(n-1,n-1) - A->topRightCorner(n-1,1)*A->bottomLeftCorner(1,n-1)  / (*A)(n-1,n-1)).eval();
        }
        A->conservativeResize(n-1,n-1);
    } 

    decltype(auto) StrongCoupling::tryZeroToFull(const unsigned f, const RealT fac, const RealT zetap)
    {
        if(!segments.isEmpty(f)) std::make_pair(fac, false);
        bool success = false;
        if(fac > zetap)
        {
            auto res = segments.insertFullLine(f); success = true; 
            if(res)
            {
                M[f].resize(1,1);
                M[f](0,0) = 1./hybCall(0., std::nextafter(conf->beta, 0), f);
            }else{
                LOG(ERROR) << "Could not insert full line from zero line!";
                throw std::logic_error("Segment insertion failed");
            }
            if(_CONFIG_HYB_DEBUG){
                runningDet[f] = 1./M[f](0,0);
            }
        }
        return std::make_pair(fac, success);
    }

    decltype(auto) StrongCoupling::tryFullToZero(const unsigned f, const RealT fac, const RealT zetap)
    {
        if(!segments.hasFullLine(f)) std::make_pair(fac, false);
        bool success = false;
        if(fac > zetap)
        {
            M[f].resize(0,0);
            segments.deleteSegment(0, f); success = true; 
            if(_CONFIG_HYB_DEBUG){
                runningDet[f] = 0;
            }
        }
        return std::make_pair(fac, success);
    }

    decltype(auto) StrongCoupling::tryInc(const RealT t, const RealT tp, const unsigned int f, const RealT fac, const RealT zetap)
    {
        bool success = false;
        if(tp <= t) LOG(WARNING) << "Trying to insert invalid segment.";
        RealT sign = tp > conf->beta ? -1 : 1;
        const unsigned int n = M[f].rows();
        RealT Spi = hybCall(t, tp, f);
        if(n > 0)
        {
            MatrixT     P = M[f];
            RowVectorT  R(n);
            VectorT     Q(n);
            for(int i = 0; i < n; i ++)
            {
                //TODO:TO  auto sec = segmentCache[f].at(i);
                auto sec = segments.getTimeOrdered(i,f);
                R(i) = hybCall(sec.first, tp, f);                                     //M_ij = Delta(tp_i - t_j) 
                Q(i) = hybCall(t, sec.second, f);
            }
            Spi = Spi - R*P*Q;
            //VLOG(5) << "Spi = " << Spi << " = " << hybCall(t, tp, f) << " - " << R << " * " << P << " * " << Q;
            if(zetap > 0) VLOG(5) << "Acceptance rate for insertion at n > 0: " << sign*fac*Spi << " = " << sign << " * " << fac << " * " << Spi  << " <?> zetap: " << zetap;
            if(zetap < 0 || std::abs(fac*Spi) > zetap)
            {
                int index = segments.insertSegment(t, tp, f);
                if(zetap > 0) VLOG(5) << "insert result: " << index;
                if(index < 0) LOG(WARNING) << "trying to insert invalid segment, this is a bug.";
                else
                {
                    MInc(&(M[f]), R, Q, 1./Spi, index);
                    //segmentCache[f].emplace_back(t,tp);
                    if(zetap > 0) VLOG(5) << "\033[32mInserted\033[0m segment of flavor " << f << " at order " << n;
                    success = true;
                    if(_CONFIG_HYB_DEBUG){
                        runningDet[f] = runningDet[f]*Spi;
                    }
                }
            }
        }
        else
        {
            if(zetap > 0) VLOG(5) << "Acceptance rate for insertion at n = 0: " << sign*fac*Spi << " = " << sign << " * " << fac << " * " << Spi << " <?> zetap: " << zetap;
            if(zetap < 0 || std::abs(fac*Spi) > zetap)
            {
                M[f].resize(1,1);
                M[f](0,0) = 1./Spi;
                int index = segments.insertSegment(t, tp, f);
                if(zetap > 0) VLOG(5) << "\033[32mInserted\033[0m segment of flavor " << f << " at order 0, with index " << index;
                success = true;
                if(_CONFIG_HYB_DEBUG){
                    runningDet[f] = Spi;
                }
            }
        }
        return std::make_pair(fac*Spi*sign, success);
    }

    decltype(auto) StrongCoupling::tryDec(const unsigned int row, const unsigned int f, const RealT fac, const RealT zetap)
    { 
        bool success = false;
        RealT Sp = M[f](row,row);
        auto seg= segments.getTimeOrdered(row, f);
        RealT sign = seg.second > conf->beta ? -1 : 1;
        if(zetap > 0) VLOG(5) << "Acceptance rate for removal: " << sign*fac*Sp << " = " << sign << " * " << fac << " * " << Sp  << "" " <?> zetap: " << zetap;
        if(zetap < 0 || std::abs(fac*Sp) > zetap)
        {
            MDec(&(M[f]), row);
            auto seg = segments.deleteSegment(row, f);
            sign = seg.second > conf->beta ? -1 : 1;
            if(zetap > 0) VLOG(5) << "\033[31mRemoved\033[0m segment of flavor " << f << " at order 0, row " << row << ". sign: " << sign;
            if(seg.first < 0) LOG(WARNING) << "Segment romoval failed in tryDec!";
            success = true;
            if(_CONFIG_HYB_DEBUG){
                runningDet[f] = runningDet[f]*Sp;
                if(M[f].rows() == 0) runningDet[f] = 0.;
            }
        }
        return std::make_pair(sign*fac*Sp, success);
    }

    decltype(auto) StrongCoupling::tryInsAntiSeg(RealT t_n, RealT tp_n, const unsigned int f, const RealT fac, const RealT zetap)
    {
        bool success = false;
        RealT fac3 = 0.;
        const auto Mbak = M[f].eval();
        const auto sbak = segments;
        const auto rdbak = runningDet[f];
        if(segments.hasFullLine(f))
        {
            RealT fac_d = tryDec(0, f, 1., -1.).first;
            RealT ts = tp_n;
            RealT tf = t_n + conf->beta;                                         // case: [=====] to [=== --- ===]
            if(tp_n > conf->beta)                                                // case: [=====] to [--- === ---]
            {
                ts = tp_n-conf->beta;
                tf = t_n;
            }
            RealT fac_i = tryInc(ts, tf, f, 1., -1).first; 
            fac3 = fac*fac_d*fac_i;
            VLOG(5) << "  Acceptance rate: " << fac3 << " = " << fac << " * " << fac_d << " * " << fac_i << " <?> " << zetap;
        }
        else
        {
            int index = segments.getIndex(t_n, f);
            auto old = segments.getTimeOrdered(index, f);
            //VLOG(7) << "\t old segment for removal: [" << old.first << " == ] " << t_n  << " -- " << tp_n << " [ == " << old.second << "]. index old: " << index; 
            auto ptl    = segments.timeOrder({std::make_pair(old.first, t_n), std::make_pair(tp_n, old.second)});
            VLOG(5) << "Trying to insert [" << ptl[0].first << "] == [" << ptl[0].second << "] -- [" << ptl[1].first << "] == [" 
                    << ptl[1].second << "], at index " << index << ", fac: " << fac;
            RealT fac_i  = tryDec(index, f, 1., -1.).first;
            RealT fac_d1 = tryInc(ptl[0].first, ptl[0].second, f, 1., -1.).first;
            RealT fac_d2 = tryInc(ptl[1].first, ptl[1].second, f, 1., -1.).first;
            fac3 = fac*fac_i*fac_d1*fac_d2;
            if(zetap > 0) VLOG(5) << "Acceptance rate for anti segment insertion: " << fac3 << " = " << fac << " * " <<  fac_i << " * " << fac_d1 <<  " * " << fac_d2 << " <?> zetap: " << zetap;
        }
        if(std::abs(fac3) < zetap)   // reset
        {
            M[f] = Mbak;
            segments = sbak;
            runningDet[f] = rdbak;
        }
        else
        {
            if(zetap > 0) VLOG(5) << "\033[32mInserted\033[0m anti segment of flavor " << f << " at order " << M[f].rows()-1;
            success = true;
        }
        return std::make_pair(fac3, success);
    }

    decltype(auto) StrongCoupling::tryRemAntiSeg(int index, const unsigned int f, const RealT fac, const RealT zetap)
    {
        bool success = false;
        RealT fac3 = 0.;
        const auto Mbak = M[f].eval();
        const auto sbak = segments;
        const auto rdbak = runningDet[f];
        if(false) //M[f].rows() == 0)
        {
            const RealT ts = 0.;
            const RealT tf = std::nextafter(conf->beta, 0);
            RealT Spi = hybCall(ts, tf, f);
            if(zetap > 0) VLOG(5) << "Acceptance rate for anti segment removal at order 0: " << fac*Spi << " = " << fac << " * " <<  Spi <<  " <?> zetap: " << zetap;
            fac3 = fac*Spi;
            if(std::abs(fac3) > zetap)
            {
                tryInc(ts, tf, f, 1, -1);
                if(segments.insertFullLine(f))
                {
                    M[f].resize(1,1);
                    M[f](0,0) = 1./Spi;
                    if(zetap > 0) VLOG(5) << "\033[32mRemoved\033[0m anti segment of flavor " << f << " at order 0, row " << index;
                    success = true;
                }
                else
                    LOG(WARNING) << "Insertion of full line after anti segment removal failed!";
            }
        }
        else if(M[f].rows() == 1 && !segments.hasFullLine(f))  // one segment to full line
        {
            // fac * [decrease 1] * [increase 1]
            fac3 = fac*M[f](0,0)*hybCall(0., std::nextafter(conf->beta, 0), f);
            if(zetap > 0) VLOG(5) << "Acceptance rate for anti segment removal at order 1: " << fac3 << " = " << fac << " * " <<  M[f](0,0) << " * " <<  hybCall(0., std::nextafter(conf->beta, 0), f)<<  " <?> zetap: " << zetap;
            if(std::abs(fac3) > zetap)
            {
                M[f](0,0) = 1./hybCall(0., std::nextafter(conf->beta, 0), f);
                segments.deleteSegment(0, f);
                int index = segments.insertFullLine(f);
                if(zetap > 0) VLOG(5) << "\033[31mRemoved\033[0m anti segment of flavor " << f << " at order 1, row " << index;
                runningDet[f] = 1./M[f](0,0);
                success = true;
            }
        }
        else if(M[f].rows() > 1)
        {
            auto oldS1 = segments.getTimeOrdered(index, f);
            auto oldS2 = segments.getTimeOrdered(index + 1, f);
            RealT tf = oldS2.second < oldS1.first ? oldS2.second + conf->beta : oldS2.second;
            const RealT fac_d1 = tryDec(index, f, 1., -1).first;
            int index2 = (index) % M[f].rows();
            const RealT fac_d2 = tryDec(index2, f, 1., -1.).first;
            const RealT fac_i1 = tryInc(oldS1.first, tf, f, 1., -1.).first;
            fac3 = fac*fac_d1*fac_d2*fac_i1;
            if(zetap > 0) VLOG(5) << "Acceptance rate for anti segment removal: " << fac3 << " = " << fac << " * " <<  fac_d1 << " * " << fac_d2 << " * " << fac_i1 <<  " <?> zetap: " << zetap;
            //if(fac3 < 0) LOG(WARNING) << "negative acceptance rate for anti segment removal: " << fac3;
            if(std::abs(fac3) < zetap)   // reset
            {
                M[f] = Mbak;
                segments = sbak;
                runningDet[f] = rdbak;
            }
            else
            {
                if(zetap > 0) VLOG(5) << "\033[32mInserted\033[0m anti segment of flavor " << f << " at order " << M[f].rows()-1;
                success = true;
            }
        }
        return std::make_pair(fac3, success);
    } 

    void StrongCoupling::updateContribution(void)
    {
        if(steps < burninSteps) return;	            // skip accumulation while still in burn in period
        totalSign += lastSign;
        if(_CONFIG_USE_LEGENDRE_P){
            updateLPoly();
            return;
        }
        for(int f = 0; f < _CONFIG_spins; f++)
        {
            unsigned int n = M[f].rows() - segments.hasFullLine(f);
            expOrd(n, f);
        }
        //TODO: tmp
        return;
        for(int f = 0; f < _CONFIG_spins; f++)
        {
            unsigned int n = M[f].rows() - segments.hasFullLine(f);
            //continue;
            //TODO: make legendre basis optional
            if(segments.hasFullLine(f))
            {
                RealT tp = conf->beta;
                itBins[f].at(0)(-std::abs(M[f](0,0)));
                return ;
            }
            for(int i=0; i<n; i++)
            {
                for(int j=0; j<n; j++)
                {
                    const RealT ts = segments.getTimeOrdered(i, f).first; 
                    const RealT tf = segments.getTimeOrdered(j, f).second;
                    RealT tp = std::fmod(tf, conf->beta) - ts;
                    if(ts > conf->beta) LOG(WARNING) << "start time can not wrap arround, error!";
                    RealT s = lastSign;
                    if(tp < 0)
                    {
                        s *= -1;
                        tp += conf->beta;
                    }
                    //RealT tp = std::fmod(segments.getTimeOrdered(i, f).second, conf->beta) - std::fmod(segments.getTimeOrdered(j, f).first, conf->beta);
                    unsigned int index = static_cast<unsigned int>( _CONFIG_maxTBins*tp/conf->beta);
                    //if(index == _CONFIG_maxTBins - 1) index = 0;
                    if(_CONFIG_HYB_USE_MMEAS)
                    {
                        for(int n =0; n < _CONFIG_maxMatsFreq; n++)
                        {
                            mfBins[f][n] -= std::exp(ComplexT(0., mFreqS(n,conf->beta)*tp))*std::abs(M[f](i,j));
                        }
                    }else{
                        if(index == _CONFIG_maxTBins || index < 0) LOG(WARNING) << "invalid accumulation index";
                        itBins[f].at(index)(-std::abs(M[f](i,j)));//
                    }
                }
            }
        }
    }

    void StrongCoupling::updateLPoly(void)
    {
        for(int f = 0; f < _CONFIG_spins; f++)
        {
            if(segments.hasFullLine(f))
            {
                RealT x = 1.;
                for(int l=0; l<_CONFIG_maxLPoly;l++)
                {
                    gl_c[f][l] -= std::sqrt(2.*l+1.)*boost::math::legendre_p(l,x)*(std::abs(M[f](0,0)));
                }
                return;
            }
            else
            {
                int n = M[f].rows();
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<n; j++)
                    {
                        const RealT ts = segments.getTimeOrdered(i, f).first; 
                        const RealT tf = segments.getTimeOrdered(j, f).second;
                        RealT tp = ts - std::fmod(tf, conf->beta);
                        if(ts > conf->beta) LOG(WARNING) << "start time can not wrap arround, error!";
                        RealT s = lastSign;
                        if(tp < 0)
                        {
                            s *= -1;
                            tp += conf->beta;
                        }
                        RealT x = (2.*tp/conf->beta) - 1. ;
                        //LOG(ERROR) << ts << " - " << tf << " = " << tp <<" : " << conf->beta << "; " << 2*tp/conf->beta;
                        for(int l=0; l<_CONFIG_maxLPoly;l++){
                            gl_c[f][l] -= std::sqrt(2.*l+1.)*boost::math::legendre_p(l,x)*(std::abs(M[f](j,i)));
                        }
                    }
                }
            }
        }
    }

    void StrongCoupling::computeImpGF(void)
    {
        expOrd.writeResults();
        const long itCount = gImp->getItGF().rows();
        ImTG g_it = ImTG::Zero(_CONFIG_maxTBins, _CONFIG_spins);
        MatG g_mf = MatG::Zero(_CONFIG_maxMatsFreq, _CONFIG_spins);
        for(int f=0; f<_CONFIG_spins; f++)
        {
            
            for(int l = 0; l<_CONFIG_maxLPoly; l++)
            {
                gImpLPoly.setGl(l, f, gl_c[f][l]/(conf->beta*(steps-burninSteps)));
            }
            for(int t=0; t<_CONFIG_maxTBins; t++)
            {
                g_it(t, f) = boost::accumulators::sum(itBins[f].at(t));///(conf->beta*totalSign[f]);
            }
            for(int n=0; n<_CONFIG_maxMatsFreq; n++)
            {
                g_mf(n, f) = mfBins[f].at(n);
            }
        }
        
        LOG(DEBUG) << "steps: " << steps-burninSteps << ", vs total sign: " << totalSign;
        if(_CONFIG_USE_LEGENDRE_P){
            gImpLPoly.setGF();
        } else {
            if(_CONFIG_HYB_USE_MMEAS){
                g_mf = g_mf/(conf->beta*(steps-burninSteps));
                gImp->setByMFreq(g_mf);
                gImp->transformMtoT();
            }else{
                g_it = g_it/(conf->beta*(steps-burninSteps));
                gImp->setByT(g_it);
                gImp->transformTtoM();
            }
        }
        
        //gImp->markTSet();
        //gImp->setByMFreq(g_mf);
        //gImp->transformMtoT();
        //for(int f =0; f <_CONFIG_spins; f++)
        //    LOG(DEBUG) << boost::accumulators::mean(expOrd[f]) << " +- " <<  boost::accumulators::variance(expOrd[f]) << ", " << boost::accumulators::skewness(expOrd[f]) << ", " <<boost::accumulators::kurtosis(expOrd[f])  <<  "\n";
    }

    void StrongCoupling::rebuildM(const bool timeOrdered)
    {
        for(int f = 0; f < _CONFIG_spins; f++)
        {
            MatrixT M_new;
            const int n = M[f].rows();
            M_new.resize(n,n);
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j < n; j++)
                    M_new(i,j) = hybCall(segments.getTimeOrdered(j,f).first, segments.getTimeOrdered(i,f).second, f);
            }
            Minv[f] = M_new;
            M[f] = M_new.inverse();
        }
    }
//TODO: M operations to separate class
//TODO: local weight calculation to separate function
    int StrongCoupling::update(const unsigned long int iterations)
    {
        ProposalRes res = std::make_pair(0, false);
        for(unsigned long it=0l;it < iterations;it++)
        {
            //TODO: if totalSign > threshold: break;
            steps 	+= 1;
            const RealT f_n   = static_cast<int>(u(r_spin) + 0.5);
            const RealT zeta  = u(r_insert);
            const RealT zetap = u(r_accept);
            const int mSize   = M[f_n].rows() - segments.hasFullLine(f_n);
            VLOG(5) << "sign: " << lastSign << ", zeta = " << zeta << "\n" << segments.print_segments();
            if(zeta < 0.1)
            {
                VLOG(3) << "Propose \033[36mzero mode\033[0m \033[33mswap\033[0m, with line: " << f_n << " is full = " << segments.hasFullLine(f_n) << ", is empty = " << segments.isEmpty(f_n);
                const RealT d_ov = segments.overlap(0., std::nextafter(conf->beta, 0), f_n);
                const RealT Sp   = std::exp(conf->mu*conf->beta - conf->U * d_ov);
                if(segments.hasFullLine(f_n))  res = tryFullToZero(f_n, 1./Sp, zetap);
                else if(segments.isEmpty(f_n)) res = tryZeroToFull(f_n, Sp, zetap);
            }
            else if(zeta < 0.35)                                                // insert segment
            {
                RealT t_n           = conf->beta*u(r_time);
                const RealT ml      = segments.maxl(t_n, f_n);
                if(ml > 0)
                {
                    RealT tp_n            = ml*u(r_time) + t_n;
                    if(tp_n == t_n) tp_n  = std::nextafter(tp_n, 2*conf->beta);      // force minimum size of segment ( tp_n \in [0, max_len) ) 
                    const RealT d_ov      = segments.overlap(t_n, tp_n, f_n);
                    VLOG(3) << "Propose \033[34msegment\033[0m \033[33minsertion\033[0m, with max length: " << ml << " pre factor: " << std::exp(conf->mu*(tp_n - t_n) - conf->U * d_ov) <<": " << f_n << " [" << t_n << ", "  << tp_n << "]";
                    const RealT Sp = std::exp(conf->mu*(tp_n - t_n) - conf->U * d_ov)*conf->beta*ml/(mSize + 1);
                    res = tryInc(t_n, tp_n, f_n, Sp, zetap);
                }
            } else if(zeta < 0.55){                                             // insert anti-segment
                RealT t_n           = conf->beta*u(r_time);
                const RealT ml      = -segments.maxl(t_n, f_n);
                if(ml > 0)
                {
                    RealT tp_n            = ml*u(r_time) + t_n;
                    if(tp_n == t_n) tp_n  = std::nextafter(tp_n, 2*conf->beta);      // force minimum size of segment ( tp_n \in [0, max_len) ) 
                    VLOG(3) << "Propose \033[35m anti segment\033[0m \033[33minsertion\033[0m, with max length: " << ml << ": " << f_n << " [" << t_n << ", "  << tp_n << "]";
                    const RealT d_ov      = segments.overlap(t_n, tp_n, f_n);
                    const RealT Sp = std::exp(conf->mu*(t_n - tp_n) + conf->U * d_ov)*conf->beta*ml/(mSize + 1);
                    res = tryInsAntiSeg(t_n, tp_n, f_n, Sp, zetap);
                }
            }
            else if(zeta < 0.80)                                                // remove segment
            {
                if(mSize > 0)
                {
                    const int n_n  = (int)(mSize)*u(r_time);
                    auto seg       = segments.getTimeOrdered(n_n, f_n);
                    auto seg2      = segments.getTimeOrdered(n_n + 1, f_n);
                    const RealT d_ov= segments.overlap(seg.first, seg.second, f_n);
                    RealT Sp = std::exp(conf->mu*(seg.first - seg.second) + conf->U*d_ov);
                    RealT ml = seg2.second - seg.first;
                    VLOG(3) << "Propose \033[34msegment\033[0m \033[33mremoval\033[0m, with max length: " << ml << ": " << f_n << " [" << seg.first << ", "  << seg.second << "]";
                    ml = ml < 0 ? ml + conf->beta : ml;
                    if(!segments.hasFullLine(f_n))
                        Sp *= mSize/(ml*conf->beta);
                    res = tryDec(n_n, f_n, Sp, zetap);
                }
            }
            else                                                                // remove anti segment
            {
                if(mSize > 0)
                {
                    const int n_n = (int)(mSize)*u(r_time);
                    auto seg      = segments.getTimeOrdered(n_n, f_n);
                    auto seg2     = segments.getTimeOrdered(n_n + 1, f_n);
                    RealT t_n     = seg.second;
                    RealT tp_n    = seg2.first;
                    tp_n          = t_n > tp_n ? tp_n + conf->beta : tp_n;
                    RealT ml = segments.dist_to_next_end(std::nextafter(seg.second, 2*conf->beta), f_n);
                    const RealT d_ov= segments.overlap(t_n, tp_n, f_n);
                    VLOG(3) << "Propose \033[35manti segment\033[0m \033[33mremoval\033[0m, with max length: " << ml << ": " << f_n << " [" << seg.first << ", "  << seg.second << "]";
                    const RealT Sp  = std::exp(conf->mu*(tp_n - t_n) - conf->U * d_ov)*(mSize)/(conf->beta*ml);
                    VLOG(5) << "Sp: " << Sp << "= exp(" << conf->mu << " * " << tp_n - t_n << " - " << conf->U << " * " << d_ov << ") * " <<M[f_n].rows()  << " / " << conf->beta << " * " << ml;
                    res = tryRemAntiSeg(n_n, f_n, Sp, zetap);
                }
            }
            VLOG(5) << "success: " << res.second << ", ratio: " << res.first; 
            if(res.second)                          // if update got rejected, we continue to use the old sign
            {
                if(res.first < 0)    // otherwise, if the sign changed, update
                {
                    //LOG(DEBUG) << "signs do not match, potential sign problem";
                    lastSign *= -1; // (we can only track changes in sign)
                }
            }
            if(_CONFIG_HYB_DEBUG)
                debug_test(f_n);
            updateContribution();
        }
        return 0;	
    }

    void StrongCoupling::debug_test(const unsigned int f_n)
    {
        for(int f = 0; f< _CONFIG_spins; f++)
        {
            if(M[f].cols() != M[f].rows()) LOG(ERROR) << "M[" << f << "] not square!";
        }
        if(segments.size(f_n) != M[f_n].rows())
        {
            LOG(ERROR) << "Number of rows/cols in M does not match number of segments! Rows in M[" << f_n << "] " << M[f_n].rows() << ", in segments: " << segments.size(f_n);
            LOG(ERROR) << segments.print_segments();
            throw std::logic_error("Number of rows and cols in M do not match.");
        }
        auto Mbak = M[f_n];
        rebuildM(true);
        if(!M[f_n].isApprox(Mbak, 0.00001))
        {
            LOG(DEBUG) << "BEFORE rebuild \n" << Mbak << "\ninv\n" << Mbak.inverse();
            LOG(DEBUG) << "AFTER rebuild \n" << M[f_n] << "\ninv\n" << M[f_n].inverse();
            throw std::logic_error("M was computed incorrectly.");
            //exit(0);
        }
        std::array<RealT, _CONFIG_spins> perm_sign;
        for(int f = 0; f < _CONFIG_spins; f++)
        {
            fullDet[f] = Minv[f].determinant();
            if(Minv[f].rows() == 0) fullDet[f] = 0.;
            perm_sign[f] = segments.sign(f);
            if(std::abs(1. - runningDet[f]/fullDet[f] ) > 0.001)
            {
                LOG(ERROR) << "expected " << fullDet[f] << " for full det, but got " << runningDet[f] << " from running det for f = " << f;
                LOG(ERROR) << "M^-1\n" << Minv[f] << "\n M : \n" << M[f];
                throw std::logic_error("Determinant ratio incorrect");
            }
            //perm_sign[f]*
            if((2*(fullDet[f] >= 0) - 1.) != lastSign)
            {
                //LOG(ERROR) << segments.print_segments();
                //for(int i =0; i<M[f].rows();i++) LOG(INFO) << segments.getTimeOrdered(i,f).first << " to " << segments.getTimeOrdered(i,f).second << " =" <<  hybCall(segments.getTimeOrdered(i,f).first, segments.getTimeOrdered(i,f).second, f);
                //LOG(ERROR) << "M^-1\n" << Minv[f] << "\n M : \n" << M[f];
                LOG(ERROR) << "wrong sign for f = " << f << ". expected " << perm_sign[f]*(2*(fullDet[f] > 0) - 1.) << " = " << perm_sign[f] << " * " << (2*(fullDet[f] > 0) - 1.) << " but got " << lastSign;
                //LOG(ERROR) << "Full det: " << fullDet[f] << ", running Det: " << runningDet[f];
                lastSign *= -1;
                //throw std::logic_error("wrong sign");
            }
                
        }
    }

} //end namespace DMFT
