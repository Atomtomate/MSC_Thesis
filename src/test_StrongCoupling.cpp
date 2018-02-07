void test_swap(void)
{
        auto P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic >(M[0].rows());
        P.setIdentity();
        P.applyTranspositionOnTheRight(0, 1);
        P.applyTranspositionOnTheRight(1, 2);
        LOG(DEBUG) << "\n" << M[0];
        LOG(DEBUG) << "\n" << P*M[0]*(P.transpose());
        swapRows(&(M[0]),2,0);
        LOG(DEBUG) << "\n" << M[0];
        swapRows(&(M[0]),0,2);
        LOG(DEBUG) << "\n" << M[0];
        swapRows(&(M[0]),2,1);
        LOG(DEBUG) << "\n" << M[0];

        ProposalRes res = std::make_pair(0, false);
        LOG(INFO) << "1.";
        LOG(DEBUG) << "insert 1";
        tryInc(0.1, 10., 0,1.,-1);
        LOG(DEBUG) << "insert 2";
        tryInc(16,17, 0, 1.,-1);
        LOG(DEBUG) << "insert 3";
        tryInc(18.1, 19., 0,1.,-1);
        LOG(DEBUG) << "insert 4";
        tryInc(11.1, 15., 0,1.,-1);
        LOG(DEBUG) << "insert 5";
        tryInc(19.1, 19.3, 0,1.,-1);
        LOG(DEBUG) << "remove 1";
        tryDec(3, 0, 1, -1.);
        auto Mbak2 = M;
        LOG(ERROR) << "1.1";
        rebuildM(false);
        LOG(INFO) << ":\n" << segments.print_segments();
        LOG(ERROR) << "1.2";
        if(!M[0].isApprox(Mbak2[0]) || !M[1].isApprox(Mbak2[1]))
        {
            LOG(DEBUG) << "BEFORE rebuild\n" << Mbak2[0] << "\ninv\n" << Mbak2[0].inverse();
            LOG(DEBUG) << "AFTER rebuild\n" << M[0] << "\ninv\n" << M[0].inverse();
            LOG(WARNING) << "\033[31merror in M\033[0m";
        }
        LOG(DEBUG) << "remove 2";
        tryDec(1, 0, 1, -1.);
        LOG(INFO) << ":\n" << segments.print_segments();
        Mbak2 = M;
        rebuildM(false);
        if(!M[0].isApprox(Mbak2[0]) || !M[1].isApprox(Mbak2[1]))
        {
            LOG(DEBUG) << "BEFORE rebuild\n" << Mbak2[0] << "\ninv\n" << Mbak2[0].inverse();
            LOG(DEBUG) << "AFTER rebuild\n" << M[0] << "\ninv\n" << M[0].inverse();
            LOG(WARNING) << "\033[31merror in M\033[0m";
        }
        LOG(INFO) << ":\n" << segments.print_segments();
        LOG(DEBUG) << "remove 3";
        tryDec(2, 0, 1, -1.);
        Mbak2 = M;
        rebuildM(false);
        if(!M[0].isApprox(Mbak2[0]) || !M[1].isApprox(Mbak2[1]))
        {
            LOG(DEBUG) << "BEFORE rebuild\n" << Mbak2[0] << "\ninv\n" << Mbak2[0].inverse();
            LOG(DEBUG) << "AFTER rebuild\n" << M[0] << "\ninv\n" << M[0].inverse();
            LOG(WARNING) << "\033[31merror in M\033[0m";
        }
}
