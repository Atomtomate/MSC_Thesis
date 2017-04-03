void updateTests(){
                        /* method1
                        MatrixT M_new(n+1,n+1);
                        M_new(n,n)      = Sp;                                                   // (S-R[M_n Q])^-1      (8.40)
                        M_new.topRightCorner(n,1).noalias() = -(M*Q)*Sp;                // -[M_n Q] Sp          (8.41)  
                        M_new.bottomLeftCorner(1,n).noalias() = -Sp*(R*M);              // -Sp*R*M              (8.42)
                        M_new.topLeftCorner(n,n).noalias() = M + (M*Q) *Sp* (R*M);      // M + M*Q * Sp * R*M   (8.43)
                        M = M_new;
                        */

                        /*//method 3
                        MatrixT M_new(n+1,n+1);
                        M_new = (Eigen::MatrixXd(n+1,n+1)<< M + (M*Q) *Sp* (R*M), -(M*Q)*Sp, -Sp*(R*M), Sp).finished();
                        M = M_new;
                        */


                        /*// method2 2nd fastest
                        M.topLeftCorner(n-1,n-1) -= M.block(0,n-1,n-1,1)*M.block(n-1,0,1,n-1)/M(n-1,n-1);       
                        M.conservativeResize(n-1,n-1);
                        */

                        //method1 - fastest
                        MatrixT M_new(n-1,n-1);
                        M_new.noalias() = M.topLeftCorner(n-1,n-1) - M.block(0,n-1,n-1,1) * M.block(n-1,0,1,n-1) / M(n-1,n-1); //M = P-Q*R/S (8.45)
                        M = M_new;

                        /* method3
                        MatrixT M_new(n-1,n-1);
                        M_new << M.topLeftCorner(n-1,n-1) - M.block(0,n-1,n-1,1) * M.block(n-1,0,1,n-1) / M(n-1,n-1);
                        M = M_new;
                        */
}
