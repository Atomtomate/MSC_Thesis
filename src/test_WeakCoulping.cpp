#include "WeakCoupling.h"
#include <iostream>

WeakCoupling::WeakCoupling()
{
	auto seed = 1234567890;
	// auto seed = make_sprng_seed();
	// int gtype;
	init_sprng(seed,SPRNG_DEFAULT, DEFAULT_RNG_TYPE);	
	
	zeroShift = 0.1;							// TODO: get zeroShift from user?
	
}

WeakCoupling::~WeakCoupling()
{
	//ptr->free_sprng(); 							//only needed for default interface
}

double WeakCoupling::weissGreensFunct(const double tau,const int& sigma){
	return sprng();// TODO: placeholder	
}

double WeakCoupling::acceptanceR(MatrixT *M, MatrixT *U, MatrixT *G0, SConfig *config)
{
	// TODO: placeholder
	//DataStructs::MatrixT S = G0[n+1,n+1]
	//DataStructs::MatrixT Q = G0[0:n,n+1]
	//DataStructs::MatrixT R = G0[n+1,0:n]
	return 0.0;
}

int WeakCoupling::update(SConfigL& confs, MatrixT& M, MatrixT& G0, const double U,const double beta)
{
	DBGPRINT("\n\nUpdating Configuration\n");
	auto	rndReal = sprng();
	auto 	zetap 	= sprng();
	char	sigma 	= 1;
	int	n	= M.rows();//get expansion order	
	// TODO: compare speed with blocks fro mlarge matrix
	if (sprng() < 0.5 or n < 500) {							// try to insert spin
		DBGPRINT("Trying to add configuration\n");
		auto s_n 	= static_cast<int>(rndReal + 0.5);
		auto alpha 	= 0.5 + sigma*s_n*(0.5 + zeroShift);		// for 0 shift: \in {0,1}
		auto t_n  	= sprng()*beta;
		RowVectorT R(n);
		VectorT Q(n);
		for(int i=0; i < n; i++){
			Q(i) = weissGreensFunct(confs[i].first - t_n, sigma);
			R(i) = weissGreensFunct(t_n - confs[i].first, sigma);
		}
		double A = -beta*U/(n+1)*( M(n-1,n-1) - R*M*Q );			// TODO: product over sigma
		//std::printf("n: %d, confs_size %d, A: %f, zetap: %f\n",n,static_cast<int>(confs.size()),A,zetap);
		
		if (zetap < 0.5 or n < 500) { 			//compute A(n+1 <- n) (8.36) (8.39) - WeakCoupling::acceptanceR
			DBGPRINT("Accepted\n");
			pushConfig(confs,std::make_pair(t_n, s_n));
			MatrixT M_new(n+1,n+1);
			M_new(n,n) 		= 1.0/(M(n-1,n-1) * R*(M*Q));	// (S-R[M_n Q])^-1 	(8.40)
			M_new.topRightCorner(n,1) = -(M*Q)*M_new(n,n);		// -[M_n Q] Sp		(8.41)	
			M_new.bottomLeftCorner(1,n)= -M_new(n,n)*(R*M);		// -Sp*R*M		(8.42)
			M_new.topLeftCorner(n,n)= M + (M*Q) *M_new(n,n)* (R*M);	// M + M*Q * Sp * R*M	(8.43)
			M.resize(n+1,n+1);
			M = M_new;
		}
		return 1;							// inserted SConfig
	} else if (confs.size() > 0) {						// try to remove spin
		DBGPRINT("Trying to remove configuration\n");
		//std::printf("n: %d, confs_size %d, zetap: %f\n",n,static_cast<int>(confs.size()),zetap);
		int rndConfPos = static_cast<int>((10+1) * rndReal);
		if ( zetap < 0.5){//M(n-1,n-1) ) {			//compute A(n-1 <- n) (8.37) (8.44) - WeakCoupling::acceptanceR
			DBGPRINT("Accepted\n");
			MatrixT M_new(n-1,n-1);
			M_new = M.topLeftCorner(n-1,n-1) - M.block(0,n-1,n-1,1) * M.block(n-1,0,1,n-1) / M(n-1,n-1); //M = P-Q*R/S (8.45)
			M.topLeftCorner(n-1,n-1) = M_new;
			M.conservativeResize(n-1,n-1);		// TODO: test resize normal, then copy
			
			deleteConfig(confs,n);	
		}
		return -1;												// removed Sconfig
	}
	// return C,M
	return 0;
}


