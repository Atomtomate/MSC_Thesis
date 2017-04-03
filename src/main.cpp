#include "main.h"
INITIALIZE_EASYLOGGINGPP


//static const char* help = "Placeholder help string, TODO: set help string\n";
// TODO: shift mu by - 1/2 U
// TODO: switch from flat array to hash (template design)
int main(int argc,char **argv)
{
    MPI::Init(argc, argv); 
    //boost::mpi::environment env(argc, argv);
    //boost::mpi::communicator world;
    // ========== configure logging		============
    START_EASYLOGGINGPP(argc, argv);
    el::Configurations conf("debug_log.conf");
    el::Loggers::reconfigureAllLoggers(conf);
	std::cout << "give U, beta: ";
	double U = 2.5;
        double beta = 64.0;
        double mixing = 0;
        //int lattice,sc;
	std::cin >> U;
	std::cin >> beta;
	//std::cin >> mu;
        //std::cout << "mixing [0,1), 0 = no mixing: ";
        //std::cin >> mixing;
        //std::cout << "lattice: 0 for bethe, 1 for SC: ";
        //std::cin >> lattice;
        double mu = U/2.0;
        LOG(INFO) <<"Parameters: \n - beta:\t" <<  beta << "\n - U:\t" << U << "\n - mu:\t" << mu; 
        DMFT::Config config(beta, mu, U, DMFT::_CONFIG_maxMatsFreq);

        //DMFT::examples::_test_hysteresis();
	DMFT::examples::_test_SOH(config,true,mixing);

        //Config conf(beta,mu,U,_CONFIG_maxMatsFreq);
        /*if (!sc)
        {
            LOG(INFO) << "using DMFT SC condition for the bethe lattice";
	    DMFT::examples::_test_SOH(config,true, mixing);
        }
        else
        {
            LOG(INFO) << "using full DMFT loop";
	    DMFT::examples::_test_SOH(config,false, mixing);
        }*/
    MPI::Finalize();
    return 0;
}
