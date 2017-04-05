#include "main.h"
INITIALIZE_EASYLOGGINGPP

void dataGenerationComm(boost::mpi::communicator local, boost::mpi::communicator world);
void dataCollectionComm(boost::mpi::communicator local, boost::mpi::communicator world);

//TODO: generate random number stream in collector processes?
//static const char* help = "Placeholder help string, TODO: set help string\n";
int main(int argc,char **argv)
{
    //TODO: centralize IO over MPI_rank 0
    START_EASYLOGGINGPP(argc, argv);
    el::Configurations conf("debug_log.conf");
    el::Loggers::reconfigureAllLoggers(conf);
    // MPI stuff
    MPI::Init(argc, argv); 
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    //TODO: for now only world.rank() == 0 collects data. change this after some testing
    bool isGenerator = (world.rank() > 0);
    boost::mpi::communicator local = world.split(isGenerator ? 0 : 1);


    double U = 2.5;
    double beta = 64.0;
    double mixing = 0;
    double mu = U/2.0;
    // ========== configure logging		============
    if(world.rank() == 0)
    {
        //std::cout << "give U, beta: ";
        //int lattice,sc;
        //std::cin >> U;
        //std::cin >> beta;
        //std::cin >> mu;
        //std::cout << "mixing [0,1), 0 = no mixing: ";
        //std::cin >> mixing;
        //std::cout << "lattice: 0 for bethe, 1 for SC: ";
        //std::cin >> lattice;
        LOG(INFO) <<"Parameters: \n - beta:\t" <<  beta << "\n - U:\t" << U << "\n - mu:\t" << mu; 
    }
    DMFT::Config config(beta, mu, U, DMFT::_CONFIG_maxMatsFreq, local, world, isGenerator);

    
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
