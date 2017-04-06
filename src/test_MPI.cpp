// =====================================================================================
//
//       Filename:  test_MPI.cpp
//
//    Description:  test cases for MPI communication of sample data. 
//                  compile with mpic++ -std=c++14 test_MPI.cpp -lm -lboost_mpi -lboost_serialization.
//                  important :
//                      - for loops with isend will override buffer
//                      - remember to close receiver before opening new one
//          
//
//        Version:  1.0
//        Created:  06.04.2017 11:12:15
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Julian Stobbe (jst), stobbe@itp.uni-frankfurt.de
//   Organization:  Goethe University Frankfurt
//
// =====================================================================================

// Copyright (C) 2006 Douglas Gregor <doug.gregor@gmail.com>

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// An example using Boost.MPI's split() operation on communicators to
// create separate data-generating processes and data-collecting
// processes.
#include <boost/mpi.hpp>
#include <iostream>
#include <cstdlib>
#include <boost/serialization/vector.hpp>


namespace mpi = boost::mpi;

enum message_tags { msg_data_packet, msg_broadcast_data, msg_finished };

void generate_data(mpi::communicator local, mpi::communicator world, bool imm, int c)
{
    using std::srand;
    using std::rand;

    // The rank of the collector within the world communicator
    int master_collector = 0;//local.size();

    srand(time(0) + world.rank());

    // Send out several blocks of random data to the collectors.
    int num_data_blocks = rand() % 3 + 1;
    for (int block = 0; block < num_data_blocks; ++block) {
        // Generate some random data
        int num_samples = 5;
        std::vector<int> data = {};
        for (int i = 0; i < num_samples; ++i) {
            data.push_back(1000000*(imm+1)  + 10000*(local.rank()+1) +  100*c + i);
        }

        // Send our data to the master collector process.
        std::cout << "Generator #" << local.rank() << " sends: ";
        for(auto el : data) std::cout << el << ", ";
        std::cout << std::endl;
        if(imm)
        {
            auto s = world.isend(master_collector, msg_data_packet, data);
        }
        else
        {
            world.send(master_collector, msg_data_packet, data);
        }
    }
}

void collect_data(mpi::communicator local, mpi::communicator world)
{
    // The rank of the collector within the world communicator
    int master_collector = 0;  

    if (world.rank() == 0) {
        while (true) {
            // Wait for a message
            mpi::status msg = world.probe();
            if (msg.tag() == msg_data_packet) {
                // Receive the packet of data
                std::vector<int> data;
                world.recv(msg.source(), msg.tag(), data);
                std::cout << "received data. format: #immediate 0 #rank+1 0 #packet 0 #sample" << std::endl;
                for(auto el : data) std::cout << el << ", ";
                std::cout << std::endl;
                // Tell each of the collectors that we'll be broadcasting some data
                //for (int dest = 1; dest < local.size(); ++dest)
                //  local.send(dest, msg_broadcast_data, msg.source());

                // Broadcast the actual data.
                //broadcast(local, data, 0);

            } else if (msg.tag() == msg_finished) {
                // Receive the message
                world.recv(msg.source(), msg.tag());

                // Tell each of the collectors that we're finished
                //for (int dest = 1; dest < local.size(); ++dest)
                //  local.send(dest, msg_finished);

                break;
            }
        }
    } /*  else {
          while (true) {
    // Wait for a message from the master collector
    mpi::status msg = local.probe();
    if (msg.tag() == msg_broadcast_data) {
// Receive the broadcast message
int originator;
local.recv(msg.source(), msg.tag(), originator);

// Receive the data broadcasted from the master collector
std::vector<int> data;
broadcast(local, data, 0);

std::cout << "Collector #" << local.rank()
<< " is processing data from generator #" << originator
<< "." << std::endl;
} else if (msg.tag() == msg_finished) {
    // Receive the message
    local.recv(msg.source(), msg.tag());

    break;
    }
    }
    } */
}

int main(int argc, char* argv[])
{
    mpi::environment env(argc, argv);
    mpi::communicator world;

    if (world.size() < 3) {
        if (world.rank() == 0) {
            std::cerr << "Error: this example requires at least 3 processes."
                << std::endl;
        }
        env.abort(-1);
    }

    bool is_generator = world.rank() > 0;
    mpi::communicator local = world.split(is_generator? 0 : 1);
    if (is_generator) 
    {
        for(int c = 0; c < 6; c++)
        {
            generate_data(local, world, false, c);
        }
    }
    else collect_data(local, world);

    (local.barrier)();
    if (local.rank() == 0)
    world.send(0, msg_finished);

    if(world.rank() == 0) std::cout << "immediate send" << std::endl;

    if (is_generator) 
    {
        for(int c = 0; c < 6; c++)
        {
            generate_data(local, world, true, c);
        }
    }
    else collect_data(local, world);

    (local.barrier)();
    if (local.rank() == 0)
    world.send(0, msg_finished);

    return 0;
}
