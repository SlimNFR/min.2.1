//MAIN MODULE - SOURCE FILE///
//////////////////////////////


//----------------------------------------------------------------APPEALING NAME FOR ME CODE: BLABLA V.1.0-------------------------------------------------------//


//C++ libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <fenv.h>
#include <chrono>

//Code defined libraries
#include "input.h"
#include "screenout.h"
#include "init.h"
#include "energy.h"
#include "utils.h"
#include "structure.h"
#include "test.h"
#include "gradient.h"
#include "solver.h"
#include "sim.h"

int main()
{
//feenableexcept(FE_INVALID | FE_OVERFLOW);
 auto start = std::chrono::high_resolution_clock::now();
 input::read();
 init::alloc_memory();
 init::sys();
 solver::read();
 //test::gradient();
 //test::energy();
    test::spin();
    test::lambda();
//    test::gradient_no_parallel_comp();
    //test::steepest_descent();
 sim::energy_barrier(input::v0_constraint_z,
                     1.0,-1.0,0.01);


 auto finish = std::chrono::high_resolution_clock::now();
 std::chrono::duration<double> elapsed = finish - start;
 std::cout<<" Elapsed time:"<<elapsed.count() <<" s"<<"\n";

 return 0;
}
