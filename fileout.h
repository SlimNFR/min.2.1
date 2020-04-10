//FILE OUTPUT MODULE - HEADER FILE//
////////////////////////////////////



#ifndef FILEOUT_H
#define FILEOUT_H

//C++ libraries
#include <vector>
#include <fstream>

//Code defined libraries

namespace fileout{


//namespace variables


//namespace functions

 //output needed for energy surfaces plots
 void en_barrier(std::ofstream& f1,
                 unsigned long int& total_steps,
                 double& total_en,
                 double& real_en,
                 double& real_v0z,
                 double& total_v0z);
 //prints spin terms and lambda params
 void spins_and_lambda(std::ofstream& f1,
                       int n_atoms);

 void gradient_and_torque(std::ofstream& f1,
                          int n_atoms);

 void energy(std::ofstream& f1);

}

#endif //end FILEOUT_H

