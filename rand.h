//RANDOM MODULE - HEADER FILE//
///////////////////////////////

#ifndef RAND_H
#define RAND_H


//C++ libraries
#include <vector>
#include <fstream>


//Code defined libraries
namespace rrandom {


//namespace variables

extern bool init_config_flag;
extern bool init_lambda_flag;

//namespace functions

//---function to generate random initial spin configs
     //---a and b-> distribution range (Should be 0:1)
double config(int n_atoms,
	      double a,
	      double b,
              std::vector<double> &sx,
              std::vector<double> &sy,
              std::vector<double> &sz);

//---function to generate random initial lambda params.
     //---n_total-> number of trials
     //---a and b-> distribution range
double lambda(int n_total,
	      double a,
	      double b,
              std::vector<double> &constraint);

}



#endif // RAND_H
