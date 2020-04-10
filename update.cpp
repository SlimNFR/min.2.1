//UPDATE MODULE - SOURCE FILE///////
////////////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>

//Code defined libraries
#include "update.h"
#include "utils.h"

namespace update{

//namespace variables


std::vector<double> lambda_dir;
std::vector<double> spin_dir_x;
std::vector<double> spin_dir_y;
std::vector<double> spin_dir_z;

//namespace functions


void lambda_dir_f(std::vector<double> &lambda,
	          std::vector<double> &update_lambda,
		  std::vector<double> &grad_lambda_param,
		  double magnitude_lambda_grad,
		  double step)
{ 
 std::fill(update_lambda.begin(), update_lambda.end(),0);
 for(int i=0; i<lambda.size() ; i++)
 { 
  update_lambda[i]=grad_lambda_param[i]/magnitude_lambda_grad;	 
  lambda[i]-=step*update_lambda[i];
 } 

//  lambda[3]=0.0;
}

void spins_dir_f(int n_atoms,
	         std::vector<double> &sx,
                 std::vector<double> &sy,
	         std::vector<double> &sz,
                 std::vector<double> &update_sx,
	         std::vector<double> &update_sy,
	         std::vector<double> &update_sz,
		 std::vector<double> &grad_tot_x,
		 std::vector<double> &grad_tot_y,
		 std::vector<double> &grad_tot_z,
		 double &magnitude_spin_grad,
	         double step)
{
 std::fill(update_sx.begin(), update_sx.end(),0);
 std::fill(update_sy.begin(), update_sy.end(),0);
 std::fill(update_sz.begin(), update_sz.end(),0);
 for(int i=0; i<n_atoms; i++)
 {
  update_sx[i]=grad_tot_x[i]/magnitude_spin_grad;
  update_sy[i]=grad_tot_y[i]/magnitude_spin_grad;
  update_sz[i]=grad_tot_z[i]/magnitude_spin_grad;

  sx[i] -= step*update_sx[i];
  sy[i] -= step*update_sy[i];
  sz[i] -= step*update_sz[i];
  
  double magnitude = utils::compute_XYZ_vec_magnitude(sx[i],sy[i],sz[i]);
  sx[i] /= magnitude;
  sy[i] /= magnitude;
  sz[i] /= magnitude;
 }
}
}



//}
