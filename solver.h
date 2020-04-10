//SOLVER METHODS MODULE - HEADER FILE//
///////////////////////////////////////

#ifndef SOLVER_H
#define SOLVER_H

//C++ libraries
#include <vector>
#include <iostream>

//Code defined libraries

//namespace for the gradient terms
namespace solver{

//namespace variables
extern bool flag_reverse;
extern bool flag_print;

extern double alpha_spin;
extern double alpha_constr;
extern double tol_torque;
extern unsigned long int max_iter_spin;
extern unsigned long int max_iter_constr;
extern double vz_dif_tol;

unsigned long extern int t_iter;
unsigned long extern int l_iter;
unsigned long extern int s_iter;
//namespace functions

void sdescent(unsigned long int &total_iter,
	      int n_atoms,
              double v0z_constraint,
              double &magnitude_spin_grad,
              double &magnitude_lambda_grad,
              std::vector<double> &grad_total_x,
              std::vector<double> &grad_total_y,
              std::vector<double> &grad_total_z,
              std::vector<double> &grad_real_x,
              std::vector<double> &grad_real_y,
              std::vector<double> &grad_real_z,
	      std::vector<double> &grad_lambda_param,
	      std::vector<double> &spin_dot_grad,
              std::vector<double> &torque_x,
              std::vector<double> &torque_y,
              std::vector<double> &torque_z,
              std::vector<double> &torque_mod,
              std::vector<double> &spin_dir_x,
              std::vector<double> &spin_dir_y,
              std::vector<double> &spin_dir_z,
              std::vector<double> &lambda_dir,
              std::vector<double> &sx,
              std::vector<double> &sy,
              std::vector<double> &sz,
	      std::vector<double> &global,
              std::vector<double> &lambda);
void read();


}

#endif //SOLVER_METHODS_H
