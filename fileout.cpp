//FILE OUTPUT MODULE - SOURCE FILE//
////////////////////////////////////



//C++ libraries
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

//Code defined libraries
#include "input.h"
#include "structure.h"
#include "energy.h"
#include "gradient.h"
#include "solver.h"
namespace fileout{


//namespace variables


//namespace functions


 void en_barrier(std::ofstream& f1,
		 unsigned long int& total_steps,
		 double& total_en,
		 double& real_en,
		 double& real_v0z,
		 double& total_v0z)
 {

 f1<<std::setprecision(10);
 f1<<total_steps<<" " //1
   <<total_en<<" "    //2
   <<real_en<<" "     //3
   <<real_v0z<<" "    //4
   <<total_v0z<<" "   //5
   <<input::v0_constraint_z<<" "//6
   <<spin::lambda[0]<<" "
   <<spin::lambda[1]<<" "
   <<spin::lambda[2]<<" "
   <<spin::lambda[3]<<" "
   <<"\n";
 }

 void spins_and_lambda(std::ofstream& f1,
		       int n_atoms)
 {
  f1<<std::setprecision(10);
 
  for(int i=0; i<n_atoms; i++)
  {
   f1<<i<<" "		 //1
     <<spin::sx[i]<<" "
     <<spin::sy[i]<<" "    
     <<spin::sz[i]<<" "
     <<sqrt(spin::sx[i]*spin::sx[i] + spin::sy[i]*spin::sy[i] + spin::sz[i]*spin::sz[i])<<" "//5
     <<spin::lambda[0]<<" "
     <<spin::lambda[1]<<" "
     <<spin::lambda[2]<<" "
     <<spin::lambda[3]<<" "//9
     <<solver::t_iter<<" "
     <<spin::global[0]<<" "//11
     <<spin::global[1]<<" "
     <<spin::global[2]<<" "//13
     <<solver::s_iter<<" "//14
     <<solver::l_iter<<" "//15
     <<input::v0_constraint_z<<"\n";//16
 
  }
  
  //f1<<"\n\n";
 }

 void gradient_and_torque(std::ofstream& f1,
		          int n_atoms)
 {
  f1<<std::setprecision(10);
  
  for(int i=0; i<n_atoms; i++)
  {
   f1<<i<<" "//1
   <<gradient::anisotropy_x[i]<<" "//2
   <<gradient::anisotropy_y[i]<<" "
   <<gradient::anisotropy_z[i]<<" "
   <<gradient::exchange_x[i]<<" "//5
   <<gradient::exchange_y[i]<<" "
   <<gradient::exchange_z[i]<<" "
   <<gradient::lagrange_x[i]<<" "//8
   <<gradient::lagrange_y[i]<<" "
   <<gradient::lagrange_z[i]<<" "
   <<gradient::lambda_param[0]<<" "//11
   <<gradient::lambda_param[1]<<" "
   <<gradient::lambda_param[2]<<" "
   <<gradient::lambda_param[3]<<" "
   <<gradient::real_x[i]<<" "//15
   <<gradient::real_y[i]<<" "
   <<gradient::real_z[i]<<" "
   <<gradient::total_x[i]<<" "//18
   <<gradient::total_y[i]<<" "
   <<gradient::total_z[i]<<" "
   <<gradient::torque_x[i]<<" "//21
   <<gradient::torque_y[i]<<" "
   <<gradient::torque_z[i]<<" "
   <<gradient::torque_mod[i]<<" "//24
   <<gradient::magnitude_spin_term<<" "//25
   <<gradient::magnitude_lambda_term<<" "//26
   <<solver::t_iter<<" "//27
   <<solver::s_iter<<" "//28
   <<solver::l_iter<<"\n";//29


  }
 
 //f1<<"\n\n";
 }
 
 void energy(std::ofstream& f1)
 {
  f1<<std::setprecision(10);
  
  f1<<energy::exchange<<" "//1
  <<energy::uniaxial_anisotropy<<" "
  <<energy::total<<" "
  <<energy::real<<" "
  <<energy::lagrange_term<<" "//5
  <<solver::t_iter<<"\n";//6

 }















/*
void print_all(std::ofstream& f1, std::ofstream& f2,std::ofstream& f3,int step)
{




for(int i = 0;i<input::no_of_atoms;i++)
{


   f1<<std::setprecision(10);
   f1<<i<< " " <<step<<" "<<spin::sx[i]<<" "<<spin::sy[i]<<" "<<spin::sz[i]<<" "<<sqrt(spin::sx[i]*spin::sx[i]+spin::sy[i]*spin::sy[i]+spin::sz[i]*spin::sz[i]) //1-6
        << " " << gradient::torque_x[i]<< " " << gradient::torque_y[i]<< " " << gradient::torque_z[i] <<" "<<gradient::torque_mod[i] //7-10
	<< " " <<gradient::total_x[i]<<" "<<gradient::total_y[i]<<" "<<gradient::total_z[i]<<" "<<gradient::lambda_param[0]<<" "<<gradient::lambda_param[1]<<" "
											        <<gradient::lambda_param[2]<<" "<<gradient::lambda_param[3]<<" "
												<<gradient::magnitude_spin_term<<" "<<gradient::magnitude_lambda_term// 11-19
	<< " " <<spin::lambda[0]<<" "<<spin::lambda[1]<<" "<<spin::lambda[2]<<" "<<spin::lambda[3]//20-23
	<< " " <<sqrt(gradient::exchange_x[i]*gradient::exchange_x[i] + gradient::exchange_y[i]*gradient::exchange_y[i] + gradient::exchange_z[i]*gradient::exchange_z[i])
	<< " " <<sqrt(gradient::anisotropy_x[i]*gradient::anisotropy_x[i] + gradient::anisotropy_y[i]*gradient::anisotropy_y[i] + gradient::anisotropy_z[i]*gradient::anisotropy_z[i])
	<< " " <<sqrt(gradient::lagrange_x[i]*gradient::lagrange_x[i] + gradient::lagrange_y[i]*gradient::lagrange_y[i] + gradient::lagrange_z[i]*gradient::lagrange_z[i])<<"\n";//24-26

*/
   /*
   f3 <<std::setprecision(10)  <<
     i <<" "<<spin::x[i]<<" "<<spin::y[i]<<" "<<spin::z[i]<<" " <<v0_constrain_z<<" "<<  // 1-5
	" "<<gradient::total[i]<<" "<<gradient::total[i+no_of_atoms]<<" "<<gradient::total[i+2*no_of_atoms]<< //6-8
	" "<< gradient::exchange_x[i]   <<" "<< gradient::exchange_y[i] <<" "<< gradient::exchange_z[i] << //9-11
	" "<< gradient::anisotropy_x[i] <<" "<< gradient::anisotropy_y[i] <<" "<< gradient::anisotropy_z[i] << //12-14
	" "<< gradient::zeeman_x[i]     <<" "<< gradient::zeeman_y[i] <<" "<< gradient::zeeman_z[i] << //15-17
	" "<< gradient::lagrangian_x[i] <<" "<< gradient::lagrangian_y[i] <<" "<< gradient::lagrangian_z[i] << //18-20
    " " << spin::lambda[0] << " "<< spin::lambda[1] << " "<< spin::lambda[2] << " "<< spin::lambda[3] << //21-24
    " " << gradient::torque_X[i] << " "<< gradient::torque_Y[i] << " "<< gradient::torque_Z[i] <<  //25-27
    " "<< step<<"\n";
*/




// f1<<spin::lambda[i]<<" "<<spin::lambda[i+1]<<" "<<spin::lambda[i+2]<<" "<<" "<<spin::lambda[i+3] << "  "<<gradient::total[i]<<" "<<gradient::total[i+no_of_atoms]<<" "<<gradient::total[i+2*no_of_atoms]<<" "<<gradient::total[i+3*no_of_atoms + 3]<<"\n";

//}
//f1<<"\n\n";


//f2<<std::setprecision(10);
//f2<<energy::real_energy<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";
//f2<<(energy::real_energy + J*(no_of_atoms-1))/(5)<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";
//f2<<(energy::real_energy + J*(no_of_atoms-1))/(no_of_atoms+0.0)<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";
//f2<<(energy::real_energy/(no_of_atoms))<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";



 //   <<gradient::total[*no_of_atoms]<<" "<<gradient::total[4*no_of_atoms + 1]<<gradient::total[4*no_of_atoms + 2];



// }


}

                                                                                                                                                                                                                                      
