//SIMULATION MODULE - SOURCE FILE//
///////////////////////////////////


//C++ libraries
#include<vector>
#include<cmath>
#include<string>

//Code defined libraries
#include "sim.h"
#include "input.h"
#include "energy.h"
#include "gradient.h"
#include "solver.h"
#include "update.h"
#include "structure.h"
#include "fileout.h"
namespace sim {


//namespace variables

//namespace functions


void energy_barrier(double &v0z,
	            double v0z_start,
	       	    double v0z_end,
	            double v0z_step
		   )
{ 
  std::string trend; //ascending or descending trend
  unsigned long int total_iter=0;//saves the total number of iterations
  int temp=fabs(v0z_start-v0z_end)/v0z_step+0.5;//+1); //the 0.5 constant makes sure there are no rounding errors 
  int npoints = temp+1;//total number of points in the energy barrier calculation;
 		       //the +1 term is necessary because the simulation would 
                       //otherwise run from v0z_start to v0z_end -1;
  int sign; //this variable will help increment/decrement the global constraint based on the "trend" variable

  if(v0z_start<=v0z_end)
  {
   trend="ascending";
   sign=1;
  }
  else if(v0z_start>v0z_end)
  {
   trend="descending";
   sign=-1;
  }
 
  //Create file pointers
  std::ofstream f1;
  std::ofstream f2;
  std::ofstream f3;

  //Create file names
  std::string f1_name = "en_barrier." + trend + ".txt";
  std::string f2_name = "spins_and_lambda." + trend + ".txt";  
  std::string f3_name = "gradient_and_torque." + trend + ".txt";

  //Open files for output 
  f1.open(f1_name, std::ofstream::out);
  f2.open(f2_name, std::ofstream::out);
  f3.open(f3_name, std::ofstream::out);




  v0z=v0z_start;//initialise the constraint 
 


  for(int i=1; i<=npoints; i++)
  {
  
    solver::sdescent(total_iter,
		     input::no_of_atoms,
		     v0z,
		     gradient::magnitude_spin_term,
		     gradient::magnitude_lambda_term,
		     gradient::total_x,
		     gradient::total_y,
		     gradient::total_z,
		     gradient::real_x,
		     gradient::real_y,
		     gradient::real_z,
		     gradient::lambda_param,
		     gradient::spin_dot_gradient,
		     gradient::torque_x,
		     gradient::torque_y,
		     gradient::torque_z,
		     gradient::torque_mod,
		     update::spin_dir_x,
		     update::spin_dir_y,
		     update::spin_dir_z,
		     update::lambda_dir,
		     spin::sx,
		     spin::sy,
		     spin::sz,
		     spin::global,
		     spin::lambda);

    energy::compute_f();//compute the total energy
    energy::total =(energy::total + (input::no_of_atoms-1))/input::no_of_atoms;
    //energy::total /=input::no_of_atoms;
    energy::real =(energy::real + (input::no_of_atoms-1))/input::no_of_atoms;
    //energy::real /= input::no_of_atoms;
    fileout::en_barrier(f1,total_iter,
			energy::total,
			energy::real,
			v0z,spin::global[2]);
    fileout::spins_and_lambda(f2,input::no_of_atoms);
    fileout::gradient_and_torque(f3,input::no_of_atoms);

    v0z=v0z + sign*v0z_step;//update v0z

  }

  //Close files
  f1.close();
  f2.close();
  f3.close();
}



}


