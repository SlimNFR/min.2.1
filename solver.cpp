//SOLVER MODULE - SOURCE FILE//
///////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

//Code defined libraries
#include "input.h"
#include "energy.h"
#include "gradient.h"
#include "utils.h"
#include "update.h"
#include "test.h"
#include "structure.h"
#include "fileout.h"
//namespace for the gradient terms
namespace solver{

//namespace variables
bool flag_reverse; //= true;
bool flag_print;

double alpha_spin;
double alpha_constr;
double tol_torque;
unsigned long int max_iter_spin;
unsigned long int max_iter_constr;
double vz_dif_tol;

unsigned long int t_iter;
unsigned long int l_iter;
unsigned long int s_iter;
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
	      std::vector<double> &lambda)
{
 
 unsigned long int iter_spin, iter_constr;
 total_iter=0;
 t_iter=0;

 //double alpha_spin=0.1;
 //double alpha_constr=0.1;
 //double tol_torque=1e-5;
 //unsigned long int max_iter_spin=10000;//1e6
 //unsigned long int max_iter_constr=10000;//1e5
 //double vz_dif_tol = 1e-3;
 double max_torque_x=0.0;
 double max_torque_y=0.0;
 double max_torque_z=0.0;
 double max_torque_mod=0.0;



// std::ofstream f_iter;
// f_iter.open("iter.txt", std::ofstream::app);
/* std::fill (sx.begin(),sx.end(),0.0);
 std::fill (sy.begin(),sy.end(),0.0);
 std::fill (sz.begin(),sz.end(),1.0);
*/
 double ax, ay, az;
 //double min.
 double vz_dif;

 double theta = 1.0*M_PI/180.0;
 double phi = 0.0*M_PI/180.0;



/* if(v0z_constraint<1.0&& v0z_constraint>0.0)
 {
 // sx[n_atoms-1]=ax;
 // sy[n_atoms-1]=ay;
 // sz[n_atoms-1]=az;

  sx[0]=ax;
  sy[0]=ay;
  sz[0]=-az;

 }*/
 bool f_reverse=false;
 double max_pair_angle; 
 std::vector<double>pair_angle;
 pair_angle.resize(n_atoms-1,0.0);
 for(int i=0; i<n_atoms-1;i++)
 {
  pair_angle[i] = (sx[i]*sx[i+1] + sy[i]*sy[i+1] + sz[i]*sz[i+1])/(1.0*1.0);   
  pair_angle[i] = acos(pair_angle[i]);
  pair_angle[i] = (pair_angle[i]*180.0)/M_PI;
  
    
 
 }

// utils::max_val_array(pair_angle,max_pair_angle);
 //if(max_pair_angle<2.0&&f_reverse==false)
 //{
   

 utils::compute_global_magnetisation(n_atoms,
                                    sx,sy,sz,
                                    global); 
 if(fabs(v0z_constraint - global[2])>=1e-5)
 {ax=sx[0]*cos(theta) + sz[0]*sin(theta);
 ay=sy[0];
 az=-sx[0]*sin(theta) + sz[0]*cos(theta);

 sx[0]=ax;
 sy[0]=ay;
 sz[0]=az;
 //f_reverse=true;
 }



 /*if(flag_reverse)
 {
 //sz[n_atoms-1]=az;
 sz[0]=az;
 }*/
 int print_step=5000;
 double print_v0z_start=-1.0;
 double print_v0z_end=1.0;
 int l;
 double tmp;
 tmp=v0z_constraint*100+0.5;
 l=(int)tmp;
 

  std::string file1 = "spin_and_lambda_vz" + std::to_string(l) + ".dat";
  std::string file2 = "grad_and_torque_vz" + std::to_string(l) + ".dat";
  std::string file3 = "energy_vz" + std::to_string(l) + ".dat";

  std::ofstream f1;
  std::ofstream f2;
  std::ofstream f3;

  std::ofstream log;
  log.open("solver_log.txt", std::ofstream::out);
  log<<"alpha spin: "<<alpha_spin<<"\n";
  log<<"alpha constr: "<<alpha_constr<<"\n";
  log<<"max iter spin: "<<max_iter_spin<<"\n";
  log<<"max iter constr: "<<max_iter_constr<<"\n";
  log<<"tol torque: "<<tol_torque<<"\n";
  log<<"tol vz dif"<<vz_dif_tol<<"\n";

/* if(v0z_constraint<=print_v0z_end && v0z_constraint>=print_v0z_start && (total_iter % print_step)==0)
 {
  f1.open(file1, std::ofstream::out);
  f2.open(file2, std::ofstream::out);
  f3.open(file3, std::ofstream::out);
 }
*/



 iter_constr=0;
 do
 {
 iter_spin=0;
  do
  {
   //print
 
   energy::compute_f();//compute the total energy
	
   gradient::total_f(n_atoms,
                     grad_total_x,
                     grad_total_y,
                     grad_total_z,
                     grad_real_x,
                     grad_real_y,
                     grad_real_z);

   gradient::remove_spin_projection(n_atoms,
                                    sx,sy,sz,
                                    grad_total_x,
                                    grad_total_y,
                                    grad_total_z,
                                    spin_dot_grad);
  /* if(v0z_constraint<=print_v0z_end && v0z_constraint >=print_v0z_start && (total_iter%print_step)==0)
   {
   fileout::spins_and_lambda(f1,n_atoms);
   fileout::gradient_and_torque(f2,n_atoms);
   fileout::energy(f3);
   }
  */
   magnitude_spin_grad=utils::compute_XYZROW_vec_magnitude(grad_total_x,
    		 		  	                   grad_total_y,
					 	 	   grad_total_z);
   if(magnitude_spin_grad<0.1){magnitude_spin_grad=1.0;}


  //std::cout<<"Computing torque.."<<"\n";
/*  gradient::torque_f(n_atoms,
                     sx,sy,sz,
                     grad_total_x,
                     grad_total_y,
                     grad_total_z,
                     torque_x,
                     torque_y,
                     torque_z,
                     torque_mod);
  utils::max_val_array(torque_x,max_torque_x);
  utils::max_val_array(torque_y,max_torque_y);
  utils::max_val_array(torque_z,max_torque_z);

  utils::max_val_array(torque_mod,max_torque_mod);
  */

   for(int i =0;i<n_atoms;i++)
   {
    double Heff_x = grad_total_x[i]/magnitude_spin_grad; //alpha_spin*grad_total_x[i]/magnitude_spin_grad;
    double Heff_y = grad_total_y[i]/magnitude_spin_grad;//alpha_spin*grad_total_y[i]/magnitude_spin_grad;
    double Heff_z = grad_total_z[i]/magnitude_spin_grad;//alpha_spin*grad_total_z[i]/magnitude_spin_grad;


    torque_x[i] = Heff_x;
    torque_y[i] = Heff_y;
    torque_z[i] = Heff_z;
    torque_mod[i] = sqrt(Heff_x * Heff_x +
			 Heff_y * Heff_y +
			 Heff_z * Heff_z);   

    /*torque_x[i] = sy[i]*Heff_z - sz[i]*Heff_y;
    torque_y[i] = sz[i]*Heff_x - sx[i]*Heff_z;
    torque_z[i] = sx[i]*Heff_y - sy[i]*Heff_x;

    torque_mod[i] = sqrt(torque_x[i]*torque_x[i] +
 		         torque_y[i]*torque_y[i] +
			 torque_z[i]*torque_z[i]);*/
 
   }
   utils::max_val_array(torque_x,max_torque_x);
   utils::max_val_array(torque_y,max_torque_y);
   utils::max_val_array(torque_z,max_torque_z); 
   utils::max_val_array(torque_mod,max_torque_mod);
  //utils::print_all(f1,f2,f3,iter_spin+total_iter);
 //std::cout<<"Update spins.."<<"\n"; 
   update::spins_dir_f(n_atoms,
                       sx,sy,sz,
                       spin_dir_x,
 		       spin_dir_y,
  		       spin_dir_z,
 		       grad_total_x,
 		       grad_total_y,
 		       grad_total_z,
		       magnitude_spin_grad,
		       alpha_spin);
  iter_spin++;
  total_iter=total_iter + 1;
  t_iter=t_iter + 1;
 }while(iter_spin<max_iter_spin &&
        max_torque_mod>tol_torque);

 utils::compute_global_magnetisation(n_atoms,
                                    sx,sy,sz,
                                    global);

 vz_dif = fabs(v0z_constraint - global[2]);
 if(vz_dif>=vz_dif_tol || 
    max_torque_mod>tol_torque)
 {

  gradient::constraints_wrt_lambda_f(n_atoms,
                                    v0z_constraint,
                                    sx,
                                    sy,
                                    sz,
                                    global,
                                    grad_lambda_param);
  magnitude_lambda_grad=utils::compute_ROW_vec_magnitude(grad_lambda_param);

  if(magnitude_lambda_grad<0.1){magnitude_lambda_grad=1.0;}

  //f_iter<<iter_spin<<" "<<magnitude_lambda_grad<<"\n";
  update::lambda_dir_f(lambda,
 		       lambda_dir,
		       grad_lambda_param,
		       magnitude_lambda_grad,
		       alpha_constr);

 
  if(fabs(lambda[0]) > 10*input::lambda_x ) {lambda[0]=input::lambda_x;}//input::lambda_x=-input::lambda_x;}//lambda[0]=input::lambda_x; }
  if(fabs(lambda[1]) > 10*input::lambda_y ) {lambda[1]=input::lambda_y;}//input::lambda_y=-input::lambda_y;}//lambda[1]=input::lambda_y; }
  if(fabs(lambda[2]) > 10*input::lambda_z ) {lambda[2]=input::lambda_z;}//input::lambda_z=-input::lambda_z;}//lambda[2]=input::lambda_z; }
  if(fabs(lambda[3]) > 10*input::l1 ) {lambda[3]=input::l1;}//input::l1=-input::l1;}//lambda[3]=input::l1; }


  iter_constr++;
  total_iter = total_iter + 1;
  t_iter = t_iter + 1;
 }
 else { total_iter = total_iter + iter_spin;
        t_iter = t_iter + iter_spin;
	std::cout<<"CORRECT";
	break;}
}while(iter_constr<max_iter_constr);//end do..while
std::cout<<"vz_dif: "<<vz_dif<<"|vz:"<<spin::global[2]<<"|v0z:"<<v0z_constraint<<"|constraint updated(x times):"<<iter_constr<<"|total_iter:"<<total_iter<<"|tq_x:"<<max_torque_x<<"|tq_y:"<<max_torque_y<<"|tq_z:"<<max_torque_z<<"|max_torque_mod:"<<max_torque_mod<<"\n";
 /*if(v0z_constraint<=print_v0z_end && v0z_constraint >=print_v0z_start && (total_iter%print_step)==0)
 {
 f1.close();
 f2.close();
 f3.close();
 }*/
 //f_iter.close();
log.close();
}//end sdescent function

void read()
{

//read Grad_descent parameters from file
 std::string s1;
 std::string input_file="SOLVER.txt";
 std::ifstream inFile(input_file);
 if(!inFile)
 {
  std::cout << std::endl << "Failed to open the file: "<< input_file <<std::endl;
  exit(1);
 }
    

inFile >> s1>>alpha_spin;
std::cout << s1 << "\t" << alpha_spin << "\t"<< std::endl;

inFile >> s1>>alpha_constr;
std::cout << s1 << "\t" << alpha_constr << "\t"<< std::endl;

inFile >> s1>>tol_torque;
std::cout << s1 << "\t" << tol_torque << "\t"<< std::endl;

inFile >> s1>>max_iter_spin;
std::cout << s1 << "\t" << max_iter_spin << "\t"<< std::endl;

inFile >> s1>>max_iter_constr;
std::cout << s1 << "\t" << max_iter_constr << "\t"<< std::endl;

inFile >> s1>>vz_dif_tol;
std::cout << s1 << "\t" << vz_dif_tol << "\t"<< std::endl;



}
}//end numerical_methods namespace
          
