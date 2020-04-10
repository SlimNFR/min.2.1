 //UTILS MODULE - SOURCE FILE//
//////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip> //setprecision
#include <iterator>//std::distance

//Code defined libraries
#include "utils.h"
#include "input.h"
#include "structure.h"
#include "gradient.h"



//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES


namespace utils {


//namespace variables

std::vector<std::vector<double> >Jmatrix;
std::vector<double>global_magnetisation = {0,0,0};
std::vector<double> sums = {0,0,0,0};
std::vector<int>i_list;
std::vector<int>start;
std::vector<int>end;




std::vector<double>temp_x;
std::vector<double>temp_y;
std::vector<double>temp_z;
std::vector<double>temp_lambda;
std::vector<double>temp_gradient;
double temp_energy;


//namespace functions
//
//
 void generate_exchange_mat(int n_materials,
			    std::vector<double>&exch_list,
			    std::vector<std::vector<double> >&Jmat)
 {
  Jmat.resize(n_materials, std::vector<double>(n_materials));

  for(int i=0;i<n_materials;i++)
  {
   for(int j=0;j<n_materials;j++)
   {
    if(abs(j-i)>1.5)//needs to be only greater than 1, but for safety..
    { 
     Jmat[i][j]=0;
     continue;
    }
    else
    {
     Jmat[i][j] = exch_list[i+j];
    }    
   }
  }


 }

 void generate_interac_list(int n_atoms,
			    std::vector<double> &x,
			    std::vector<double> &y,
			    std::vector<double> &z,
			    std::vector<int> &int_list,
			    std::vector<int> &start,
			    std::vector<int> &end)
 {
  start.resize(n_atoms);
  end.resize(n_atoms);
  double r0=1.0;//range of interactions(currently equal with the unit cell step)
  const double rtoll=r0*1.0e-5;

  for(int i=0; i<n_atoms; i++)
  {
   int count_interactions=0;//count the total interactions an atom has
   for(int j=0; j<n_atoms; j++)
    {
      if(i==j)continue;//exclude self interactions
      else
      {
      //compute the distance between atoms
      double rij=sqrt( pow((x[i]-x[j]),2.0) +
		       pow((y[i]-y[j]),2.0) +
		       pow((z[i]-z[j]),2.0) );
      if(rij-r0<rtoll)//if the distance agrees with the interaction range..
      {
       count_interactions++;
       int_list.push_back(j);//update interaction list
       int index_j=std::distance(int_list.begin(),int_list.end());
       end[i]=index_j;//update end list
      }
      } 
    }
    start[i]=end[i]-count_interactions;
 }

 }


 void compute_global_magnetisation(int n_atoms,
				   std::vector<double> &sx,
				   std::vector<double> &sy,
				   std::vector<double> &sz,
				   std::vector<double> &global)
 {
  double modulus = 0, term1 = 0, term2 = 0, term3 = 0;
  
  for(int count_atom=0;count_atom<n_atoms;count_atom++)
  {
   term1 = term1 + sx[count_atom];
   term2 = term2 + sy[count_atom];
   term3 = term3 + sz[count_atom];

  }

   modulus = sqrt( pow(term1,2.0) + pow(term2,2.0) + pow(term3,2.0)  );
   global[0] = term1/n_atoms;
   global[1] = term2/n_atoms;
   global[2] = term3/n_atoms;
  


 }


 double compute_XYZROW_vec_magnitude(std::vector<double> &x,
				     std::vector<double> &y,
				     std::vector<double> &z)
 {

  double sum=0.0;	
  for(int i=0;i<x.size();i++)
  {
   sum += x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
  }
  return sqrt(sum);
  
 }

 double compute_XYZ_vec_magnitude(double &x,
				  double &y,
				  double &z)
 {

  return sqrt(x*x + y*y + z*z);
 }

 double compute_ROW_vec_magnitude(std::vector<double> &vec)
 {
  double sum=0.0;
  for(int i=0; i<vec.size();i++)
  {
   sum +=vec[i]*vec[i];
  }
   return sqrt(sum);
 }




 double max_val_array(std::vector<double> &vec, double& max_val)
 {

  max_val = fabs(vec[0]);

  for (int i = 0; i<vec.size(); i++)
  {
   if(max_val<fabs(vec[i]))max_val=fabs(vec[i]);
  }

  return max_val;
 }
}
