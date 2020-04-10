//RANDOM MODULE - SOURCE FILE//
///////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

//Code defined libraries
#include "rand.h"

//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES


namespace rrandom {


//namespace variables

bool init_config_flag=false;
bool init_lambda_flag=false;

//namespace functions

double config(int n_atoms,
	      double a,
	      double b,
	      std::vector<double> &sx,
	      std::vector<double> &sy,
	      std::vector<double> &sz)
{

    std::ofstream f1;
    f1.open("random_configs.txt",std::ios::out);
    int total = n_atoms;//sets the total number of spin orientations
    double seed;
    double theta,phi;
    double x,y,z;


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    seed=rd();//Save the seed
    f1<<"SEED: "<<seed<<"\n";//Output seed

    std::uniform_real_distribution<double> distribution(a,b);//Create real distribution data type
    std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with "seed"


    std::cout<<seed<<"\n";

   //Generate the random spin configurations
   for (int i=0; i<total; ++i)
    {
        theta = 2*M_PI*distribution(gen); 
        phi = acos(1-2*distribution(gen)); 
        x=sin(phi)*cos(theta);
        y=sin(phi)*sin(theta);
        z=cos(phi);
	sx[i]=x;
	sy[i]=y;
	sz[i]=z;
	f1<<theta<<" "<<phi<<" "<<x<<" "<<y<<" "<<z<<" "<<sqrt(x*x + y*y + z*z)<<"\n";
     }
}

double lambda(int n_total,
	      double a,
	      double b,
	      std::vector<double> &constraint)
{

 std::ofstream f1;
 f1.open("random_lambda.txt",std::ios::out);
 
 int total = n_total ;//sets the total sets of LAMBDA parameters to be generated
 double lx,ly,lz,l1;//lambda parameters
 double seed;

 std::random_device rd;  //Will be used to obtain a seed for the random number engine
 seed=rd();//Save the seed
 f1<<"SEED: "<<seed<<"\n";//Output seed

 std::uniform_real_distribution<double> distribution(a,b);//Create real distribution data type
 std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with "seed"


 std::cout<<seed<<"\n";


 for (int n=0; n<total; ++n)
 {
  constraint[0]=distribution(gen);
  constraint[1]=distribution(gen);
  constraint[2]=distribution(gen);
  constraint[3]=distribution(gen);

  f1<<constraint[0]<<" "<<constraint[1]<<" "<<constraint[2]<<" "<<constraint[3]<<" "<<"\n";

  }
  
  f1.close();
}


}//end of rand namespace
