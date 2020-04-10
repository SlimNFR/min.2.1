//INPUT MODULE - SOURCE FILE//
//////////////////////////////


//C++ libraries
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
//Code defined libraries
#include "input.h"
#include "energy.h"
#include "utils.h"
#include "structure.h"



//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES

namespace input{
//------------SYSTEM PARAMETERS-------------//
int nx;
int ny;
int nz;
int no_of_atoms;
//------------------------------------------//


//------------MATERIAL PARAMETERS-----------//
int n_materials; 
std::vector<int>natoms_material;
std::vector<double>exchange_list;
std::vector<double>spin_moment_list;
std::vector<double>anis_const_list;
std::vector<double> anis_dir_x;
std::vector<double> anis_dir_y;
std::vector<double> anis_dir_z;
double norm_factor;
//------------------------------------------//


//------------APPLIED FIELD-----------------//
double hx; 
double hy;
double hz;
double H;
//------------------------------------------//


//------------SPIN'S INITIAL ORIENTATION----//
std::vector<double>s0x;
std::vector<double>s0y;
std::vector<double>s0z;
//------------------------------------------//


//------------LAGRANGE PARAMETERS-----------//
double lambda_x;
double lambda_y;
double lambda_z;
double l1;
//------------------------------------------//


//------------ORIENTATION CONSTRAINT--------//
double v0_constraint_x;
double v0_constraint_y;
double v0_constraint_z;
//------------------------------------------//

void read()
{

 std::string s1;
 std::string input_file="INPUT.txt";
 std::ifstream inFile(input_file);


if(!inFile)
    {
        std::cout << std::endl << "Failed to open the file: "<< input_file <<std::endl;
        exit(1);
    }

//------------SYSTEM PARAMETERS-------------//

    inFile >> s1>> nx;
    std::cout << s1 << "\t" << nx << "\t"<< std::endl;
    inFile >> s1>> ny;
    std::cout << s1 << "\t" << ny << "\t"<< std::endl;
    inFile >> s1>> nz;
    std::cout << s1 << "\t" << nz << "\t"<< std::endl;

    no_of_atoms = nx*ny*nz;
    std::cout<<"no_of_atoms\t"<<no_of_atoms<<std::endl;
//------------------------------------------//


//------------MATERIAL PARAMETERS-----------// 
    inFile >>s1>>n_materials;
    std::cout << s1 << "\t" << n_materials << "\t"<< std::endl;

    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    int tmp;
    inFile>>tmp;
    natoms_material.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< natoms_material[i] << "\n";
    }




    inFile>>s1;
    for(int i=0;i<(2*n_materials-1);i++)
    {
    double tmp;
    inFile>>tmp;
    exchange_list.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< exchange_list[i] << "\n";
    }

    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    spin_moment_list.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< spin_moment_list[i] << "\n";
    }

    inFile>>s1;
   
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    anis_const_list.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< anis_const_list[i] << "\n";
    }

    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    anis_dir_x.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< anis_dir_x[i] << "\n";
    }

    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    anis_dir_y.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< anis_dir_y[i] << "\n";
    }
 
    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    anis_dir_z.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< anis_dir_z[i] << "\n";
    }
   
    inFile>>s1>>norm_factor;
    std::cout << s1 << "\t" << norm_factor << "\t"<< std::endl;
//------------------------------------------//


//------------APPLIED FIELD-----------------//
    inFile>>s1>>hx;
    std::cout << s1 << "\t" << hx << "\t"<< std::endl;


    inFile>>s1>>hy;
    std::cout << s1 << "\t" << hy << "\t"<< std::endl;

    inFile>>s1>>hz;
    std::cout << s1 << "\t" << hz << "\t"<< std::endl;

    inFile>>s1>>H;
    std::cout << s1 << "\t" << H << "\t"<< std::endl;
//------------------------------------------//


//------------SPIN'S INITIAL ORIENTATION----// 

    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    s0x.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< s0x[i] << "\n";
    }
   
    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    s0y.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< s0y[i] << "\n";
    }


    inFile>>s1;
    for(int i=0;i<n_materials;i++)
    {
    double tmp;
    inFile>>tmp;
    s0z.push_back(tmp);
    std::cout << s1 << "["<<i<<"]\t "<< s0z[i] << "\n";
    }
    


//------------------------------------------//


//------------LAGRANGE PARAMETERS-----------//
   inFile >> s1>> lambda_x;
   std::cout << s1 << "\t" << lambda_x << "\t"<< std::endl;
   inFile >> s1>> lambda_y;
   std::cout << s1 << "\t" << lambda_y << "\t"<< std::endl;
   inFile >> s1>> lambda_z;
   std::cout << s1 << "\t" << lambda_z << "\t"<< std::endl;
   inFile >> s1>> l1;
   std::cout << s1 << "\t" << l1 << "\t"<< std::endl;
//------------------------------------------//


//------------ORIENTATION CONSTRAINT--------//
   inFile >> s1>> v0_constraint_x; 
   std::cout << s1 << "\t" << v0_constraint_x << "\t"<< std::endl;
   inFile >> s1>> v0_constraint_y;
   std::cout << s1 << "\t" << v0_constraint_y << "\t"<< std::endl;
   inFile >> s1>> v0_constraint_z;
   std::cout << s1 << "\t" << v0_constraint_z << "\t"<< std::endl;
//------------------------------------------//



    
}



}
