#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"  //Parallel Random Number Generator
		     //NY University


int main(){


/***************************************************************************************/
/* CREATE THE RANDOM GENERATOR */
  Random rnd;
  int seed[4];
  int p1, p2;
  std::ifstream Primes("../../Parallel Random Number Generator/Primes");
  if(Primes.is_open()){
    Primes >> p1 >> p2;
  }
  else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
  Primes.close();
  std::ifstream input("../../Parallel Random Number Generator/seed.in");
  std::string property;
  if(input.is_open()){
    while(!input.eof()){
      input >> property;
      if(property == "RANDOMSEED"){
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed,p1,p2);
      }
    }
    input.close();
  }
  else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
/***************************************************************************************/


//Variables declaration
  int M=pow(10, 4);  //Number of Random Walks (RW)
  int N_steps=pow(10, 2);  //Number of discrete time steps
  			   //for each RW
  double a=1.0;  //Lattice spacing
  double r[3]={0, 0, 0};  //Position of the walker 
  			  //at each discrete time step
  			  //this is the 3-vector (x, y, z)
  double r2=0;
  int dir=0;
  int vers=0;
  double err=0;
  std::vector<double> r2_ave(N_steps);  //vector of the average square position
  				        //of the walker at each discrete time step
  std::vector<double> r2_ave2(N_steps);  //to use in the blocking method			  
  std::ofstream outfile, RW;

//Point 1)
//3D RW in a Cubic Lattice
  outfile.open("RW_lattice.dat");
  RW.open("path_lattice.dat");
  
  //choose a particular walk
  int w=int(rnd.Rannyu(0, 100));
  
  for(int j=0; j<M; j++){
  
    //starting from the origin
    r[0]=0;  
    r[1]=0;
    r[2]=0;
    
    for(int k=0; k<N_steps-1; k++){
      dir=int(rnd.Rannyu(0, 3));  //choose a direction randomly
      if(rnd.Rannyu()<0.5)
        vers=+1;  //going forward
      else
        vers=-1;  //going backward
      
      //save the particular walk
      if(j==w){
        RW << r[0] << " " << r[1] << " " << r[2] << std::endl;
      }
      
      r[dir]+=a*vers;
      r2=pow(r[0], 2)+pow(r[1], 2)+pow(r[2], 2);  //|r_i|^2
      r2_ave[k+1]+=r2;
      r2_ave2[k+1]+=pow(r2, 2);
    }
  }
  
  //Blocking Method
  for(int l=0; l<N_steps; l++){
    r2_ave[l]=r2_ave[l]/M;
    r2_ave2[l]=r2_ave2[l]/M;
    err=(r2_ave2[l]-pow(r2_ave[l], 2))/float(M-1);
    if(l==0)
      outfile << sqrt(r2_ave[l]) << " " << 0. << std::endl;
    else
      outfile << sqrt(r2_ave[l]) << " " << 0.5*sqrt(err/r2_ave[l]) << std::endl;
  }
  
  outfile.close();
  RW.close();
 
 
//Point 2)
//Variables declaration
  //Solid angle in 3D
  float theta=0;  //ϑ angle in [0, π]
  float phi=0;  //φ angle in [0, 2π]
  r2_ave.clear();
  r2_ave2.clear();
  r2_ave.resize(N_steps);
  r2_ave2.resize(N_steps);

//3D RW in the Continuum  
  outfile.open("RW_continuum.dat");
  RW.open("path_continuum.dat");
  
  for(int j=0; j<M; j++){
    
    //starting from the origin
    r[0]=0;  
    r[1]=0;
    r[2]=0;
    
    for(int k=0; k<N_steps-1; k++){
      theta=rnd.Theta();  //choose a direction (solid angle) randomly
      phi=rnd.Rannyu(0, 2*M_PI);
      
      //save the particular walk
      if(j==w){
        RW << r[0] << " " << r[1] << " " << r[2] << std::endl;
      }
      
      r[0]+=a*sin(theta)*cos(phi);
      r[1]+=a*sin(theta)*sin(phi);
      r[2]+=a*cos(theta);
      r2=pow(r[0], 2)+pow(r[1], 2)+pow(r[2], 2);  //|r_i|^2
      r2_ave[k+1]+=r2;
      r2_ave2[k+1]+=pow(r2, 2);
    }
  }
  
  //Blocking Method
  for(int l=0; l<N_steps; l++){
    r2_ave[l]=r2_ave[l]/M;
    r2_ave2[l]=r2_ave2[l]/M;
    err=(r2_ave2[l]-pow(r2_ave[l], 2))/float(M-1);
    if(l==0)
      outfile << sqrt(r2_ave[l]) << " " << 0. << std::endl;
    else
      outfile << sqrt(r2_ave[l]) << " " << 0.5*sqrt(err/r2_ave[l]) << std::endl;
  }
  
  outfile.close();
  RW.close();
  rnd.SaveSeed();  
  
return 0;
}  
