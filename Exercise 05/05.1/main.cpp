#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"  //Parallel Random Number Generator
		     //NY University


double pdf_gs(double, double, double);
double pdf_2p(double, double, double);


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
  int M=pow(10, 6);  //Number of MC throws
  int N_blk=pow(10, 2);  //Number of blocks
  int n=int(M/N_blk);  //Length of the single block
  double delta_gs=1.2;  // \delta = 1.2*a0
  double delta_2p=3.0;  // \delta = 3.0*a0
  double delta_gs_gauss=0.7;  // \delta = 0.7*a0
  double delta_2p_gauss=1.9;  // \delta = 1.9*a0
  double sum_r=0.0, sum_r_gauss=0.0;
  int accept=0, accept_gauss=0;
  std::ofstream outfile;
  std::ofstream out3D;


//Point 1)
//Ground State
//Variables declaration
  double r_gs[3]={1., 1., 1.};
  double r_gs_gauss[3]={1., 1., 1.};
  double x=0, y=0, z=0, x_gauss=0, y_gauss=0, z_gauss=0;
  double r_gs_ave=0, r2_gs_ave=0, err_gs=0;
  double r_gs_ave_gauss=0, r2_gs_ave_gauss=0, err_gs_gauss=0; 
  
//Metropolis Algorithm
//Starting from {a0, a0, a0}
//Blocking Method
  std::cout << "\nHydrogen Atom Ground State and 2p Orbital Sampling" << std::endl;
  std::cout << "The programe uses the Metropolis algorithm with a Uniform and Normal Transition Probability" << std::endl;
  std::cout << "The step of the Transition Probability is chosen so as to have 50% of acceptance" << std::endl;
  std::cout << "The program uses Bohr radius units for distances" << std::endl;
  std::cout << "\nSimulation\nGround State\n" << "============================" << std::endl;
  outfile.open("ground_state.dat");
  out3D.open("gs_3D_distr.dat");
  //out3D << r_gs[0] << " " << r_gs[1] << " " << r_gs[2] << " " << r_gs_gauss[0] << " " << r_gs_gauss[1] << " " << r_gs_gauss[2] << std::endl; 
  for(int j=1; j<=N_blk; j++){
    sum_r=0.0;
    sum_r_gauss=0.0;
    accept=0;
    accept_gauss=0;
    for(int k=0; k<n; k++){
      //Uniform Transition Probability
      x=rnd.Rannyu(r_gs[0]-delta_gs, r_gs[0]+delta_gs);
      y=rnd.Rannyu(r_gs[1]-delta_gs, r_gs[1]+delta_gs);
      z=rnd.Rannyu(r_gs[2]-delta_gs, r_gs[2]+delta_gs); 
      if(rnd.Rannyu() < std::min(1., pdf_gs(x, y, z)/pdf_gs(r_gs[0], r_gs[1], r_gs[2]))){
        r_gs[0]=x;  
        r_gs[1]=y;
        r_gs[2]=z;
        accept++;
      }
      sum_r+=sqrt(pow(r_gs[0], 2)+pow(r_gs[1], 2)+pow(r_gs[2], 2));
     
      //Normal Transition Probability
      x_gauss=rnd.Gauss(r_gs_gauss[0], delta_gs_gauss);
      y_gauss=rnd.Gauss(r_gs_gauss[1], delta_gs_gauss);
      z_gauss=rnd.Gauss(r_gs_gauss[2], delta_gs_gauss);  
      if(rnd.Rannyu() < std::min(1., pdf_gs(x_gauss, y_gauss, z_gauss)/pdf_gs(r_gs_gauss[0], r_gs_gauss[1], r_gs_gauss[2]))){
        r_gs_gauss[0]=x_gauss;  
        r_gs_gauss[1]=y_gauss;
        r_gs_gauss[2]=z_gauss;
        accept_gauss++;
      }
      sum_r_gauss+=sqrt(pow(r_gs_gauss[0], 2)+pow(r_gs_gauss[1], 2)+pow(r_gs_gauss[2], 2));
      //Write on file all the sampled points
      out3D << r_gs[0] << " " << r_gs[1] << " " << r_gs[2] << " " << r_gs_gauss[0] << " " << r_gs_gauss[1] << " " << r_gs_gauss[2] << std::endl;
    }
    
    r_gs_ave+=sum_r/(n*1.0);
    r2_gs_ave+=pow(sum_r/(n*1.0), 2);
    r_gs_ave_gauss+=sum_r_gauss/(n*1.0);
    r2_gs_ave_gauss+=pow(sum_r_gauss/(n*1.0), 2);
    if(j-1==0){
      err_gs=0;
      err_gs_gauss=0;
      outfile << j << " " << r_gs_ave/(j*1.0) << " " << err_gs << " " << r_gs_ave_gauss/(j*1.0) << " " << err_gs_gauss <<std::endl;
    }
    else{
      err_gs=sqrt((r2_gs_ave/(j*1.0)-pow(r_gs_ave/(j*1.0), 2))/(1.0*(j-1)));
      err_gs_gauss=sqrt((r2_gs_ave_gauss/(j*1.0)-pow(r_gs_ave_gauss/(j*1.0), 2))/(1.0*(j-1)));
      outfile << j << " " << r_gs_ave/(j*1.0) << " " << err_gs << " " << r_gs_ave_gauss/(j*1.0) << " " << err_gs_gauss <<std::endl;
    }
    std::cout << "Completed Block " << j << std::endl; 
    std::cout << "Acceptance ratio in the block (Uniform case) = " << (accept*1.0)/(n*1.0) << std::endl;
    std::cout << "Acceptance ratio in the block (Normal case) = " << (accept_gauss*1.0)/(n*1.0) << std::endl << std::endl; 
  }
  
  outfile.close();
  out3D.close();
  

//Point 2)
//2p Orbital
//Variables declaration
  double r_2p[3]={1., 1., 1.};
  double r_2p_gauss[3]={1., 1., 1.};
  double r_2p_ave=0, r2_2p_ave=0, err_2p=0;
  double r_2p_ave_gauss=0, r2_2p_ave_gauss=0, err_2p_gauss=0;
  //sum_r=0.0;
  //sum_r_gauss=0.0;
  
//Metropolis Algorithm
//Starting from {a0, a0, a0}
//Blocking Method
  std::cout << "\nSimulation\n2p Orbital\n" << "============================" << std::endl;
  outfile.open("2p_orbital.dat");
  out3D.open("2p_3D_distr.dat");
  //out3D << r_2p[0] << " " << r_2p[1] << " " << r_2p[2] << " " << r_2p_gauss[0] << " " << r_2p_gauss[1] << " " << r_2p_gauss[2] << std::endl;
  for(int j=1; j<=N_blk; j++){
    sum_r=0.0;
    sum_r_gauss=0.0;
    accept=0;
    accept_gauss=0;
    for(int k=0; k<n; k++){
      //Uniform Transition Probability
      x=rnd.Rannyu(r_2p[0]-delta_2p, r_2p[0]+delta_2p);
      y=rnd.Rannyu(r_2p[1]-delta_2p, r_2p[1]+delta_2p);
      z=rnd.Rannyu(r_2p[2]-delta_2p, r_2p[2]+delta_2p);  
      if(rnd.Rannyu() < std::min(1., pdf_2p(x, y, z)/pdf_2p(r_2p[0], r_2p[1], r_2p[2]))){
        r_2p[0]=x;  
        r_2p[1]=y;
        r_2p[2]=z;
        accept++;
      }
      sum_r+=sqrt(pow(r_2p[0], 2)+pow(r_2p[1], 2)+pow(r_2p[2], 2));
      
      //Normal Transition Probability
      x_gauss=rnd.Gauss(r_2p_gauss[0], delta_2p_gauss);
      y_gauss=rnd.Gauss(r_2p_gauss[1], delta_2p_gauss);
      z_gauss=rnd.Gauss(r_2p_gauss[2], delta_2p_gauss);  
      if(rnd.Rannyu() < std::min(1., pdf_2p(x_gauss, y_gauss, z_gauss)/pdf_2p(r_2p_gauss[0], r_2p_gauss[1], r_2p_gauss[2]))){
        r_2p_gauss[0]=x_gauss;  
        r_2p_gauss[1]=y_gauss;
        r_2p_gauss[2]=z_gauss;
        accept_gauss++;
      }
      sum_r_gauss+=sqrt(pow(r_2p_gauss[0], 2)+pow(r_2p_gauss[1], 2)+pow(r_2p_gauss[2], 2));
      //Write on file all the sampled points
      out3D << r_2p[0] << " " << r_2p[1] << " " << r_2p[2] << " " << r_2p_gauss[0] << " " << r_2p_gauss[1] << " " << r_2p_gauss[2] << std::endl;
    }
    
    r_2p_ave+=sum_r/(n*1.0);
    r2_2p_ave+=pow(sum_r/(n*1.0), 2);
    r_2p_ave_gauss+=sum_r_gauss/(n*1.0);
    r2_2p_ave_gauss+=pow(sum_r_gauss/(n*1.0), 2);
    if(j-1==0){
      err_2p=0;
      err_2p_gauss=0;
      outfile << j << " " << r_2p_ave/(j*1.0) << " " << err_2p << " " << r_2p_ave_gauss/(j*1.0) << " " << err_2p_gauss <<std::endl;
    }
    else{
      err_2p=sqrt((r2_2p_ave/(j*1.0)-pow(r_2p_ave/(j*1.0), 2))/(1.0*(j-1)));
      err_2p_gauss=sqrt((r2_2p_ave_gauss/(j*1.0)-pow(r_2p_ave_gauss/(j*1.0), 2))/(1.0*(j-1)));
      outfile << j << " " << r_2p_ave/(j*1.0) << " " << err_2p << " " << r_2p_ave_gauss/(j*1.0) << " " << err_2p_gauss <<std::endl;
    }
    std::cout << "Completed Block " << j << std::endl; 
    std::cout << "Acceptance ratio in the block (Uniform case) = " << (accept*1.0)/(n*1.0) << std::endl;
    std::cout << "Acceptance ratio in the block (Normal case) = " << (accept_gauss*1.0)/(n*1.0) << std::endl << std::endl;
  }
  
  outfile.close();
  out3D.close();
  rnd.SaveSeed();

return 0;
}


//Ground State Probability Density Function
//Square Modulus of the WaveFunction
double pdf_gs(double x, double y, double z){

  return (1./M_PI)*exp(-2.*sqrt(pow(x, 2)+pow(y, 2)+pow(z, 2)));

}


//2p Orbital Probability Density Function
//Square Modulus of the WaveFunction
double pdf_2p(double x, double y, double z){

  return (1./(32*M_PI))*pow(z, 2)*exp(-sqrt(pow(x, 2)+pow(y, 2)+pow(z, 2)));

}
