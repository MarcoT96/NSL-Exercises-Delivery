#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "random.h"  //Parallel Random Number Generator
		     //NY University


/*1d Quantum Particle*/
double psi(double, double, double);  // \psi_T^{\sigma, \mu} (x)
double pdf_T(double, double, double);  // | \psi_T^{\sigma, \mu} (x)|^2
double derivative_2(double, double, double);  //∂^2/∂x^2 \psi_T^{\sigma, \mu} (x)
double V(double);  // V(x) = x^4 - 5/2*x^2
double f(double, double, double);  //Integrand to evaluate
				   // (H_T*\psi_T^{\sigma, \mu} (x))/(\psi_T^{\sigma, \mu} (x))


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
  //Simulation
  int M=pow(10, 4);  //Number of MC throws
  int N_blk=pow(10, 2);  //Number of blocks
  int n=int(M/N_blk);  //Length of the single block
  double mu=0.5;       //Parameters of the Variational
  double sigma=1.5;    //Wave Function
  std::ofstream Energy, pdf;
  
  //Metropolis
  //double delta=1.9;  //HO
  double delta=2.7;  //moving parameter for M(RT)^2
  int accept=0;
  double x=0.0;  //sampled point
  double xstart=0.0;  //starting point in M(RT)^2
  double xnew=0.0;
  
  //Blocking
  double ene_blk=0.0;  //energy average in the single block
  double ene_ave=0.0;  //progressive energy average 
  double ene_ave2=0.0;  //progressive energy^2 average
  double ene_err=0.0;  //energy uncertainties

  
//Metropolis Algorithm
//Starting from xstart
//Blocking Method

  std::cout << "\nVariational Monte Carlo code for a single 1d Quantum Particle" << std::endl;
  std::cout << "The particle is confined by the potential V(x) = x^4 - 5/2 x^2" << std::endl;
  std::cout << "The code use natural units: hbar = 1, m = 1" << std::endl;
  std::cout << "Sampling the Variational pdf with the Metropolis algorithm" << std::endl << std::endl;

  //Energy.open("energy_HO.dat");
  //pdf.open("pdf_HO.dat");
  Energy.open("energy.dat");
  pdf.open("pdf.dat");
  
  x=xstart;
  std::cout << "Simulation\n" << "================================" << std::endl;
  for(int j=1; j<=N_blk; j++){
  
    std::cout << "Block " << j << std::endl;
    ene_blk=0.0;
    accept=0;
    for(int k=0; k<n; k++){
    
      //Sample with Metropolis
      xnew=rnd.Rannyu(x-delta, x+delta);  //Uniform Transition Probability
      if(rnd.Rannyu() < std::min(1., pdf_T(sigma, mu, xnew)/pdf_T(sigma, mu, x))){
        x=xnew;  
        accept++;
      }
      pdf << x << std::endl;
      ene_blk+=f(sigma, mu, x);
      
    }

    std::cout << "Acceptance ratio in the block = " << (accept*1.0)/(n*1.0) << std::endl << std::endl; 
    ene_ave+=ene_blk/double(n);
    ene_ave2+=pow(ene_blk/double(n), 2);
    if(j-1==0){
      ene_err=0;
      Energy << j << " " << ene_ave/double(j) << " " << ene_err << std::endl;
    }
    else{
      ene_err=sqrt((ene_ave2/double(j)-pow(ene_ave/double(j), 2))/double((j-1)));
      Energy << j << " " << ene_ave/double(j) << " " << ene_err << std::endl;
    }
    
  }

  Energy.close();
  pdf.close();
  rnd.SaveSeed();

return 0;
}


//Trial Wave Function
double psi(double sigma, double mu, double x){

  double a = pow(x-mu, 2)/(2*pow(sigma, 2));
  double b = pow(x+mu, 2)/(2*pow(sigma, 2));
  return exp(-a)+exp(-b);
  
}


//Trial Probability Density Function
//Square Modulus of the Trial WaveFunction
//Variational Monte Carlo (1D)
double pdf_T(double sigma, double mu, double x){

  return pow(psi(sigma, mu, x), 2);  //Real Wave Function

}


//External potential for the 1D
//quantum particle
double V(double x){

  return pow(x, 4) -  (5./2.)*(x*x);
  //return (1.0/2.0)*x*x;  //HO

}


double derivative_2(double sigma, double mu, double x){

  double e1 = exp((-1.0*(x-mu)*(x-mu))/(2.0*sigma*sigma));
  double e2 = exp((-1.0*(x+mu)*(x+mu))/(2.0*sigma*sigma));
  double p1 = ((x-mu)*(x-mu))/(sigma*sigma*sigma*sigma) - 1.0/(sigma*sigma);
  double p2 = ((x+mu)*(x+mu))/(sigma*sigma*sigma*sigma) - 1.0/(sigma*sigma);
  
  return e1*p1 + e2*p2;

}


//Integrand Function
double f(double sigma, double mu, double x){
  
  return ((-1.0)/(2*psi(sigma, mu, x)))*derivative_2(sigma, mu, x) + V(x);
  
}
