#include <iostream>
#include <cmath>
#include <fstream>
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
  double sigma=0.0, mu=0.0;
  double mu_opt=0.0;       //Optimized Parameters of the
  double sigma_opt=0.0;    //Variational Wave Function
  double energy_min=INT_MAX;
  double energy=0.0;
  int flag=0;
  std::ofstream outfile, out_pdf, out_ene;
  std::ifstream ReadPar;  //Read the parameters grid
  
  //Metropolis
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
//Starting from x = xstart
//Blocking Method

  std::cout << "\nVariational Monte Carlo code for a single 1d Quantum Particle" << std::endl;
  std::cout << "The particle is confined by the potential V(x) = x^4 - 5/2 x^2" << std::endl;
  std::cout << "The code use natural units: hbar = 1, m = 1" << std::endl;
  std::cout << "Sampling the Variational pdf with the Metropolis algorithm" << std::endl << std::endl;


/*Optimization*/
  
  //Read file
  ReadPar.open("var_par.dat");
  out_ene.open("var_ene.dat");  //Python --> 3D plot of <H>_T (σ, µ) 
  
  if(!ReadPar){
    std::cerr << "Error: file could not be opened" << std::endl << std::endl;
    exit(1);
  }  
   
  x=xstart;
  std::cout << "Simulation\n" << "Search for Optimal Parameters\n" << "================================" << std::endl;
  std::cout << "\nRead a grid of points of the type (σ, µ) from file var_par.dat" << std::endl;
  out_ene << "sigma\t\t" << "mu\t\t" << "energy" << std::endl << std::endl;

  while(!ReadPar.eof()){  //keep reading until end-of-file
    
    flag++;
    if(flag%1000==0){
      std::cout << "Read & Optimize 1000 lines of (σ, µ)" << std::endl;
    }
    ReadPar >> std::setprecision(4) >> sigma >> mu;
    
    ene_ave=0.0;
    for(int j=0; j<M; j++){
    
      //Sample with Metropolis
      xnew=rnd.Rannyu(x-delta, x+delta);  //Uniform Transition Probability
      if(rnd.Rannyu() < std::min(1., pdf_T(sigma, mu, xnew)/pdf_T(sigma, mu, x))){
        x=xnew;  
      }
      ene_ave+=f(sigma, mu, x);
      
    }
    energy = ene_ave/M*1.0;
    //Saving each (σ, µ) and last average value
    //for the corresponding energy
    out_ene << sigma << "\t\t" << mu << "\t\t" << energy << std::endl;
    
    //Optimization
    //Check the minimum energy
    if(energy < energy_min){
      //Update
      energy_min=energy;
      sigma_opt=sigma;
      mu_opt=mu;
    }
    
  }
  
  ReadPar.close();
  out_ene.close();
  std::cout << "End-of-file reached..." << std::endl;
  std::cout << "Optimization Complete" << std::endl;
  std::cout << "Optimal Parameters found \t (σ, µ) = (" << sigma_opt << ", " << mu_opt << ")" << std::endl << std::endl;
  

/*Variational Ground State Energy*/
  
  outfile.open("ground_state.dat");
  out_pdf.open("pdf.dat");
  
  outfile << sigma_opt << " " << mu_opt << std::endl;
  x=xstart;
  std::cout << "Simulation\n" << "Variational Ground State Energy\n" << "================================" << std::endl;

  ene_ave=0.0;
  for(int j=1; j<=N_blk; j++){
    std::cout << "Block " << j << std::endl;
    ene_blk=0.0;
    accept=0;
    for(int k=0; k<n; k++){
      //Sample with Metropolis
      xnew=rnd.Rannyu(x-delta, x+delta);  //Uniform Transition Probability
      if(rnd.Rannyu() < std::min(1., pdf_T(sigma_opt, mu_opt, xnew)/pdf_T(sigma_opt, mu_opt, x))){
        x=xnew;  
        accept++;
      }
      ene_blk+=f(sigma_opt, mu_opt, x);
      out_pdf << x << std::endl;
    }
    std::cout << "Acceptance ratio in the block = " << (accept*1.0)/(n*1.0) << std::endl << std::endl; 
    ene_ave+=ene_blk/double(n);
    ene_ave2+=pow(ene_blk/double(n), 2);
    if(j-1==0){
      ene_err=0;
      outfile << j << " " << ene_ave/double(j) << " " << ene_err << std::endl;
    }
    else{
      ene_err=sqrt((ene_ave2/double(j)-pow(ene_ave/double(j), 2))/double((j-1)));
      outfile << j << " " << ene_ave/double(j) << " " << ene_err << std::endl;
    }  
  }
  outfile.close();
  out_pdf.close();
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
