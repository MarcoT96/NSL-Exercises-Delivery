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
  double S_0=pow(10, 2);  //Asset price at t=0
  double T=1.;  //Delivery time
  double K=pow(10, 2);  //Strike price
  double r=0.1;  //Risk-free interest rate
  double sigma=0.25;  //Volatility
  int M=pow(10, 4);  //Number of asset prices
  int N_blk=pow(10, 2);  //Number of blocks
  int n=int(M/N_blk);  //Length of the single block
  std::ofstream outfile;
  
  
//Point 1)
//Sampling directly S(T)
//Variables declaration
  double W=0, sum_C=0, sum_P=0, S_T=0;
  double C_prog=0, C2_prog=0, C_ave=0, C2_ave=0, err_C=0;
  double P_prog=0, P2_prog=0, P_ave=0, P2_ave=0, err_P=0;
  
//Blocking Method to estimate C[S(0),0] and P[S(0),0]
  outfile.open("direct_sampling.dat");
  
  for(int j=1; j<=N_blk; j++){
    sum_C=0;
    sum_P=0;
    
    for(int l=0; l<n; l++){
      
      //Call-option price in the block
      W=rnd.Gauss(0, 1.0);
      S_T=S_0*exp((r-0.5*pow(sigma, 2))*T + sigma*W*sqrt(T));
      sum_C+=exp(-r*T)*std::max(0., S_T-K);
      
      //Put-option price in the block
      sum_P+=exp(-r*T)*std::max(0., -S_T+K);
    }
    
    C_prog+=(sum_C/n);
    C2_prog+=pow((sum_C/n), 2);
    C_ave=C_prog/j;
    C2_ave=C2_prog/j;
    P_prog+=(sum_P/n);
    P2_prog+=pow((sum_P/n), 2);
    P_ave=P_prog/j;
    P2_ave=P2_prog/j;
    if(j-1==0){
      err_C=0;
      err_P=0;
    }
    else{
      err_C=sqrt((C2_ave-pow(C_ave, 2))/(j-1));
      err_P=sqrt((P2_ave-pow(P_ave, 2))/(j-1));
    }
    
    outfile << j << " " << C_ave << " " << err_C << " " << P_ave << " " << err_P << std::endl;
  }
  
  outfile.close();
 
 
//Point 2)
//Sampling the discretized GBM(r, σ^2) of S(T)
//Variables declaration
  int time_steps=pow(10, 2);
  double dt=(T*1.)/(time_steps*1.);
  double S_t_C=0;  //Asset price at discrete time t
  		   //for the Call-option price
  double S_t_P=0;  //Asset price at discrete time t
  		   //for the Put-option price		     
  C_prog=0, C2_prog=0, C_ave=0, C2_ave=0, err_C=0;
  P_prog=0, P2_prog=0, P_ave=0, P2_ave=0, err_P=0;
  
  //choose a particular GBM
  int p=int(rnd.Rannyu(0, 100));
  std::ofstream GBM;

//Blocking Method to estimate C[S(0),0] and P[S(0),0]
  outfile.open("discrete_sampling.dat");
  GBM.open("GBM.dat");
  
  for(int j=1; j<=N_blk; j++){
    sum_C=0;
    sum_P=0;
    for(int l=0; l<n; l++){
      S_t_C=S_0;
      S_t_P=S_0;
      
      //Call-option price in the block
      //Discretized path
      for(int t=0; t<time_steps; ++t){
        W=rnd.Gauss(0, 1.0);
        S_t_C=S_t_C*exp((r-0.5*pow(sigma, 2))*dt + sigma*W*sqrt(dt));
        if(j==p){
          GBM << S_t_C << std::endl;
        }
      }
      
      //Put-option price in the block
      //Discretized path
      for(int t=0; t<time_steps; ++t){
        W=rnd.Gauss(0, 1);
        S_t_P=S_t_P*exp((r-0.5*pow(sigma, 2))*dt + sigma*W*sqrt(dt));
      }
      
      sum_C+=exp(-r*T)*std::max(0., S_t_C-K);
      sum_P+=exp(-r*T)*std::max(0., K-S_t_P);
    }
    
    C_prog+=(sum_C/n);
    C2_prog+=pow((sum_C/n), 2);
    C_ave=C_prog/j;
    C2_ave=C2_prog/j;
    P_prog+=(sum_P/n);
    P2_prog+=pow((sum_P/n), 2);
    P_ave=P_prog/j;
    P2_ave=P2_prog/j;
    if(j-1==0){
      err_C=0;
      err_P=0;
    }
    else{
      err_C=sqrt((C2_ave-pow(C_ave, 2))/(j-1));
      err_P=sqrt((P2_ave-pow(P_ave, 2))/(j-1));
    }
    outfile << j << " " << C_ave << " " << err_C << " " << P_ave << " " << err_P << std::endl;
  }
  
  outfile.close();
  GBM.close();
  rnd.SaveSeed();

return 0;
}
