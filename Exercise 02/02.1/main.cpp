#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"  //Parallel Random Number Generator
		     //NY University
		     

//Integrand Function
float f(float arg){
  return (M_PI/2.)*cos(arg*M_PI/2.);
}


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
  int M=pow(10, 4);  //Number of MC throws
  int N_blk=pow(10, 2);  //Number of blocks
  int n=int(M/N_blk);  //Length of the single block
  double a=0, b=1;  //Domain of integration [a,b]
  double x=0;
  double I_blk=0, I_ave=0, I_ave2=0, I_err=0;
  std::ofstream outfile;

//Blocking Method to estimate I
//Point 1)
//Using p(x) as Uniform Distribution
  outfile.open("I_uniform.dat");
  
  for(int j=1; j<=N_blk; j++){  //j is the number of blocks used
    I_blk=0;
    for(int l=0; l<n; l++){
      x=rnd.Rannyu();  //RV Uniformly distributed in [0,1)
      I_blk+=f(x);
    }
    
    //Effective observables in the single block
    I_ave+=((b-a)*I_blk/n);  //Sum{k=1,j} <f>_k
    I_ave2+=pow((b-a)*I_blk/n, 2);  //Sum{k=1,j} <f>_k * <f>_k
    
    //Statistical estimation
    if(j-1==0){
      I_err=0;
    }
    else{
      I_err=sqrt((I_ave2/j-pow(I_ave/j, 2))/(j-1));
    }
    
    outfile << j << " " << I_ave/j << " " << I_err << std::endl;  
  }
  
  outfile.close();

  
//Point 2)
//Using Importance Sampling
//Using p(x)=2*(1-x)
  I_blk=0, I_ave=0, I_ave2=0, I_err=0;
  outfile.open("I_importance.dat");
  
  for(int j=1; j<=N_blk; j++){  //j is the number of blocks used
    I_blk=0;
    for(int l=0; l<n; l++){
      x=1-pow(1-rnd.Rannyu(), 1/2.);  //RV distributed as p(x)
      		       		      //Using method of the inverse of the 
      		       		      //cumulative function 
      I_blk+=f(x)/(2*(1-x));  //Using the new integrand
      		  	    //according to Importance Sampling
    }
    
    //Effective observables in the single block
    I_ave+=((b-a)*I_blk/n);  //Sum{k=1,j} <f>_k
    I_ave2+=pow((b-a)*I_blk/n, 2);  //Sum{k=1,j} <f>_k * <f>_k
    
    //Statistical estimation
    if(j-1==0){
      I_err=0;
    }
    else{
      I_err=sqrt((I_ave2/j-pow(I_ave/j, 2))/(j-1));
    }
    
    outfile << j << " " << I_ave/j << " " << I_err << std::endl;
  }
  
  outfile.close();
  rnd.SaveSeed();  
  
return 0;
}
