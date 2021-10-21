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
  // std::ifstream input("../../Parallel Random Number Generator/seed.in");
  std::ifstream input("../../Parallel Random Number Generator/new_seed.in");  //choose a different seed
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


//Point 1) and 2)
//Variables declaration
  int M=pow(10, 4);  //Number of MC throws
  int N_blk=pow(10, 2);  //Number of blocks
  int n=int(M/N_blk);  //Length of the single block
  double r=0;  //uniformly sampled variable
  double r_blk=0, r_ave=0, r2_ave=0, r_err=0;
  double sigma_blk=0, sigma_ave=0, sigma2_ave=0, sigma_err=0;
  std::ofstream outfile;

//Blocking Method to estimate <r> and σ^2  
  outfile.open("uniform.dat");
  //outfile.open("uniform_seed.dat");  //choose a different seed
  
  for(int j=1; j<=N_blk; j++){  //j is the number of blocks used
    r_blk=0, sigma_blk=0;
    for(int l=0; l<n; l++){
      r=rnd.Rannyu();  //Uniform Distribution in [0,1)
      r_blk+=r;
      sigma_blk+=pow((r-0.5),2);
    }
    
    //Effective observables in the single block
    r_ave+=r_blk/n;  //Sum{k=1,j} <r>_k  
    sigma_ave+=sigma_blk/n;  //Sum{k=1,j} σ^2_k
    r2_ave+=pow(r_blk/n, 2);  //Sum{k=1,j} <r>_k * <r>_k
    sigma2_ave+=pow(sigma_blk/n, 2);  //Sum{k=1,j} (σ^2)_k * (σ^2)_k

    //Statistical estimation
    if(j-1==0){
      r_err=0;
      sigma_err=0;
    }
    else{
      r_err=sqrt((r2_ave/j-pow(r_ave/j, 2))/(j-1));
      sigma_err=sqrt((sigma2_ave/j-pow(sigma_ave/j, 2))/(j-1));
    }
    outfile << j << " " << r_ave/j << " " << r_err << " " << sigma_ave/j << " " << sigma_err << std::endl;
    
  }
  
  outfile.close();


//Point 3)
//Variables declaration
  //Redefining M and n wrt the previous points
  M=pow(10,2);  //Number of sub-intervals in [0,1)
  n=pow(10, 4);  //Number of throws
  int exp=100;  //Number of experiments
  double th=n*1./M;  //Expected number of events
  		     //in each sub-interval (theory)
  double chi=0;  //Chi square
  double a=0, b=0;
  std::vector<int> n_i(M);  //Vector of RV fallen in each
  			    //sub-interval during a single
  			    //experiment
  
//Chisquare Test
  outfile.open("chisquare.dat"); 
  
  for(int j=1; j<=exp; j++){
    n_i.clear();  //cleaning up n_i at each experiment
    n_i.resize(M);
    chi=0;
    for(int k=0; k<n; k++){
      r=rnd.Rannyu();  //Uniform Distribution in [0,1)
      for(int l=0; l<M; l++){
        a=l*1./M;
        b=(l+1)*1./M;
        if(r>=a && r<b){
          n_i[l]+=1;
        }
      }
    }
    
    for(auto el : n_i){
      chi+=pow((el-th), 2)/th;
    }
    outfile << j << " " << chi << std::endl;
    
  }
  
  outfile.close();
  rnd.SaveSeed();

return 0;
}
