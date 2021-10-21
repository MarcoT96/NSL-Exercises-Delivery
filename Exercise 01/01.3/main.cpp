#include <iostream>
#include <fstream>
#include <cmath>
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
  int N_thr=pow(10, 4);  //Number of MC throws
  int N=pow(10, 2);  //Number of blocks
  int n=int(N_thr/N);  //Length of the single block
  double L= 2.0;  //Length of the needle (cm)
  double d= 2.5;  //Distance between the straight lines (cm)
  	          //d > L (but not d >> L)
  double dist=0, x=0, y=0, theta=0;  
  int N_hit=0;  //Number of intersections between
  		//the needle and the line
  double pi_blk=0, pi2_blk=0, pi_ave=0, pi2_ave=0, pi_err=0;
  std::ofstream outfile;
  
//Blocking Method to estimate π
  outfile.open("Buffon.dat");
  
  for(int j=1; j<=N; j++){  //j is the number of blocks used
    N_hit=0;
    for(int l=0; l<n; l++){
      dist=rnd.Rannyu(0, d/2.);  //Distance between the midpoint of the needle
      				 //and the top of a straight line
      				 //Uniformly distributed in [0,d/2)
      do{
        x=rnd.Rannyu(0, 1);  //cartesian coordinate in the plane R^2
        y=rnd.Rannyu(0, 1);
      }while((pow(x, 2)+pow(y, 2)>1));  //Sampling an angle Uniformly in 
      					//[0, π/2) without using π itself
      					//Rejection technique
      theta=acos(x/sqrt(pow(x, 2)+pow(y, 2)));  //Angle between the needle and a straight line
      
      if(dist<=(L/2.)*sin(theta)){
        N_hit++;
      }
      
    }
    
    //Effective observables in the single block
    pi_blk+=(2*L/d)*(double(n)/double(N_hit));  //Compute <π> in each single block
    pi2_blk+=pow((2*L/d)*(double(n)/double(N_hit)), 2);
    
    //Statistical estimation
    pi_ave=pi_blk/j;  //1/j * Sum{k=1, j} <π>_k
    pi2_ave=pi2_blk/j;  //1/j * Sum{k=1, j} <π>_k * <π>_k
    if(j-1==0){
      pi_err=0;
    }
    else{
      pi_err=sqrt((pi2_ave-pow(pi_ave, 2))/(j-1));
    }
    outfile << j << " " << pi_ave << " " << pi_err << std::endl;
    
  }	      
  
  outfile.close();	      
  rnd.SaveSeed();

return 0;
}
