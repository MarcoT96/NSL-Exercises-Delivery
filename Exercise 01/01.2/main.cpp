#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"  //Parallel Random Number Generator


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
  int M=pow(10, 4);  //Number of realizations
  int N[4]={1, 2, 10, 100};  //Number of rolling in each realization
  double u=0, e=0, l=0;  //RV for the 3 "dice"
  double S_u=0, S_e=0, S_l=0;  //Sum_{i=1, N} x_i
  std::ofstream outfile;
  
//Rolling the dice! 
  outfile.open("histo.dat");
  
  for(int j=0; j<4; j++){
    for(int k=0; k<M; k++){
      S_u=0, S_e=0, S_l=0;
      for(int m=1; m<=N[j]; m++){
        u=rnd.Rannyu();  //Uniform in [0,1)
        e=rnd.Exp(1.);  //Exponential with λ=1
        l=rnd.Lorentzian(1., 0);  //Lorentzian Γ=1, µ=0
        S_u+=u;
        S_e+=e;
        S_l+=l;
      }
      outfile << S_u/N[j] << " " << S_e/N[j] << " " << S_l/N[j] << "\n";
    }
    outfile << std::endl;
  } 
  
  outfile.close();
  rnd.SaveSeed();

return 0;
}
