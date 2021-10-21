#include <iostream>
#include <iomanip>
#include <fstream>

#define sigma_0 0.1
#define sigma_f 2.0 
#define mu_0 0.0
#define mu_f 2.0
#define N 100


int main(){

  /*
  
    Create a grid of points to 
    minimize with the Variational
    Monte Carlo GS Energy
  
  */

  double sigma_step = (sigma_f - sigma_0)/100.0;
  double mu_step = (mu_f - mu_0)/100.0;
  std::cout << "\nsigma_step = " << sigma_step << std::endl;
  std::cout << "mu_step = " << mu_step << std::cout << std::endl << std::endl;
  std::ofstream grid;
  
  grid.open("var_par.dat");
  //grid << "sigma" << std::setw(12) << "mu" << std::endl << std::endl;
  for(int j=0; j<N; j++){
    for(int k=0; k<N; k++){
      grid << std::fixed << std::setprecision(3) << sigma_0 + double(j)*sigma_step << " " <<  mu_0 + double(k)*mu_step << std::endl;
    }
  }
  grid.close();
  
return 0;
}
