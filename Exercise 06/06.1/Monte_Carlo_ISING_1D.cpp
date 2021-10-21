/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"


using namespace std;


int main(int argc, char** argv){

  if(argc != 3){
    cerr << endl <<"Error!!\nUsage: " << argv[0] << "  <path_to_input>"  << "  <path_to_save>" << endl << endl;
    return -1;
  }
  string path_input = argv[1];
  string path_save = argv[2];
  Input(path_input, path_save);  //Inizialization

  cout << endl << "Simulation" << endl << "=======================" << endl;
  cout << "Blocking Method to estimate observables averages and uncertainties" << endl;
  for(int iblk=1; iblk <= nblk; ++iblk){

    Reset(iblk);  //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){

      Move(metro);  //Move spins
      Measure(path_input, path_save);
      Accumulate();  //Update block averages

    }
    Averages(iblk, path_input, path_save);  //Print results for current block

  }
  ConfFinal(path_input, path_save);  //Write final configuration

return 0;
}


/**********************/
/*FUNCTIONS DEFINITION*/
/**********************/
//Prepare all stuff for the simulation
void Input(string folder, string folder1){
  ifstream ReadInput, ReadConf;

  cout << endl;
  cout << "Classic 1D Ising model" << endl;
  cout << "Monte Carlo simulation" << endl;
  cout << "Nearest neighbour interaction" << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl;
  cout << "The program uses mu_B = k_B = 1 units" << endl << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("../../Parallel Random Number Generator/Primes");
  if(Primes.is_open()){

    Primes >> p1 >> p2;
  }
  else cerr << "PROBLEM: Unable to open Primes" << endl; 
  Primes.close();
  ifstream input("../../Parallel Random Number Generator/seed.in");
  if(input.is_open()){

    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);

  }
  else cerr << "PROBLEM: Unable to open seed.in" << endl;
  input.close();
  
  //Parameter for the simulation
  ReadInput.open(folder + string("input.dat"));  //Read input
  if(ReadInput.is_open()){

    ReadInput >> temp;
    beta = 1.0/temp;
    ReadInput >> nspin;
    ReadInput >> J;
    ReadInput >> h;
    ReadInput >> metro;  //if=1 Metropolis else Gibbs
    ReadInput >> nblk;
    ReadInput >> nstep;
    ReadInput >> restart;  //if=1 restart from
                           //a previous spin configuration

  }
  else cerr << "PROBLEM: Unable to open the input file" << endl;
  			 
  cout << "Details of the system" << endl;
  cout << "=======================" << endl;
  cout << "Temperature = " << temp << endl;
  cout << "Number of spins = " << nspin << endl;
  cout << "Exchange interaction = " << J << endl;
  cout << "External field = " << h << endl << endl;
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  //Prepare arrays for measurements
  iu = 0;  //Energy
  ic = 1;  //Heat capacity
  im = 2;  //Magnetization
  ix = 3;  //Magnetic susceptibility
  n_props = 4;  //Number of observables

  //Initial configuration
  //Choose whether to equilibrate or not
  if(restart==0){

    cout << "No Equilibration" << endl;
    cout << "Random Initial configuration" << endl;
    for (int i=0; i<nspin; ++i){

      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;

    }

  }
  else{

    cout << "Equilibration" << endl;
    cout << "Read initial configuration from file config.final" << endl;
    ReadConf.open(folder + folder1 + string("config.final"));
    if(ReadConf.is_open()){

      for(int i=0; i<nspin; ++i)
        ReadConf >> s[i];

    }
    else cerr << "PROBLEM: Unable to open the configuration file" << endl;

  }
  
  //Evaluate observables of the initial configuration
  Measure(folder, folder1);
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  cout << "Heat Capacity = " << (beta*beta*((walker[ic]/(double)nspin)-pow(walker[iu]/(double)nspin, 2)))/double(nspin) << endl;
  cout << "Magnetization = " << walker[im]/(double)nspin << endl;
  cout << "Magnetic Susceptibility = " << (beta*walker[ix])/(double)nspin << endl;

}


//Move spin configuration
void Move(int metro){
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i){

    //Select randomly a particle
    //for C++ syntax, 0 <= o <= nspin-1
    o = (int)(rnd.Rannyu()*nspin);
    
    if(metro==1){  //Metropolis

      sm = s[o];
      energy_old = Boltzmann(sm, o);
      sm = -s[o];
      energy_new = Boltzmann(sm, o);
      
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu()){

        s[o] = sm;
        accepted++;

      }
      attempted++;

    }
    else{  //Gibbs sampling

      energy_up = Boltzmann(+1,o);
      energy_down = Boltzmann(-1,o); 

      p = exp(-beta*energy_up)/(exp(-beta*energy_up)+exp(-beta*energy_down));
      if(p >= rnd.Rannyu()){

        s[o] = 1;

      }
      else{

        s[o] = -1;

      }
      accepted++;
      attempted++;

    }

  }

}


//Compute the Boltzmann Weight
double Boltzmann(int sm, int ip){
 
  return -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;

}


//Properties measurement
void Measure(std::string folder, std::string folder1){
  double u = 0.0, m = 0.0;
  ofstream Ene, Mag;

  Ene.open(folder + folder1 + string("ene_istant.dat"), ios::app);
  Mag.open(folder + folder1 + string("mag_istant.dat"), ios::app);

  for (int i=0; i<nspin; ++i){

    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];

  }

  walker[iu] = u;
  walker[ic] = pow(u, 2);
  walker[im] = m;
  walker[ix] = pow(m, 2);

  Ene << u << std::endl;
  Mag << m << std::endl;

  Ene.close();
  Mag.close();

}


//Reset block averages
void Reset(int iblk){
   
  if(iblk == 1){

    for(int i=0; i<n_props; ++i){

      glob_av[i] = 0;
      glob_av2[i] = 0;

    }

  }

  for(int i=0; i<n_props; ++i){

    blk_av[i] = 0;

  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;

}


//Update block averages
void Accumulate(void){

  for(int i=0; i<n_props; ++i){

    blk_av[i] = blk_av[i] + walker[i];

  }

  blk_norm++;

}


//Print results for current block
void Averages(int iblk, string folder, string folder1){  
  ofstream Ene, Heat, Mag, Chi;
  const int wd=12;
    
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl;
    
  Ene.open(folder + folder1 + string("output.ene.dat"), ios::app);
  Heat.open(folder + folder1 + string("output.heat.dat"), ios::app);
  Mag.open(folder + folder1 + string("output.mag.dat"), ios::app);
  Chi.open(folder + folder1 + string("output.chi.dat"), ios::app);
  
  //Energy  
  stima_u = blk_av[iu]/blk_norm/(double)nspin;
  glob_av[iu] += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);

  //Heat capacity
  stima_c = blk_av[ic]/blk_norm; 
  stima_c = beta*(stima_c - blk_av[iu]/blk_norm * blk_av[iu]/blk_norm)/temp;
  stima_c = stima_c/(double)nspin;
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);

  //Magnetization
  stima_m = blk_av[im]/blk_norm/(double)nspin; 
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);

  //Susceptibility
  stima_x = blk_av[ix]/blk_norm;         
  stima_x = beta*stima_x/(double)nspin;
  glob_av[ix] += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x = Error(glob_av[ix],glob_av2[ix],iblk);
    
  Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
  Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
  Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;

  Ene.close();
  Heat.close();
  Mag.close();
  Chi.close();

}


void ConfFinal(string folder, string folder1){
  ofstream WriteConf;

  cout << endl << "Print final configuration to file config.final" << endl << endl;
  WriteConf.open(folder + folder1 + string("config.final"));
  for (int i=0; i<nspin; ++i){

    WriteConf << s[i] << endl;

  }
  WriteConf.close();
  rnd.SaveSeed();

}


//Algorithm for periodic boundary conditions
int Pbc(int i){
    
  if(i >= nspin) i = i - nspin;
  else if(i < 0) i = i + nspin;
  
  return i;

}


double Error(double sum, double sum2, int iblk){
    
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));

}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
