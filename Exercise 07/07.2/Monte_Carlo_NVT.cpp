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
#include <iomanip>
#include <fstream>
#include "Monte_Carlo_NVT.h"


using namespace std;


int main(int argc, char** argv){

  if(argc != 2){
    cerr <<"Error!!\nUsage: " << argv[0] << "  <path_to_save>" << endl;
    return -1;
  }
  string path = argv[1];
   
  Input(path);  //Inizialization
  
  //Equilibration Phase
  cout << "\nEquilibration phase started" << endl;
  cout << "................" << endl;
  for(int eq_step=0; eq_step<equilibration_steps; eq_step++){
    Move();
  }
  cout << "Equilibration phase finished" << endl;

  //Actual Simulation
  cout << "\nActual simulation started" << endl;
  cout << "Simulation" << endl;
  cout << "Blocking Method to estimate observables averages and uncertainties" << endl;
  cout << "===============================================================================" << endl;
  for(int iblk=1; iblk <= N_blk; ++iblk){
  
    Reset(iblk);  //Reset block averages
    for(int istep=1; istep <= N_step; ++istep){
    
      Move();
      Measure(path);
      Accumulate();  //Update block averages
    
    }
    Averages(iblk, path);  //Print results for current block
    
  }
  ConfFinal(path);  //Write final configuration
 

return 0;
}


/**********************/
/*FUNCTIONS DEFINITION*/
/**********************/
//Prepare all stuff for the simulation
void Input(string path){
  ifstream ReadInput, ReadConf;

  cout << endl;
  cout << "Classic Lennard-Jones fluid" << endl;
  cout << "Monte Carlo simulation" << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T" << endl;
  cout << "The program uses Lennard-Jones units " << endl << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("../../Parallel Random Number Generator/Primes");
  if(Primes.is_open()){
    Primes >> p1 >> p2;
  }
  else cerr << "PROBLEM: Unable to open Primes" << std::endl;
  Primes.close();

  ifstream input("../../Parallel Random Number Generator/seed.in");
  string property;
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
  else cerr << "PROBLEM: Unable to open seed.in" << endl;
  
  //Parameters for the simulation
  ReadInput.open(path + string("input.dat"));
  if(ReadInput.is_open()){
    ReadInput >> temp;
    beta = 1.0/temp;
    ReadInput >> npart;
    ReadInput >> rho;
    vol = (double)npart/rho;
    box = pow(vol,1.0/3.0);
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> N_blk;
    ReadInput >> N_step;
    ReadInput >> restart;
    ReadInput >> equilibration_steps;
  }
  else cerr << "PROBLEM: Unable to open the input file" << endl;
  
  cout << "Details of the system" << endl;
  cout << "=======================" << endl;
  cout << "Temperature = " << temp << endl; 
  cout << "Number of particles = " << npart << endl;
  cout << "Density of particles = " << rho << endl;
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;
  cout << "Cutoff of the interatomic potential = " << rcut << endl;
    
  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial = " << ptail << endl; 
  
  cout << "\nDetails of the Simulation" << endl;
  cout << "=======================" << endl;
  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << N_blk << endl;
  cout << "Length of the blocks = " << N_step << endl;
  cout << "Number of total Monte Carlo steps = " << int(N_blk*N_step) << endl;
  cout << "Number of initial Monte Carlo equilibration steps = " << equilibration_steps << endl;
  ReadInput.close();

  //Prepare arrays for measurements
  iv = 0;  //Potential energy
  iw = 1;  //Virial
  n_props = 2;  //Number of observables

  igofr = 2;  //measurement of g(r)
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;  //for r \in [0, L/2]
  
  cout << "Number of bins fot the calculation of g(r) = " << nbins << endl;
  cout << "Linear dimension of each bin = " << bin_size << endl << endl;

  //Choose whether to equilibrate or not
  if(restart==0){
    cout << "No equilibration" << endl;
    cout << "Read initial configuration from file config.0" << endl;
    ReadConf.open(path + string("config.0"));
    if(ReadConf.is_open()){
      for(int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = Pbc( x[i] * box );
        y[i] = Pbc( y[i] * box );
        z[i] = Pbc( z[i] * box );
      }
    }
    else cerr << "PROBLEM: Unable to open the configuration file" << endl;
    ReadConf.close();
  }
  else if(restart==1){
    cout << "Equilibration" << endl;
    cout << "Read initial configuration from file old.0" << endl;
    ReadConf.open(path + string("old.0"));
    if(ReadConf.is_open()){
      for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = Pbc( x[i] * box );
        y[i] = Pbc( y[i] * box );
        z[i] = Pbc( z[i] * box );
      }
    }
    else cerr << "PROBLEM: Unable to open the configuration file" << endl;
    ReadConf.close();
  }
  
  //Evaluate potential energy and virial
  //of the initial configuration
  Measure(path);

  //Print initial values for
  //the potential energy and virial
  cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  cout << "Virial (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  cout << "Pressure (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}


//Move particles with Metropolis algorithm
void Move(void){
  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;

  for(int i=0; i<npart; ++i){
    //Select randomly a particle 
    //for C++ syntax, 0 <= o <= npart-1
    o = (int)(rnd.Rannyu()*npart);  
    				    
    //Old configuration
    xold = x[o];
    yold = y[o];
    zold = z[o];
    energy_old = Boltzmann(xold,yold,zold,o);

    //New configuration
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );
    energy_new = Boltzmann(xnew,ynew,znew,o);

    //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu()){
      //Update
      x[o] = xnew;
      y[o] = ynew;
      z[o] = znew;
      
      accepted++;
    }
    attempted++;
  }
}


//Compute the Boltzmann Weight
double Boltzmann(double xx, double yy, double zz, int ip){
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      //distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut){
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

return 4.0*ene;
}


//Properties measurement
void Measure(string path){
  int bin;
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Pres;

  //reset the hystogram of g(r)
  for(int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
      //distance i-j in pbc  
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);
      
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      
      //update of the histogram of g(r)
      bin = igofr + (int)(dr/bin_size);
      if(bin<(igofr + nbins))
        walker[bin] += 2.0;

      if(dr < rcut){
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
        //contribution to energy and virial
        v += vij;
        w += wij;
      }
    }          
  }
  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;
  
  Epot.open(path + string("epot.dat"),ios::app);
  Pres.open(path + string("pres.dat"),ios::app);
  Epot << walker[iv]/(double)npart + vtail << endl;
  //cout << "Virial (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  Pres << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl;
  
  Epot.close();
  Pres.close();
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
  blk_norm = blk_norm + 1.0;
}


//Print results for current block
void Averages(int iblk, string path){
  double r;
  ofstream Gofr, Gave, Epot, Pres;
  const int wd=12;
    
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
  Epot.open(path + string("output.epot.0"),ios::app);
  Pres.open(path + string("output.pres.0"),ios::app);
  Gofr.open(path + string("output.gofr.0"),ios::app);
  Gave.open(path + string("output.gave.0"),ios::app);
    
  stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
  stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
  glob_av[iw] += stima_pres;
  glob_av2[iw] += stima_pres*stima_pres;
  err_press=Error(glob_av[iw],glob_av2[iw],iblk);

  //Potential energy per particle
  Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
  
  //Pressure
  Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;

  //g(r)
  for(int ibin=igofr; ibin<igofr+nbins; ibin++){
    r = (ibin - igofr)*bin_size;
    stima_gofr = blk_av[ibin]/blk_norm;
    stima_gofr *= 3.0/((4.0)*M_PI * (pow(r + bin_size,3) - pow(r,3)) * rho * (double)npart);
    glob_av[ibin] += stima_gofr;
    glob_av2[ibin] += stima_gofr*stima_gofr;
    err_gofr = Error(glob_av[ibin],glob_av2[ibin],iblk);
    
    Gofr << setw(wd) << stima_gofr << setw(wd) << blk_norm << setw(wd) << r << endl;
    if(iblk == N_blk)
      Gave << scientific << r << "  " << glob_av[ibin]/(double)iblk << "  " << err_gofr << endl;
  }

  Epot.close();
  Pres.close();
  Gofr.close();
  Gave.close();
}


void ConfFinal(string path){
  ofstream WriteConf;
  ofstream WriteOld;

  //Write final configuration
  cout << "Print final configuration to file config.final" << endl;
  WriteConf.open(path + string("config.final"));
  for(int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  //Enable the possibility to (re)start reading
  //an old spatial configuration
  cout << "Saving final configuration also in old.0" << endl << endl;
  WriteOld.open(path + string("old.0"));
  for (int j=0; j<npart; ++j){
    WriteOld << x[j]/box << "   " <<  y[j]/box << "   " << z[j]/box << endl;
  }
  WriteOld.close();
}


/*
//Write configuration in .xyz format
void ConfXYZ(int nconf){ 
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}
*/


//Algorithm for periodic boundary conditions with side L=box
double Pbc(double r){
  return r - box * rint(r/box);
}


double Error(double sum, double sum2, int iblk){
 
  if( iblk == 1 ) return 0.0;
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
