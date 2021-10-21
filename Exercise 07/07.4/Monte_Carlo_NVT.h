/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __NVT__
#define __NVT__


#include <cmath>
/*
  Random numbers
  NY University
*/
#include "random.h"
int seed[4];
Random rnd;


//parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];


//Blocking Averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, stima_gofr, err_pot, err_press, err_gofr;


/*
  Thermodynamical state and
  Configuration variables
*/
const int m_part=108;
double x[m_part], y[m_part], z[m_part];
int npart;
double beta, temp, vol, rho, box, rcut;


//Simulation
int restart;  //whether to equilibrate or not
int equilibration_steps;  //Number of MC-steps of equilibration
			  //before measuring any observable
//const int M = 5*pow(10, 5);  //Monte Carlo steps
int N_step;  //Number of MC steps in each block
int N_blk;   //Number of blocks
int nconf;  //Number of saved configuration
double delta;


//pigreco
const double pi=M_PI;


//Functions
void Input(std::string);
void Reset(int);
void Accumulate(void);
void Averages(int, std::string);
void Move(void);
void ConfFinal(std::string);
void ConfXYZ(int, std::string);
void Measure(std::string);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);


#endif
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
