/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


#ifndef __ISING__
#define __ISING__


//Random numbers
#include "random.h"
int seed[4];
Random rnd;


//Observables
//const int m_props=1000;
const int m_props = 10;
int n_props, iu, ic, im, ix;
double walker[m_props];


//Blocking Variables
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_u, stima_c, stima_m, stima_x;
double err_u, err_c, err_m, err_x;


//Configuration
const int m_spin=50;
double s[m_spin];


//Thermodynamical state variables
int nspin;
int restart;
double beta, temp, J, h;


//Simulation variables
int nstep, nblk, metro;


//Functions
void Input(std::string, std::string);
void Reset(int);
void Accumulate(void);
void Averages(int, std::string, std::string);
void Move(int);
void ConfFinal(std::string, std::string);
void Measure(std::string, std::string);
double Boltzmann(int, int);
int Pbc(int);
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
