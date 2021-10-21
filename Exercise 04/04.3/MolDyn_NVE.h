/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


const int m_props=4;
int n_props;  //Number of observables
int iv, ik, it, ie;
double stima_pot, stima_kin, stima_etot, stima_temp;  //Istantaneous values
int nconf;


//Configuration variables
const int m_part=108;  //Number of particles in the system
double x[m_part], y[m_part], z[m_part];
double xold[m_part], yold[m_part], zold[m_part];  //this is for the restart option
double vx[m_part], vy[m_part], vz[m_part];


//Thermodynamical state variables
int npart;  //Number of the total particles
	    //m_part to use in a simulation
double energy, temp, vol, rho, box, rcut;
int restart;  //if = 1 Equilibrate
int blocking;  //if = 1 Blocking Method


//Simulation variables
int nstep, iprint, seed;  //nstep --> number of integration steps
double delta;  //Infinitesimal time step
	       //for the integration of the 
	       //equations of motion

	     
//Blocking variables
int n_blks, counter;
double sum_epot, sum_ekin, sum_etot, sum_temp;
double prog_epot=0, prog_ekin=0, prog_etot=0, prog_temp=0;
double prog2_epot=0, prog2_ekin=0, prog2_etot=0, prog2_temp=0;
double ave_epot, ave_ekin, ave_etot, ave_temp;
double ave2_epot, ave2_ekin, ave2_etot, ave2_temp;
double err_epot, err_ekin, err_etot, err_temp;


//Functions
void Input(std::string);
void Move(void);
void ConfFinal(std::string);
void ConfXYZ(int, std::string);
void Measure(std::string);
double Force(int, int);
double Pbc(double);
void Accumulate(void);
void Blocking(int, std::string);


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
