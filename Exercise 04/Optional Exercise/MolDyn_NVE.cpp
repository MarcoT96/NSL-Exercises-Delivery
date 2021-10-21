/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"


using namespace std;


int main(int argc, char** argv){

  if(argc != 2){
    cerr <<"Error!!\nUsage: " << argv[0] << "  <path_to_save>" << endl;
    return -1;
  }
  string path = argv[1];
  Input(path);  //Inizialization
  
  if(blocking==0){
  
    cout << endl << "Simulation" << endl << "=======================" << endl;
    cout << "Computing only istantaneous observables during the equilibration phase" << endl;
    for(int istep=1; istep <= nstep; ++istep){
    
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      Move();  //Move particles with Verlet algorithm
      Measure(path);  //Istantaneous properties measurement
        
    }
    
  }
  else{
  
    int n = nstep/n_blks;  //Length of the single block
    cout << endl << "Simulation" << endl << "=======================" << endl;
    cout << "Blocking Method to estimate observables averages and uncertainties" << endl;
    for(int j=1; j<=n_blks; ++j){

      sum_epot = 0.0;
      sum_ekin = 0.0;
      sum_etot = 0.0;
      sum_temp = 0.0;
      sum_pres = 0.0;
      
      for(int k=0; k<nbins; ++k)
        sum_gofr[k]=0;
      
      counter = 0;
      for(int istep=1; istep <= n; ++istep){

        Move();  //Move particles with Verlet algorithm
        Measure(path);
        counter++;
        Accumulate();  //Accumulate block observables
      
      }
      Blocking(j, path);  //Compute observables with the 
       		          //Blocking Method
      cout << "Completed block " << j << endl;

    } 
    
  }
   
  ConfFinal(path);  //Write final configuration to restart
  Delete_g();

return 0;
}


/**********************/
/*FUNCTIONS DEFINITION*/
/**********************/
//Prepare all stuff for the simulation
void Input(string path){
  ifstream ReadInput, ReadConf;
  
  //Units of measure LJ --> SI
  //Use only to convert the parameters
  double SIGMA=0.34*pow(10, -9);  //[𝜎] = m
  double EPS=120;  //[𝜖/𝑘b] = K
  double Kb=1.38064852*pow(10, -23);  //[Kb] = m^2 kg s^(-2) K^(-1)
  double M=6.6335209*pow(10, -26); //[m] = kg --> m=39.948 amu
  double unit_of_energy=EPS*Kb;
  double unit_of_time=SIGMA*sqrt(M/unit_of_energy);

  cout << endl;
  cout << "Classic Lennard-Jones system" << endl;
  cout << "Molecular dynamics simulation in NVE ensemble" << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
  cout << "The program uses Lennard-Jones units" << endl << endl;

  seed = 1;  //Set seed for random numbers
  srand(seed);  //Initialize random number generator
  
  //Parameter for the simulation  
  ReadInput.open(path + string("input.dat"));  //Read input
  if(ReadInput.is_open()){

    ReadInput >> temp;
    ReadInput >> npart;
    ReadInput >> rho;
    vol = (double)npart/rho;
    box = pow(vol,1.0/3.0);
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> n_blks;
    ReadInput >> iprint;
    ReadInput >> restart;
    ReadInput >> blocking;

  }
  else cerr << "PROBLEM: Unable to open the input file" << endl;
  
  cout << "Details of the system" << endl;
  cout << "=======================" << endl;
  cout << "Target Temperature = " << temp << "\t  ( " << EPS*temp << "K )" << endl;
  cout << "Number of particles = " << npart << endl;
  cout << "Density of particles = " << rho << "\t  ( " << rho*pow(SIGMA, -3) << "m^(-3) )" << endl;
  cout << "Volume of the simulation box = " << vol << "\t ( " << vol*pow(SIGMA, 3) << "m^3 )" << endl;
  cout << "Edge of the simulation box = " << box << "\t ( " << box*SIGMA << "m )" << endl;
  cout << "Cutoff of the interatomic potential = " << rcut << "\t ( " << rcut*SIGMA << "m )" << endl;
  
  //Tail corrections for potential energy and pressure
  vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << "\t ( " << vtail*unit_of_energy << "J )" << endl;
  cout << "Tail correction for the virial = " << ptail << endl;
  
  cout << "\nDetails of the Simulation" << endl;
  cout << "=======================" << endl;
  cout << "The program integrates Newton equations with the Verlet method" << endl;
  cout << "Time step = " << delta << "\t ( " << delta*unit_of_time << "s )" << endl;
  cout << "Number of integration steps = " << nstep << endl;
  cout << "Number of blocks = " << n_blks << endl;
  ReadInput.close();

  //Prepare array for measurements
  iv = 0;  //Potential energy
  ik = 1;  //Kinetic energy
  ie = 2;  //Total energy
  it = 3;  //Temperature
  iw = 4;  //Virial
  n_props = 5;  //Number of observables
  
  //measurement of g(r)
  nbins = 100;
  bin_size = (box/2.0)/(double)nbins;
  gofr = new double[nbins];
  sum_gofr = new double[nbins];
  ave_gofr = new double[nbins];
  ave2_gofr = new double[nbins];
   
  cout << "Number of bins for the calculation of g(r) = " << nbins << endl;
  cout << "Linear dimension of each bin = " << bin_size << "\t ( " << bin_size*SIGMA << "m )" << endl << endl;

  //Choose whether to equilibrate or not
  if(restart==0){

    cout << "No equilibration" << endl;
    cout << "Read initial configuration from file config.0" << endl;
    ReadConf.open(path + string("config.0"));
    if(ReadConf.is_open()){

      for (int i=0; i<npart; ++i){

        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;

      }

    }
    else cerr << "PROBLEM: Unable to open the configuration file" << endl;
    ReadConf.close();

    //Prepare initial velocities randomly
    cout << "Prepare random velocities (avoiding drift motion)" << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){

      vx[i] = rand()/double(RAND_MAX) - 0.5;
      vy[i] = rand()/double(RAND_MAX) - 0.5;
      vz[i] = rand()/double(RAND_MAX) - 0.5;
   
      //Calculate the velocity center of mass
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];

    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){

      //P_tot=0
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2]; 

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];

    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   //fs = velocity scale factor 
    for (int i=0; i<npart; ++i){

      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);

    }

  }
  else if(restart==1){

    cout << "Equilibration" << endl;
    cout << "Read initial configuration from file old.0" << endl;
    cout << "Read not only r(t) but also r(t-dt) from a previous simulation" << endl; 
    ReadConf.open(path + string("old.0"));
    if(ReadConf.is_open()){

      //Point 1)
      for (int i=0; i<npart; ++i){

        ReadConf >> x[i] >> y[i] >> z[i] >> xold[i] >> yold[i] >> zold[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
      
        xold[i] = xold[i] * box;
        yold[i] = yold[i] * box;
        zold[i] = zold[i] * box;

      }

    }
    else cerr << "PROBLEM: Unable to open the configuration file" << endl;
    ReadConf.close();
    
    //Point 2)
    Move();  //After this we have x[i]=r_x(t+dt)
    	     //and xold[i]=r(t)
    	     //and the same for the other two components
    double sumv2 = 0.0, fs;
    for(int j=0; j<npart; ++j){  //v[j]=v_j(t+dt/2)

      vx[j] = (x[j]-xold[j])/(delta*1.0);
      vy[j] = (x[j]-xold[j])/(delta*1.0);
      vz[j] = (x[j]-xold[j])/(delta*1.0);
      sumv2 += pow(vx[j], 2) + pow(vy[j], 2) + pow(vz[j], 2); 

    }
    sumv2 /= (double)npart;
    
    //Point 3)
    fs = sqrt(3 * temp / sumv2);   //fs = velocity scale factor 
    for (int i=0; i<npart; ++i){

      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
      
      //Point 4) & 5)
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);   

    }

  }

}


//Move particles with Verlet algorithm
void Move(void){ 
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ 

    //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);

  }

  //Verlet integration
  for(int i=0; i<npart; ++i){ 

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;

  }

}


//Compute forces as -Grad_ip V(r)
double Force(int ip, int idir){ 
  double f=0.0;
  double dvec[3], dr;

  //distance ip-i in pbc
  for (int i=0; i<npart; ++i){

    if(i != ip){

      dvec[0] = Pbc( x[ip] - x[i] );  
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){

        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)

      }

    }

  }

return f;

}


//Properties measurement
void Measure(string path){ 
  int bin;
  double v, t, vij, w, wij;
  double dx, dy, dz, dr;
  double dx_g, dy_g, dz_g, dr_g;
  ofstream Epot, Ekin, Etot, Temp, Gofr, Pres;

  Epot.open(path + string("output_epot.dat"),ios::app);
  Ekin.open(path + string("output_ekin.dat"),ios::app);
  Etot.open(path + string("output_etot.dat"),ios::app);
  Temp.open(path + string("output_temp.dat"),ios::app); 
  Gofr.open(path + string("output_gofr.dat"),ios::app); 
  Pres.open(path + string("output_pres.dat"),ios::app); 

  //reset observables
  v = 0.0; 
  t = 0.0;
  w = 0.0;
  
  //reset the hystogram of g(r)
  for (int k=0; k<nbins; ++k) gofr[k]=0.0;

  //Cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){

    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );  //here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] );  //to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] );  //=> EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     
     dx_g = Pbc( x[i] - x[j] );
     dy_g = Pbc( y[i] - y[j] );
     dz_g = Pbc( z[i] - z[j] );

     dr_g = dx_g*dx_g + dy_g*dy_g + dz_g*dz_g;
     dr_g = sqrt(dr_g);

     //update of the histogram of g(r)
     bin = (int)(dr_g/bin_size);
     if(bin < nbins)
       gofr[bin] += 2.0;
     
     if(dr < rcut){

       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 4.0/pow(dr,12) - 2.0/pow(dr,6);
       v += vij;  //Potential Energy
       w += wij;  //Virial

     }

    }          

  }

  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  //Istantaneous Observables 
  stima_pot = v/(double)npart + vtail;  //Potential energy per particle
  stima_kin = t/(double)npart;  //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart;  //Temperature
  stima_etot = (t+v)/(double)npart;  //Total energy per particle
  stima_pres = rho * temp + (w + (double)npart*ptail)/vol;  //Pressure

  //Write on files
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  Pres << stima_pres << endl;
  
  for(int k=0; k<nbins; ++k)
    Gofr << k*bin_size << " " << gofr[k] << endl;
  
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Gofr.close();
  Pres.close();

}


void ConfFinal(string path){ 
  ofstream WriteConf;
  ofstream WriteOld;

  //Write final configuration
  cout << endl << "Print final configuration to file config.final" << endl;
  WriteConf.open(path + string("config.final"));
  for (int i=0; i<npart; ++i){

    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;

  }
  WriteConf.close();
  
  //Write also r(t-dt)
  //to improve the MD code by enabling  
  //the possibility to (re)start reading
  //an old spatial configuration [𝑟⃗ (𝑡−𝑑𝑡)]  
  cout << "Print r(t) and r(t-dt) for faster equilibration" << endl;
  cout << "Saving r(t) and r(t-dt) in old.0" << endl << endl;
  WriteOld.open(path + string("old.0"));
  for (int j=0; j<npart; ++j){

    WriteOld << x[j]/box << "   " <<  y[j]/box << "   " << z[j]/box << "   " <<  xold[j]/box << "   " <<  yold[j]/box << "   " << zold[j]/box << endl;

  }
  WriteOld.close();

}


/*
//Write configuration in .xyz format
void ConfXYZ(int nconf, string folder){ 
  ofstream WriteXYZ;

  WriteXYZ.open(folder + string("frames/config_") + to_string(nconf) + ".xyz");
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


//Accumulate block variables
void Accumulate(){

  sum_epot += stima_pot;
  sum_ekin += stima_kin;
  sum_etot += stima_etot;
  sum_temp += stima_temp;
  sum_pres += stima_pres;
  
  for(int k=0; k<nbins; ++k)
    sum_gofr[k] += gofr[k];
  
}


//Blocking Method
void Blocking(int blk_used, string path){
  ofstream Epot, Ekin, Etot, Temp, Gave, Pres;
  
  Epot.open(path + string("ave_epot.dat"),ios::app);
  Ekin.open(path + string("ave_ekin.dat"),ios::app);
  Etot.open(path + string("ave_etot.dat"),ios::app);
  Temp.open(path + string("ave_temp.dat"),ios::app);
  Gave.open(path + string("ave_gofr.dat"),ios::app);
  Pres.open(path + string("ave_pres.dat"),ios::app);

  //Sum of the averages in
  //each single block
  prog_epot += sum_epot/(counter*1.0);
  prog_ekin += sum_ekin/(counter*1.0);
  prog_etot += sum_etot/(counter*1.0);
  prog_temp += sum_temp/(counter*1.0);  
  prog_pres += sum_pres/(counter*1.0);
  
  //Sum of the square of the averages
  //in each single block
  prog2_epot += pow(sum_epot/(counter*1.0), 2);
  prog2_ekin += pow(sum_ekin/(counter*1.0), 2);
  prog2_etot += pow(sum_etot/(counter*1.0), 2);
  prog2_temp += pow(sum_temp/(counter*1.0), 2);
  prog2_pres += pow(sum_pres/(counter*1.0), 2);

  //Block observables averages
  ave_epot = prog_epot/(blk_used*1.0);
  ave_ekin = prog_ekin/(blk_used*1.0);
  ave_etot = prog_etot/(blk_used*1.0);
  ave_temp = prog_temp/(blk_used*1.0);
  ave_pres = prog_pres/(blk_used*1.0);
  ave2_epot = prog2_epot/(blk_used*1.0);
  ave2_ekin = prog2_ekin/(blk_used*1.0);
  ave2_etot = prog2_etot/(blk_used*1.0);
  ave2_temp = prog2_temp/(blk_used*1.0);
  ave2_pres = prog2_pres/(blk_used*1.0);
  
  //Block observables uncertainties
  if(blk_used-1==0){

      err_epot=0.0;
      err_ekin=0.0;
      err_etot=0.0;
      err_temp=0.0;
      err_pres=0.0;

  }
  else{

    err_epot=sqrt((ave2_epot-pow(ave_epot, 2))/(blk_used-1));
    err_ekin=sqrt((ave2_ekin-pow(ave_ekin, 2))/(blk_used-1));
    err_etot=sqrt((ave2_etot-pow(ave_etot, 2))/(blk_used-1));
    err_temp=sqrt((ave2_temp-pow(ave_temp, 2))/(blk_used-1));
    err_pres=sqrt((ave2_pres-pow(ave_pres, 2))/(blk_used-1));

  }
  
  //measure of g(r)
  for (int k=0; k<nbins; ++k){

    r = k * bin_size;
    stima_gofr = sum_gofr[k]/double(nstep/n_blks);
    stima_gofr *= 3.0/((4.0)*M_PI * (pow(r + bin_size,3) - pow(r,3)) * rho * (double)npart);
    ave_gofr[k] += stima_gofr;
    ave2_gofr[k] += stima_gofr*stima_gofr;
    if(blk_used == n_blks)

      Gave << scientific << r << " " << ave_gofr[k]/double(blk_used) << " " << sqrt(((ave2_gofr[k]/blk_used)-pow(ave_gofr[k]/blk_used, 2))/(blk_used-1)) << endl;

  }
  
  Epot << ave_epot << " " << err_epot << endl;
  Ekin << ave_ekin << " " << err_ekin << endl;
  Etot << ave_etot << " " << err_etot << endl;
  Temp << ave_temp << " " << err_temp << endl;
  Pres << ave_pres << " " << err_pres << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Gave.close();
  Pres.close();

}


void Delete_g(void){

  delete [] gofr;
  delete [] sum_gofr;
  delete [] ave_gofr;
  delete [] ave2_gofr;

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
