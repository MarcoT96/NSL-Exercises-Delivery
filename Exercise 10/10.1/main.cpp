#include <iostream>
#include <fstream>
#include "random.h"
#include "TSP.h"


/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/*
  Random numbers
  NY University
*/
int seed[4];
Random rnd;


/*

  Global Variable of the program

*/
//TSP Problem
unsigned int N_city;  //number of city of the TSP
unsigned int d;  //dimensionality of the TSP
unsigned int geometry;  //positions of the cities
std::vector <std::vector <double>> cities;  //vector of the positions of TSP cities

//The Annealing Schedule
unsigned int N_T;  //Length of the Annealing Schedule
double temp_step;  //amount of subtracted temperature in cooling
double T_0;  //Starting temperature of SA (high)
double T_f;  //Final temperature of SA (low)
unsigned int N_step;
//std::vector <double> temperatures={};  
//std::vector <unsigned int> n_steps={};
std::map <double, unsigned int> ann_sched={};

//SA Algorithm
double prob_PP, prob_CP, prob_IN;  //call probability of mutation operators
chromosome initial_path;


/*

  Functions

*/
void readInput();
void Initialize();
void Measure(annealer&, unsigned int);
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


int main(){

 
  system("./clean.sh");
  /***************************************************************************************/
  /*Initialization*/
  readInput();
  Initialize();
  /***************************************************************************************/

  annealer SA(initial_path, ann_sched);
  std::cout << std::endl;
  SA.print_schedule();
  std::cout << std::endl;
  
  /***************************************************************************************/
  /*Optimization via Simulated Annealing*/
  
  std::cout << "Start Optimization via Simulated Annealing" << std::endl;
  for(int t=0; t<N_T; t++){
    
    if(t%100==0)
      std::cout << "***" << std::endl; 
    SA.annealing(t, prob_PP, prob_CP, prob_IN);
    Measure(SA, t);
  
  }
  
  std::cout << "Finish Optimization via Simulated Annealing" << std::endl;
  SA.check();
  std::cout <<  std::endl;
  
  /***************************************************************************************/

  rnd.SaveSeed();


return 0;
}


/***********************************************************************************************/
/***********************************************************************************************/
/***************************************UTILITY FUNCTIONS***************************************/
/***********************************************************************************************/
/***********************************************************************************************/
void readInput(){

  std::ifstream input_file("input.dat");
  char* string_away = new char[60];

  if(input_file.is_open()){
  
    input_file >> string_away >> N_city;
    input_file >> string_away >> d;
    input_file >> string_away >> geometry;
    
    input_file >> string_away >> temp_step;
    input_file >> string_away >> T_0;
    input_file >> string_away >> T_f;
    input_file >> string_away >> N_step;
    
    input_file >> string_away >> prob_PP;
    input_file >> string_away >> prob_CP;
    input_file >> string_away >> prob_IN;
    
  }
  else std::cout << "PROBLEM: Unable to open the parameters file" << std::endl;
  input_file.close();

  delete [] string_away;

}


void Initialize(){

  /***************************************************************************************/
  /* CREATE THE RANDOM GENERATOR */
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

  std::cout << "\nThe Traveling Salesman Problem (TSP)" << std::endl;
  std::cout << "Dimensionality of the TSP = " << d << std::endl;
  if(geometry==0)
    std::cout << N_city << " cities randomly placed on a circumference" << std::endl << std::endl;
  else if(geometry==1)
    std::cout << N_city << " cities randomly placed inside a square" << std::endl << std::endl;
  else
    std::cerr << "!! Geometry of cities not allowed !!" << std::endl << std::endl;
  
  
  //Create the Annealing Schedule
  double beta=0;
  double temp=T_0;
  while(temp>T_f+temp_step){

    //temperatures.push_back(temp);
    beta=1/temp;
    ann_sched.insert(std::pair <double, unsigned int> (beta, N_step));
    temp-=temp_step;
  
  }
  ann_sched.insert(std::pair <double, unsigned int> (1/0.01, N_step));
  ann_sched.insert(std::pair <double, unsigned int> (1/0.001, N_step));
  N_T = ann_sched.size();
  
  
  std::cout << "The program solves the TSP using Simulated Annealing" << std::endl;
  std::cout << "It cool the system " << N_T << " times from T = " << T_0 << " to T = " << T_f << std::endl;
  std::cout << "Number of steps in each Metropolis sampling of the Annealing Schedule = " << N_step << std::endl << std::endl;
  
  std::ofstream outfile;
  
  /*
    Create cities randomly placed
    if geometry = 0 --> on a circumference
    if geometry = 1 --> inside a square
  */
  if(geometry==0){
  
    outfile.open("cities_on_circumference.dat");  //save on file the positions
    double r=1.0;  //radius of the circumference
    double theta=0.0;  //angle in d=2 dimension
    
    for(int j=0; j<N_city; j++){
    
      theta=rnd.Rannyu(0, 2*M_PI);
      cities.push_back({r*cos(theta), r*sin(theta)});
      outfile << "City " << j+1 << " " << r*cos(theta) << " " << r*sin(theta) << std::endl;
    
    }
    
  }
  else if(geometry==1){
  
    outfile.open("cities_in_square.dat");  //save on file the positions
    double edge=1.0;  //side of the square
    		      //centered in (0, 0)
    double x=0.0;  //coordinate of the city
    double y=0.0;
    
    for(int j=0; j<N_city; j++){
    
      x=rnd.Rannyu((-1.0*edge)/(2.0), (1.0*edge)/(2.0));
      y=rnd.Rannyu((-1.0*edge)/(2.0), (1.0*edge)/(2.0));
      cities.push_back({x, y});
      outfile << "City " << j+1 << " " << x << " " << y << std::endl;
    
    }
  
  }
  else std::cerr << "Error!\nGeometry of the TSP not allowed!" << std::endl;
  outfile.close();
  std::cout << "Create the Annealing Schedule" << std::endl;
  
  //Create the random starting path
  std::cout << "Create a random starting path between the cities" << std::endl;
  initial_path=chromosome(cities.size(), cities);

}


void Measure(annealer& SA, unsigned int index_temp){

  std::ofstream L, best;
  
  L.open("L.dat", std::ios::app);
  best.open("optimized_path.dat", std::ios::app);
  
  //save the best path of each generation
  L << std::showpoint << std::setprecision(3) << SA.get_temperature(index_temp) << " ";
  L << std::setprecision(4) << SA.get_path().get_L() << " ";
  L << std::showpoint << std::setprecision(4) << SA.get_acceptance() << std::endl;
  best << std::showpoint << std::setprecision(3) << SA.get_temperature(index_temp) << " ";
  if(SA.get_path().get_length()==0)
    best << "!! The chromosome does not contain any Genetic Material !!" << std::endl;
  else{
  
    for(unsigned int g=0; g<N_city; g++){
    
      if(g==N_city-1)
        best << SA.get_path().get_gene(g).get_allele() << "" << std::endl;
      else
        best << SA.get_path().get_gene(g).get_allele() << " ";
  
    }
  }

  L.close();
  best.close();

}
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
