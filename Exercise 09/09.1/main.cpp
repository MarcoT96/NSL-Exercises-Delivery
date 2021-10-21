#include <iostream>
#include <fstream>
#include "random.h"
#include "GA_TSP.h"


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
unsigned int N_city;  //number of city of the TSP
unsigned int d;  //dimensionality of the TSP
unsigned int geometry;  //positions of the cities
unsigned int N_steps;  //number of steps for the
  		       //optimization via Genetic Algorithm
  		       //i.e. the number of Generations

unsigned int N_individual;  //number of individual making up a population
unsigned int generation;
double p;  //Selection exponent
double prob, prob_PP, prob_CP, prob_IN;  //call probability of mutation operators
double prob_CR;  //call probability of the crossover operator

std::vector <std::vector <double>> cities;  //vector of the positions of TSP cities
std::vector <chromosome> population;  //random starting population


/*

  Functions

*/
void readInput();
void createPopulation();
void Measure(evolution&, unsigned int);
void print(evolution&, unsigned int);
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


int main(){

 
  system("./clean.sh");
  /***************************************************************************************/
  /*Initialization*/
  readInput();
  createPopulation();
  /***************************************************************************************/

  std::cout << "Population created" << std::endl;
  evolution nature(population);
  nature.check_population();
  Measure(nature, 0);
  //print(nature, 0);
  
  /***************************************************************************************/
  /*Optimization via Genetic Algorithm*/
  
  std::cout << "Start Optimization via Genetic Algorithm" << std::endl;
  for(generation=1; generation<=N_steps; generation++){
    
    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << "Generation " << generation << std::endl;
    if(generation%100==0)
      std::cout << "*****" << std::endl;
    nature.new_generation(p, prob_PP, prob_CP, prob_IN, prob_CR);
    Measure(nature, generation);
  
  }
  
  std::cout << "Finish Optimization via Genetic Algorithm" << std::endl;
  std::cout << "Evolution up to generation " << generation-1 << std::endl;
  nature.check_population();
  std::cout <<  std::endl;
  
  /***************************************************************************************/

  //print(nature, generation-1);
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
    input_file >> string_away >> N_steps;
    input_file >> string_away >> N_individual;
    input_file >> string_away >> p;
    input_file >> string_away >> prob_PP;
    input_file >> string_away >> prob_CP;
    input_file >> string_away >> prob_IN;
    input_file >> string_away >> prob_CR;
  }
  else std::cout << "PROBLEM: Unable to open the parameters file" << std::endl;
  input_file.close();

  delete [] string_away;

}


void createPopulation(){

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
    std::cout << N_city << " cities randomly placed on a circumference" << std::endl;
  else if(geometry==1)
    std::cout << N_city << " cities randomly placed inside a square" << std::endl;
  else
    std::cerr << "!! Geometry of cities not allowed !!" << std::endl;
  std::cout << "The program solves the TSP using a Genetic Algorithm" << std::endl;
  std::cout << "Random starting population of " << N_individual << " individuals" << std::endl << std::endl;
  
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
  
  //Create the random starting population
  for(int times=0; times<N_individual; times++)
    population.push_back(chromosome(N_city, cities));

}


void Measure(evolution& nature, unsigned int gen){

  std::ofstream L, L_ave, best;
  double L_sum=0.0;
  
  L.open("L.dat", std::ios::app);
  L_ave.open("L_ave.dat", std::ios::app);
  best.open("optimized_path.dat", std::ios::app);
  
  //save the best path of each generation
  L << gen << " " << nature.get_best().get_L() << std::endl;
  best << gen << " ";
  if(nature.get_best().get_length()==0)
    best << "!! The chromosome does not contain any Genetic Material !!" << std::endl;
  else{
  
    for(unsigned int g=0; g<N_city; g++){
    
      if(g==N_city-1)
        best << nature.get_best().get_gene(g).get_allele() << "" << std::endl;
      else
        best << nature.get_best().get_gene(g).get_allele() << " ";
  
    }
  }
  
  //save L averaged on the best half of the population
  nature.restore_fitness();
  for(int best_half=0; best_half<(N_individual/2); best_half++)  
    L_sum+=nature.get_individual(best_half).get_L();
  L_ave << gen << " " << (L_sum)/((N_individual*1.0/2.0)) << std::endl;

  L.close();
  L_ave.close();
  best.close();

}


void print(evolution& nature, unsigned int generation){
  
  if(generation==0){
    
    std::cout << "Starting Population\n-----------------------------------\n" << std::endl;
    nature.print_population();
    std::cout << std::endl;
    nature.check_population();
    std::cout << "Chromosome with the best path: Individual ";
    std::cout << nature.where_is_the_individual(nature.get_best())+1;
    std::cout << "  -->  ";
    nature.get_best().print_chromosome();
    std::cout << "Chromosome with the worst path: Individual ";
    std::cout << nature.where_is_the_individual(nature.get_worst())+1;
    std::cout << "  -->  ";
    nature.get_worst().print_chromosome();
    std::cout << "\n-----------------------------------\n" << std::endl;
    
  }
  else{
  
    std::cout << "Generation " << generation << " Population\n-----------------------------------\n" << std::endl;
    nature.print_population();
    std::cout << std::endl;
    nature.check_population();
    std::cout << "Chromosome with the best path: Individual ";
    std::cout << nature.where_is_the_individual(nature.get_best())+1;
    std::cout << "  -->  ";
    nature.get_best().print_chromosome();
    std::cout << "Chromosome with the worst path: Individual ";
    std::cout << nature.where_is_the_individual(nature.get_worst())+1;
    std::cout << "  -->  ";
    nature.get_worst().print_chromosome();
    std::cout << "\n-----------------------------------\n" << std::endl;
  
  }
  
}
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
