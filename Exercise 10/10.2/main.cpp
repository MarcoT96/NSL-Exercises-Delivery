#include <iostream>
#include <fstream>
#include "random.h"
#include "TSP.h"
#include "mpi.h"


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
unsigned int N_migr;  //best individuals exchange period

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
void Initialize(int);
void Measure(evolution&, unsigned int, int);
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


int main(int argc, char* argv[]){

 
  system("./clean.sh");
  
  /*
    Set the MPI stage
  */
  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat1, stat2;
  MPI_Request req;
  
  /*
    Bidirectional Communications Variables
  */
  int rank1, rank2;  //the two random communicating Processes
  int itag1_allele=1;  //flag for rank1 --> rank2 communications
  int itag2_allele=2;  //flag for rank2 --> rank1 communications
  std::vector <int> imsg1_alleles;  //The messages to be exchange
  std::vector <int> imsg2_alleles;
  chromosome* imsg1;  //Supporting chromosomes
  chromosome* imsg2;
  
  /***************************************************************************************/
  /*Initialization*/
  
  readInput();
  Initialize(rank);
  if(rank==0){
  
    std::cout << "MPI Parallel GA searches of the optimal path by using " << size << " Processes (called Continents)" << std::endl;
    if(N_migr==0)
      std::cout << "Perform every GA searches without communications between the Processes" << std::endl;
    else
      std::cout << "Every " << N_migr << " generations two random Continents should exchange their best individuals randomly" << std::endl;
    std::cout << "Evolution up to generation " << N_steps << std::endl << std::endl;
  
  }
  
  /***************************************************************************************/
  
  /***************************************************************************************/
  /*Optimization via Genetic Algorithm*/
  
  //Start time
  double tstart = MPI_Wtime();
    
  /*
    Each Process initialize a random starting population
  */
  for(int times=0; times<N_individual; times++)
    population.push_back(chromosome(N_city, cities));
  evolution Continents(population);
  
  /*
    Each Process perform an independent Genetic Search
  */
  if(rank==0)
    std::cout << "Start Parallel GA Optimization" << std::endl << std::endl;
  for(generation=1; generation<=N_steps; generation++){
    
    if(N_migr!=0){
      if(generation%N_migr==0){
    
      
        /*
          Create the messages to exchange by taking only the order
          of cities from the best chromosomes of all Processes
        */
        Continents.restore_fitness();  //from now on BEST=P[0]
        imsg1_alleles.clear();
        imsg2_alleles.clear();
        for(int city=0; city<N_city; city++){
        
          imsg1_alleles.push_back(Continents.get_individual(0).get_gene(city).get_allele());
          imsg2_alleles.push_back(Continents.get_individual(0).get_gene(city).get_allele());
        
        }
        if(imsg1_alleles.size()!=N_city)
          std::cerr << "PROBLEM in the size of the message 1 to be exchanged." << std::endl;
        if(imsg2_alleles.size()!=N_city)
          std::cerr << "PROBLEM in the size of the message 2 to be exchanged." << std::endl;
      
      
        /*
          Choose the two Continents that must 
          perform the communication
        */
        rank1=rnd.Rannyu_INT(0, 3);
        rank2=rnd.Rannyu_INT(0, 3);
        while(rank1==rank2){
          rank2=rnd.Rannyu_INT(0, 3);
        }
        if(rank==0)
          std::cout << "Communication " << (int)(generation/N_migr) << "  between Process " << rank1 << " & " << rank2 << std::endl;
      
        /*
          Perform Bidirectional Communication
          between the two selected Continents
        */
        if(rank==rank1){
      
          MPI_Isend(&imsg2_alleles[0], N_city, MPI_INTEGER, rank2, itag1_allele, MPI_COMM_WORLD, &req);
          MPI_Recv(&imsg1_alleles[0], N_city, MPI_INTEGER, rank2, itag2_allele, MPI_COMM_WORLD, &stat2);
          imsg1=new chromosome(imsg1_alleles, cities);  //imsg1 is the best chromosome from the Continent rank2
          Continents.set_individual(0, *imsg1);  //Exchange best chromosomes
          Continents.restore_fitness();
          delete imsg1;  //delete the supporting chromosome
        
        
        }
        else if(rank==rank2){
      
          MPI_Send(&imsg1_alleles[0], N_city, MPI_INTEGER, rank1, itag2_allele, MPI_COMM_WORLD);
          MPI_Recv(&imsg2_alleles[0], N_city, MPI_INTEGER, rank1, itag1_allele, MPI_COMM_WORLD, &stat1);
          imsg2=new chromosome(imsg2_alleles, cities);  //imsg2 is the best chromosome from the Continent rank2
          Continents.set_individual(0, *imsg2);  //Exchange best chromosomes
          Continents.restore_fitness();
          delete imsg2;  //delete the supporting chromosome
        
        
        }
      
    
      }
    
    }
    
    
    /*
      Moving forward with evolution
    */
    Continents.new_generation(p, prob_PP, prob_CP, prob_IN, prob_CR);
    Measure(Continents, generation, rank);
   
      
  }
    
  /***************************************************************************************/
  
  //Stop time
  double tend = MPI_Wtime();
  double dt = tend - tstart;
  if(rank==0){
    
    std::cout << "\nFinish Optimization via Parallel GA Searches" << std::endl;
    std::cout << "Evolution up to generation " << generation-1 << std::endl;
  }
  std::cout << "Process " << rank << "\t Time = " << dt << " seconds" << std::endl;
  //Continents.check_population();
  
  rnd.SaveSeed();
  MPI_Finalize();

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
      input_file >> string_away >> N_migr;
    
    }
    else std::cout << "PROBLEM: Unable to open the parameters file" << std::endl;
    input_file.close();

    delete [] string_away;

}


void Initialize(int rank){

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

  if(rank==0){
    
    std::cout << "\nThe Traveling Salesman Problem (TSP)" << std::endl;
    std::cout << "Dimensionality of the TSP = " << d << std::endl;
    if(geometry==0)
      std::cout << N_city << " cities randomly placed on a circumference" << std::endl;
    else if(geometry==1)
      std::cout << N_city << " cities randomly placed inside a square" << std::endl;
    else
      std::cerr << "!! Geometry of cities not allowed !!" << std::endl;
    std::cout << "The program solves the TSP using a Genetic Algorithm" << std::endl;
    std::cout << "Random starting population of " << N_individual << " individuals" << std::endl;
  
  }
  
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
  
}


void Measure(evolution& Continents, unsigned int gen, int rank){

  std::ofstream L, L_ave, best;
  double L_sum=0.0;
  
  if(rank==0){
  
    L.open("L_rank_0.dat", std::ios::app);
    L_ave.open("L_ave_rank_0.dat", std::ios::app);
    best.open("optimized_path_rank_0.dat", std::ios::app);
  
  }
  if(rank==1){
  
    L.open("L_rank_1.dat", std::ios::app);
    L_ave.open("L_ave_rank_1.dat", std::ios::app);
    best.open("optimized_path_rank_1.dat", std::ios::app);
  
  }
  if(rank==2){
  
    L.open("L_rank_2.dat", std::ios::app);
    L_ave.open("L_ave_rank_2.dat", std::ios::app);
    best.open("optimized_path_rank_2.dat", std::ios::app);
  
  }
  if(rank==3){
  
    L.open("L_rank_3.dat", std::ios::app);
    L_ave.open("L_ave_rank_3.dat", std::ios::app);
    best.open("optimized_path_rank_3.dat", std::ios::app);
  
  }
  
  //save the best path of each generation
  L << gen << " " << Continents.get_best().get_L() << std::endl;
  best << gen << " ";
  if(Continents.get_best().get_length()==0)
    best << "!! The chromosome does not contain any Genetic Material !!" << std::endl;
  else{
  
    for(unsigned int g=0; g<N_city; g++){
    
      if(g==N_city-1)
        best << Continents.get_best().get_gene(g).get_allele() << "" << std::endl;
      else
        best << Continents.get_best().get_gene(g).get_allele() << " ";
  
    }
  }
  
  //save L averaged on the best half of the population
  Continents.restore_fitness();
  for(int best_half=0; best_half<(N_individual/2); best_half++)  
    L_sum+=Continents.get_individual(best_half).get_L();
  L_ave << gen << " " << (L_sum)/((N_individual*1.0/2.0)) << std::endl;

  L.close();
  L_ave.close();
  best.close();

}
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
