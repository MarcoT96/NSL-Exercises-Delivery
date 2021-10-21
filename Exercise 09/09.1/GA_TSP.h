#ifndef __GA_TSP_h__
#define __GA_TSP_h__


#include <cmath>
#include <vector>
#include <algorithm>
#include <random>


/*

  The declaration of the Classes to manage the optimization
  of the TSP in two dimensions appears below.
  A brief explanation appears on the top of each of them.
  The code does not represent the maximum efficiency, at least for now,
  but it works (I hope!).

*/


/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/*
  
  The Class gene manages the elementary object of the TSP,
  i.e. the single city, characterized by an integer that identifies it
  (the 'allele' using a genetic language), an integer d that expresses
  the dimensionality of the problem and a vector position that indicates
  its location in the d-dimensional space.
  
*/
class gene{

  private:
    
    unsigned int _allele;  //the index that identifies the city	  	   
    unsigned int _d;  //dimensionality of the TSP problem
    std::vector <double> _position;  //cartesian coordinate of the city
    		  
  public:
    
    /*
    creation methods
    */
    gene() {_allele=0, _d=0; _position={};}  //trivial constructor
    gene(unsigned int city, unsigned int dim) {_allele=city, _d=dim, _position.resize(dim);}  //non-trivial constructor, no position
    gene(unsigned int, std::vector <double>);  //non-trivial constructor
    /*
      Copy Constructor
      The Copy Constructor is for creating a new object.
      It copies an existing object to a newly constructed object.
      The copy constructor is used to initialize a new instance
      from an old instance.
    */
    gene(const gene&);
    /*
      Copy Assignment Operator
      The assignment operator is to deal with an already
      existing object. The assignment operator is used to
      change an existing instance to have the same values 
      as the rvalue, which means that the instance has to
      be destroyed and re-initialized if it has internal dynamic memory.
    */
    gene& operator=(const gene&);
    //trivial destructor
    ~gene() {};  
    /*
    set methods
    */
    void set_allele(unsigned int city) {_allele=city;}
    void set_d(unsigned int dim);
    void set_position(std::vector <double>);
    /*
    get methods
    */
    unsigned int get_allele() const {return _allele;}
    unsigned int get_d() const {return _d;}
    std::vector <double> get_position() const {return _position;};
    bool operator==(const gene&) const;
    bool operator!=(const gene&) const;
    double distance(const gene&) const;  //euclidean L^2 norm
    /*
    print methods
    */
    void print_gene() const;

};
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/*

  The Class chromosome manages the single individual of the 
  population, i.e. a group of cities (the genes), whose order
  of appearance represents the order in which the salesman
  visits them.
  
  I choose to express a possible path via a 1D vector representation
  imposing the following constraints for each chromosome: the salesman
  must visit one and only one time every city and must be back to the 
  first city in the end of the path. Furthermore, I set (as suggested)
  the first city to always appear at position 1 of the vector-data_member,
  so as to reduce the degeneration of the shortest route to be 2, 
  which corresponds to walking the shortest route in clockwise and
  anti-clockwise directions.
  
  This Class is characterized by an integer number which expresses the length
  of the path, i.e. the number of cities of the TSP, a vector of genes as explained
  before, and a double-value Cost Function, which must be minimized through
  the Genetic Algorithm (GA) and represents the length of the particular path
  described by the chromosome.
  
*/
class chromosome{

  private:
    
    unsigned int _length;  //number of genes making up the chromosome
    std::vector <gene> _path;  //the route between cities
    double _L;  //Cost Function to be minimized
  
  public:
  
    /*
    creation methods
    */
    chromosome() {_length=0, _L=0; _path={};}  //trivial constructor
    chromosome(unsigned int, std::vector <std::vector <double>>);  //non-trivial constructor
    chromosome(const chromosome&);  //Copy constructor
    chromosome& operator=(const chromosome&);  //Copy Assignment Operator
    ~chromosome() {}  //trivial destructor
    /*
    set methods
    */
    void set_gene(unsigned int, const gene&);
    void set_length(unsigned int l) {_length=l;}
    void recompute_L();  //the set method for _L
    void refill(std::vector <std::vector <double>>);  //random refill of the chromosome
    /*
    get methods
    */
    unsigned int get_length() const {return _length;}
    gene get_gene(unsigned int) const;
    std::vector <gene> get_path() const {return _path;}
    double get_L() const {return _L;}
    bool operator==(const chromosome&) const;
    bool operator!=(const chromosome&) const;
    /*
    print methods
    */
    void print_chromosome() const;
    void print() const;  //complete decription of the chromosome
    bool check() const;
    /*
    mutation operator
    */
    void swap_gene(unsigned int, unsigned int);
    bool is_the_gene_inside(unsigned int, unsigned int, const gene&) const;  //true if the target gene is
    									     //inside [ind1, ind2] in the path
    unsigned int where_is_the_gene(const gene&) const;  //the position in the path of the target gene

};
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/*
  
  The Class evolution manages the Genetic Algorithm (GA).
  The methods defined below allow selection, mutation 
  and crossover operators to be applied to a population of individuals.
  
*/
class evolution{

  private:

    unsigned int _Nindividual;  //number of chromosome making up the population
    std::vector <chromosome> _p;  //The population of the paths
    chromosome _best;  //the individual with the shortest path in the population
    chromosome _worst;  //the individual with the longer path in the population
  
  public:
  
    /*
    creation methods
    */
    evolution() {_Nindividual=0, _best=chromosome(), _worst=chromosome(); _p={};}  //trivial constructor
    evolution(std::vector <chromosome>);  //non-trivial constructor
    ~evolution() {}  //destructor
    /*
    set methods
    */
    void set_individual(unsigned int, const chromosome&);
    void set_population(const std::vector <chromosome>&);  //this does not recalculate _best and _worst
    void reset_best_worst();
    /*
    get methods
    */
    unsigned int get_Nindividual() const {return _Nindividual;}
    chromosome get_best() const {return _best;}
    chromosome get_worst() const {return _worst;}
    std::vector <chromosome> get_population() const {return _p;}
    chromosome get_individual(unsigned int) const;
    /*
    print methods
    */
    void print_population() const;
    void print() const;  //complete description of the population
    void check_population() const;  //check if all the bounds of the TSP are satisfied
    /*
    selection operator
    */
    void restore_fitness();
    unsigned int selection(double);
    bool is_the_individual_inside(const chromosome&) const;  //true if the target individual is
    							     //inside the population
    unsigned int where_is_the_individual(const chromosome&) const;  //the position in the population
    								    //of the target individual
    /*
    genetic-mutation operator
    */
    void pair_permutation();  //multiple random applications of the PP
    void PP(unsigned int);  //single application of the PP
    void pair_permutation(chromosome&);  //overload for new_generation
    chromosome pair_permutation(unsigned int);  //overload for new_generation
    
    void contiguous_permutation();  //multiple random applications of the CP
    void CP(unsigned int);  //single application of the CP 
    void contiguous_permutation(chromosome&);  //overload for new_generation
    chromosome contiguous_permutation(unsigned int);  //overload for new_generation
    
    void inversion();  //multiple random applications of the IN
    void IN(unsigned int);  //single applications of the IN
    void inversion(chromosome&);  //overload for new_generation
    chromosome inversion(unsigned int);  //overload for new_generation
    
    void crossover(double);  //multiple random applications of the CR 
    void CR(unsigned int, unsigned int);  //single application of the CR
    std::vector <chromosome> crossover(unsigned int, unsigned int);  //overload for new_generation
    
    /*
    evolution process
    */
    void new_generation(double, double, double, double, double);
};
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


#endif
