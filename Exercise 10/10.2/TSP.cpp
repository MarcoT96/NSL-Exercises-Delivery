#include <iostream>
#include "TSP.h"


/***********************************************************************************************/
/***********************************************************************************************/
/*******************************************CLASS GENE******************************************/
/***********************************************************************************************/
/***********************************************************************************************/
gene::gene(unsigned int city, std::vector <double> pos){

  _allele=city;
  _d=pos.size();
  
  for(auto x : pos)
    _position.push_back(x);

}


//Copy Constructor
gene::gene(const gene& copy){

  _allele=copy.get_allele();
  _d=copy.get_d();
  _position=copy.get_position();
  
  /*
  _position.clear();
  _position.resize(_d);
  for(int dim=0; dim<_d; dim++){
    _position[dim]=copy.get_position()[dim];
  }
  */

}


//Copy Assignment
gene& gene::operator=(const gene& copy){

  _allele=copy.get_allele();
  _d=copy.get_d();
  
  //delete memory
  _position.clear();
  for(auto pos : copy.get_position())
    _position.push_back(pos);
    
  return *this;

}


void gene::set_d(unsigned int dim){

  _d=dim;
  
  _position.clear();
  _position.resize(dim);

}


void gene::set_position(std::vector <double> pos){

  if(_d != pos.size())
    _d=pos.size();
  
  _position.clear();
  for(auto new_pos : pos)
    _position.push_back(new_pos);

}


bool gene::operator==(const gene& rhs) const{

  if(_allele==rhs.get_allele() && _d==rhs.get_d() && _position==rhs.get_position())
    return true;
  else
    return false;

}


bool gene::operator!=(const gene& rhs) const{

  if(_allele!=rhs.get_allele() || _d!=rhs.get_d() || _position!=rhs.get_position())
    return true;
  else
    return false;

}


double gene::distance(const gene& g) const{

  double dist=0;

  if(_d==g.get_d()){
    for(int j=0; j<_d; j++){
      dist+=pow(_position[j]-g.get_position()[j], 2);
    }
  }
  else{
    std::cerr << "Error: the two genes are incompatible!\nThey have different dimension!" << std::endl;
  }
  return dist;

}


void gene::print_gene() const{

  std::cout << "City " << get_allele() << " with cartesian coordinates (";
  for(int j=0; j<get_d(); j++){
  
    if(j==get_d()-1)
      std::cout << get_position()[j] << ")" << std::endl;
    else
      std::cout << get_position()[j] << ", ";
  }

}
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


/***********************************************************************************************/
/***********************************************************************************************/
/***************************************CLASS CHROMOSOME****************************************/
/***********************************************************************************************/
/***********************************************************************************************/
chromosome::chromosome(unsigned int l, std::vector <std::vector <double>> positions){

  _length=l;
  int r=0;
  double L=0.0;
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(2, _length);
  
  
  //set the first city to always appear at position 1
  _path.push_back(gene(1, positions[0]));
  
  //create a path between the cities randomly
  for(int j=0; j<(_length-1); j++){
    
    r=rnd(gen);
    while(std::find(_path.begin(), _path.end(), gene(r, positions[r-1])) != _path.end()){  //in this case std::find() returns
      r=rnd(gen);			   						   //true if gene(r) is already present 
    }  									   		   //in the chromosome
    
    _path.push_back(gene(r, positions[r-1]));
  
  }
  
  //calculate the length of the path using 
  //euclidean L^2 norm
  for(auto it=_path.begin(); it<_path.end(); it++){
    
    //the salesman must be back to the 
    //first city in the end of the path
    if(it==_path.end()-1)
      L+=(*it).distance(_path[0]);
    else
      L+=(*it).distance(*(it+1));
      
  }
  _L=L;

}


chromosome::chromosome(std::vector <int> alleles, std::vector <std::vector <double>> positions){

  /*
    Explicitly pass the sequence of integers that defines the path
    This is convenient in the procedure of MPI Bidirectional Communication
  */

  _length=alleles.size();
  double L=0.0;

  //check the first city to always appear at position 1
  if(alleles[0]!=1)
    std::cerr << " !! Error: the first city must always appear at position 1 !!" << std::endl;
  //check dimensions matching
  if(alleles.size()!=positions.size())
    std::cerr << " PROBLEM: the size of the list of cities labels does not match with the size of cities positions." << std::endl;

  //create the argument path between the cities
  for(int j=0; j<_length; j++)
    _path.push_back(gene(alleles[j], positions[alleles[j]-1]));

  //calculate the length of the path using
  //euclidean L^2 norm
  for(auto it=_path.begin(); it<_path.end(); it++){

    //the salesman must be back to the
    //first city in the end of the path
    if(it==_path.end()-1)
      L+=(*it).distance(_path[0]);
    else
      L+=(*it).distance(*(it+1));

  }
  _L=L;

}


//Copy Constructor
chromosome::chromosome(const chromosome& copy){

  _length=copy.get_length();
  _L=copy.get_L();
  _path=copy.get_path();
  /*
  _path.clear();
  _path.resize(_length);
  for(int g=0; g<_length; g++){
    _path[g]=copy.get_path()[g];
  }
  */
}


//Copy Assignment
chromosome& chromosome::operator=(const chromosome& copy){

  _length=copy.get_length();
  _L=copy.get_L();  //useless to recompute it
  		    //I'm making a copy
  		    
  //delete memory
  _path.clear();
  for(auto gene : copy.get_path())
    _path.push_back(gene);

  return *this;
  
}


void chromosome::set_gene(unsigned int index, const gene& g){

  if(index >= _length)
    std::cerr << "Error!\nIndex out of range!" << std::endl;
  else{
    //Copy Assignment of Class gene
    _path[index]=g;
  }

}


void chromosome::recompute_L(){

  double L=0.0;
  
  for(auto it=_path.begin(); it<_path.end(); it++){
    
    //the salesman must be back to the 
    //first city in the end of the path
    if(it==_path.end()-1)
      L+=(*it).distance(_path[0]);
    else
      L+=(*it).distance(*(it+1));
      
  }
  _L=L;

}


void chromosome::refill(std::vector <std::vector <double>> positions){
  
  _length=positions.size();
  int r=0;
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(2, _length);
  
  
  _path.clear();
  //set the first city to always appear at position 1
  _path.push_back(gene(1, positions[0]));
  
  //create a path between the cities randomly
  for(int j=1; j<_length; j++){
    
    r=rnd(gen);
    while(std::find(_path.begin(), _path.end(), gene(r, positions[r-1])) != _path.end()){  //in this case std::find() returns
      r=rnd(gen);							   		   //true if gene(r) is already present 
    }									   		   //in the chromosome
    
    _path.push_back(gene(r, positions[r-1]));
  }
  
  //re-calculate the length of the path using 
  //euclidean L^2 norm
  recompute_L();
  
}


gene chromosome::get_gene(unsigned int g) const{

  if(g >= get_length()){
    std::cerr << "Error!\nThe element " << g+1 << " does not exist in the cities path!" << std::endl;
    return gene();
  }
  else
    return _path[g];

}


bool chromosome::operator==(const chromosome& rhs) const{

  if(_length==rhs.get_length() && _L==rhs.get_L() && _path==rhs.get_path())
    return true;
  else
    return false;

}


bool chromosome::operator!=(const chromosome& rhs) const{

  if(_length!=rhs.get_length() || _L!=rhs.get_L() || _path!=rhs.get_path())
    return true;
  else
    return false;

}


void chromosome::print_chromosome() const{
  
  if(_length==0)
    std::cerr << "!! The chromosome does not contain any Genetic Material !!" << std::endl;
  else{
  
    std::cout << "[";
  
    for(auto it=_path.begin(); it<_path.end(); it++){
    
      if(it==_path.end()-1)
        std::cout << (*it).get_allele() << "]" << std::endl;
      else
        std::cout << (*it).get_allele() << ", ";
  
    }
  }

}


void chromosome::print() const{

  if(_length==0)
    std::cerr << "!! The chromosome does not contain any Genetic Material !!" << std::endl;
  else{

    for(auto g : _path)
      g.print_gene();
      
  }
    
}


bool chromosome::check() const{

  unsigned int equals=0;
  unsigned int ok=0;
  unsigned int city_not_allowed=0;

  if(_length==0){
    
    std::cerr << "!! The chromosome does not contain any Genetic Material !!" << std::endl;
    return false;
  
  }
  else{
  
    for(int j=0; j<_length; j++){
    
      equals=0;
      for(int k=j+1; k<_length; k++){
      
        if(_path[j]==_path[k])
          equals++;
      
      }
      
      if(equals==0)
        ok++;
        
      if(_path[j].get_allele()<1 || _path[j].get_allele()>_length)
        city_not_allowed++;
    
    
    }
    
    if(ok==_length && _path[0].get_allele()==1 && city_not_allowed==0)
      return true;
    else
      return false;
      
  }

}


void chromosome::swap_gene(unsigned int g1, unsigned int g2){

  if(g1!=0 && g2!=0 && g1<_length && g2<_length){
  
    gene temp=_path[g1];  //Copy Constructor
    _path[g1]=_path[g2];  //Copy Assignment
    _path[g2]=temp;  //Copy Assignment
    
  }
  
  //re-calculate the length of the path using 
  //euclidean L^2 norm
  //after pair permutation
  recompute_L();
  
}


bool chromosome::is_the_gene_inside(unsigned int ind1, unsigned int ind2, const gene& target_gene) const{
  
  unsigned int found=0;					 
  for(int index_gene=ind1; index_gene<=ind2; index_gene++){
  
    if(_path[index_gene]==target_gene)
      found++;
  
  }
  
  if(found != 0)
    return true;
  else
    return false;

}


unsigned int chromosome::where_is_the_gene(const gene& target_gene) const{

  unsigned int target_index=0;
  if(is_the_gene_inside(0, _length-1, target_gene)){
  
    auto target_it=std::find(_path.begin(), _path.end(), target_gene);
    target_index=std::distance(_path.begin(), target_it);
    return target_index;
  
  }
  else{
  
    std::cerr << "!! Gene not present inside the chromosome !!" << std::endl;
    return -1;
  
  }

}
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


/***********************************************************************************************/
/***********************************************************************************************/
/****************************************CLASS EVOLUTION****************************************/
/***********************************************************************************************/
/***********************************************************************************************/
evolution::evolution(std::vector <chromosome> population){
  
  _Nindividual=population.size();
  
  for(auto c : population)
    _p.push_back(c);

  auto best=std::min_element(_p.begin(), _p.end(), 
  			     [&](const chromosome& lhs, const chromosome& rhs){
    				 return lhs.get_L()<rhs.get_L();
    			     }
    				);
  auto worst=std::max_element(_p.begin(), _p.end(), 
  			      [&](const chromosome& lhs, const chromosome& rhs){
    				  return lhs.get_L()<rhs.get_L();
    			      }
    				);
  _best=*best;
  _worst=*worst;
    				
}


void evolution::set_individual(unsigned int index, const chromosome& c){

  if(index >= _Nindividual)
    std::cerr << "Error!\nIndex out of range!" << std::endl;
  else{
    //Copy Assignment of Class chromosome
    _p[index]=c;
  }

}


void evolution::set_population(const std::vector <chromosome>& population){

  _Nindividual=population.size();

  _p.clear();
  for(auto c : population)
    _p.push_back(c);

}


void evolution::reset_best_worst(){

  auto best=std::min_element(_p.begin(), _p.end(), 
  			     [&](const chromosome& lhs, const chromosome& rhs){
    				 return lhs.get_L()<rhs.get_L();
    			     }
    				);
  
  auto worst=std::max_element(_p.begin(), _p.end(), 
  			      [&](const chromosome& lhs, const chromosome& rhs){
    				  return lhs.get_L()<rhs.get_L();
    			      }
    				);
    				
  _best=*best;
  _worst=*worst;

}


chromosome evolution::get_individual(unsigned int index) const{

  if(index >= _Nindividual){
    std::cerr << "Error!\nThe Individual " << index+1 << " does not exist in the population!" << std::endl;
    return chromosome();
  }
  else
    return _p[index];

}


void evolution::print_population() const{

  for(int j=1; j<=_Nindividual; j++){
    std::cout << "=======================" << std::endl;
    std::cout << "Individual " << j << std::endl;
    std::cout << "=======================" << std::endl;
    _p[j-1].print_chromosome();
    std::cout << "Length of the path between the cities: " << _p[j-1].get_L();
    std::cout << std::endl;
  }

}


void evolution::print() const{

  for(int j=1; j<=_Nindividual; j++){
    std::cout << "=======================" << std::endl;
    std::cout << "Individual " << j << std::endl;
    std::cout << "=======================" << std::endl;
    _p[j-1].print();
    std::cout << "Length of the path between the cities: " << _p[j-1].get_L();
    std::cout << std::endl;
  }

}


void evolution::check_population() const{

  int good=0;
  for(int j=0; j<_Nindividual; j++){
    
    if(_p[j].check())
      good++;
      
  }
  
  if(good==_Nindividual)
    std::cout <<"The population fulfils all the bonds imposed by the TSP" << std::endl;
  else
    std::cerr <<"Error!\nThe population does not fulfill all the bonds imposed by the TSP" << std::endl;
    
}


void evolution::restore_fitness(){
  
  //order the population based on fitness
  std::sort(_p.begin(), _p.end(),
  	    [&](const chromosome& lhs, const chromosome& rhs){
    		return lhs.get_L()<rhs.get_L();
    	    }
    		);
  
  //Restore best and worst individual
  reset_best_worst();

}


/*
  This is the Roulette Wheel Selection
  It may be useful to try to implement 
  another selection method...
  N.B. Remember to order the population FIRST! 
*/
unsigned int evolution::selection(double p){
  
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> rnd(0, 1);
  
  //order the population based on fitness
  //restore_fitness();
  
  //select an individual based on fitness
  return int(_Nindividual*pow(rnd(gen), p)) + 1;
  
}


/*
  Apply mutation to a random group of chromosomes
  consecutively, choosing at random the genes
  to be exchanged
*/
void evolution::pair_permutation(){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd1(1, _Nindividual);  //to how many individuals
  						          //of the population to apply
  						          //the pair permutation, even to
  						          //the same (! remember closed interval [a, b] !)
  std::uniform_int_distribution<> rnd2(1, _p[0].get_length()-1);  //which genes to exchange
  								  //in the chosen chromosome
  unsigned int how_many=rnd1(gen);
  unsigned int choosen;
  unsigned int g1, g2;
  
  
  for(int times=0; times<how_many; times++){
    
    choosen=rnd1(gen)-1;  //choose the individual random
    g1=rnd2(gen);  //choose the two 
    g2=rnd2(gen);  //genes to swap random
    
    while(g1==g2)
      g2=rnd2(gen);
    
    if(is_the_individual_inside(_p[choosen])){
    
      _p[choosen].swap_gene(g1, g2);  //this already recalculates the fitness
				      //of the mutated chromosome
    }
    else std::cerr << "!! Individual not present inside the chromosome !!" << std::endl;
    
  }  
  
  //After the mutation restore
  //the fitness of the population
  reset_best_worst();

}


/*
  Single application of the previous operator
  to the given CHOOSEN chromosome
*/
void evolution::PP(unsigned int choosen){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(1, _p[0].get_length()-1);  //which genes to exchange in the chosen chromosome
  unsigned int g1, g2;
  
 
  g1=rnd(gen);  //choose the two 
  g2=rnd(gen);  //genes to swap random
    
  while(g1==g2)
    g2=rnd(gen);
    
  if(is_the_individual_inside(_p[choosen])){
    
      _p[choosen].swap_gene(g1, g2);  //this already recalculates the fitness
				      //of the mutated chromosome
  }
  else std::cerr << "!! Individual not present inside the chromosome !!" << std::endl;

}


void evolution::pair_permutation(chromosome& mutated){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(1, mutated.get_length()-1);  //which genes to exchange in the chosen chromosome
  unsigned int g1, g2;
  
 
  g1=rnd(gen);  //choose the two 
  g2=rnd(gen);  //genes to swap random
    
  while(g1==g2)
    g2=rnd(gen);
   
  mutated.swap_gene(g1, g2);  //this already recalculates the fitness
			      //of the mutated chromosome

}


chromosome evolution::pair_permutation(unsigned int choosen){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(1, _p[0].get_length()-1);  //which genes to exchange in the chosen chromosome
  unsigned int g1, g2;
  chromosome mutated=_p[choosen];
  
 
  g1=rnd(gen);  //choose the two 
  g2=rnd(gen);  //genes to swap random
    
  while(g1==g2)
    g2=rnd(gen);
    
  if(is_the_individual_inside(mutated)){
    
      mutated.swap_gene(g1, g2);  //this already recalculates the fitness
				      //of the mutated chromosome
  }
  else std::cerr << "!! Individual not present inside the chromosome !!" << std::endl;

  return mutated;

}


bool evolution::is_the_individual_inside(const chromosome& target_ind) const{

  auto target_it=std::find(_p.begin(), _p.end(), target_ind);  //in this case std::find() returns
							       //true if target_gene is present in [ind1, ind2]
							       //Remember: in std::find last is not taken

  if(target_it != _p.end())
    return true;
  else
    return false;

}


unsigned int evolution::where_is_the_individual(const chromosome& target_ind) const{

  unsigned int target_index=0;
  if(is_the_individual_inside(target_ind)){
  
    auto target_it=std::find(_p.begin(), _p.end(), target_ind);
    target_index=std::distance(_p.begin(), target_it);
    return target_index;
  
  }
  else{
  
    std::cerr << "!! Individual not present inside the chromosome !!" << std::endl;
    return -1;
  
  }

}


void evolution::contiguous_permutation(){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd1(1, _Nindividual);  //to how many individuals
  						          //of the population to apply
  						          //the contiguous permutation, even to
  						          //the same
  std::uniform_int_distribution<> rnd2(1, _p[0].get_length()/2);  //the contiguous cities in the chosen chromosome
  std::uniform_int_distribution<> rnd3(1, _p[0].get_length()-1);				       
  unsigned int how_many=rnd1(gen);
  unsigned int choosen=0, start_gene1=0, start_gene2=0, m=0;
  

  for(int times=0; times<how_many; times++){
  
    choosen=rnd1(gen)-1;  //the choosen chromosome
    start_gene1=rnd3(gen);  //the 1st-first of the 1st m-contiguous cities
    start_gene2=rnd3(gen);  //the 2nd-first of the 2nd m-contiguous cities
    m=rnd2(gen);  //the number of contiguous cities
    
    while(start_gene1+m >= _p[choosen].get_length() || start_gene2+m >= _p[choosen].get_length() || start_gene2 <= start_gene1+m){
      
      m=rnd2(gen);
      start_gene1=rnd3(gen);
      start_gene2=rnd3(gen);
    
    }
    
    if(is_the_individual_inside(_p[choosen])){
    
      for(int count=0; count<m; count++)
        _p[choosen].swap_gene(start_gene1+count, start_gene2+count);
    
    }
  
  }
  
  //Restore best and worst individual
  reset_best_worst();
  
}


void evolution::CP(unsigned int choosen){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd1(1, _p[0].get_length()/2);  //the contiguous cities in the chosen chromosome
  std::uniform_int_distribution<> rnd2(1, _p[0].get_length()-1);				       
  unsigned int start_gene1=0, start_gene2=0, m=0;
  

  start_gene1=rnd2(gen);  //the 1st-first of the 1st m-contiguous cities
  start_gene2=rnd2(gen);  //the 2nd-first of the 2nd m-contiguous cities
  m=rnd1(gen);  //the number of contiguous cities
    
  while(start_gene1+m >= _p[choosen].get_length() || start_gene2+m >= _p[choosen].get_length() || start_gene2 <= start_gene1+m){
      
    m=rnd1(gen);
    start_gene1=rnd2(gen);
    start_gene2=rnd2(gen);
    
  }
    
  if(is_the_individual_inside(_p[choosen])){
    
    for(int count=0; count<m; count++)
      _p[choosen].swap_gene(start_gene1+count, start_gene2+count);
    
  }
  
}


void evolution::contiguous_permutation(chromosome& mutated){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd1(1, mutated.get_length()/2);  //the contiguous cities in the chosen chromosome
  std::uniform_int_distribution<> rnd2(1, mutated.get_length()-1);				       
  unsigned int start_gene1=0, start_gene2=0, m=0;
  

  start_gene1=rnd2(gen);  //the 1st-first of the 1st m-contiguous cities
  start_gene2=rnd2(gen);  //the 2nd-first of the 2nd m-contiguous cities
  m=rnd1(gen);  //the number of contiguous cities
    
  while(start_gene1+m >= mutated.get_length() || start_gene2+m >= mutated.get_length() || start_gene2 <= start_gene1+m){
      
    m=rnd1(gen);
    start_gene1=rnd2(gen);
    start_gene2=rnd2(gen);
    
  }
    
  for(int count=0; count<m; count++)
    mutated.swap_gene(start_gene1+count, start_gene2+count);
  
}


chromosome evolution::contiguous_permutation(unsigned int choosen){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd1(1, _p[0].get_length()/2);  //the contiguous cities in the chosen chromosome
  std::uniform_int_distribution<> rnd2(1, _p[0].get_length()-1);				       
  unsigned int start_gene1=0, start_gene2=0, m=0;
  chromosome mutated=_p[choosen];
  

  start_gene1=rnd2(gen);  //the 1st-first of the 1st m-contiguous cities
  start_gene2=rnd2(gen);  //the 2nd-first of the 2nd m-contiguous cities
  m=rnd1(gen);  //the number of contiguous cities
    
  while(start_gene1+m >= mutated.get_length() || start_gene2+m >= mutated.get_length() || start_gene2 <= start_gene1+m){
      
    m=rnd1(gen);
    start_gene1=rnd2(gen);
    start_gene2=rnd2(gen);
    
  }
    
  if(is_the_individual_inside(mutated)){
    
    for(int count=0; count<m; count++)
      mutated.swap_gene(start_gene1+count, start_gene2+count);
    
  }
  
  return mutated;
  
}


void evolution::inversion(){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd1(1, _Nindividual);  //to how many individuals
  						          //of the population to apply
  						          //the contiguous shift, even to
  						          //the same
  std::uniform_int_distribution<> rnd2(1, _p[0].get_length()-1);  //the contiguous cities in the chosen chromosome					       
  unsigned int how_many=rnd1(gen);
  unsigned int choosen=0, start_gene=0, m=0, start=0, end=0;
  
  
  for(int times=0; times<how_many; times++){
  
    choosen=rnd1(gen)-1;  //the choosen chromosome
    start_gene=rnd2(gen);  //the first of the m-contiguous cities
    m=rnd2(gen);  //the number of contiguous cities
    
    while(start_gene+m >= _p[choosen].get_length()){
    
      m=rnd2(gen);
      start_gene=rnd2(gen);
      
    }
    
    if(is_the_individual_inside(_p[choosen])){
      		          
      start=start_gene;
      end=start_gene+m-1;
    
      while(start<end){
    
        _p[choosen].swap_gene(start, end);
        start++;
        end--;
    
      }
      
    }  
    
  }
  
  //Restore best and worst individual
  reset_best_worst();

}


void evolution::IN(unsigned int choosen){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(1, _p[0].get_length()-1);  //the contiguous cities in the chosen chromosome					       
  unsigned int start_gene=0, m=0, start=0, end=0;
  

  start_gene=rnd(gen);  //the first of the m-contiguous cities
  m=rnd(gen);  //the number of contiguous cities
    
  while(start_gene+m >= _p[choosen].get_length()){
    
    m=rnd(gen);
    start_gene=rnd(gen);
      
  }
    
  if(is_the_individual_inside(_p[choosen])){
      		          
    start=start_gene;
    end=start_gene+m-1;
    
    while(start<end){
    
      _p[choosen].swap_gene(start, end);
      start++;
      end--;
    
    }
      
  }  

}


void evolution::inversion(chromosome& mutated){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(1, mutated.get_length()-1);  //the contiguous cities in the chosen chromosome					       
  unsigned int start_gene=0, m=0, start=0, end=0;
  

  start_gene=rnd(gen);  //the first of the m-contiguous cities
  m=rnd(gen);  //the number of contiguous cities
    
  while(start_gene+m >= mutated.get_length()){
    
    m=rnd(gen);
    start_gene=rnd(gen);
      
  }
    
  start=start_gene;
  end=start_gene+m-1;
    
  while(start<end){
    
    mutated.swap_gene(start, end);
    start++;
    end--;
    
  }
   
}


chromosome evolution::inversion(unsigned int choosen){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rnd(1, _p[0].get_length()-1);  //the contiguous cities in the chosen chromosome					       
  unsigned int start_gene=0, m=0, start=0, end=0;
  chromosome mutated=_p[choosen];
  

  start_gene=rnd(gen);  //the first of the m-contiguous cities
  m=rnd(gen);  //the number of contiguous cities
    
  while(start_gene+m >= mutated.get_length()){
    
    m=rnd(gen);
    start_gene=rnd(gen);
      
  }
    
  if(is_the_individual_inside(mutated)){
      		          
    start=start_gene;
    end=start_gene+m-1;
    
    while(start<end){
    
      mutated.swap_gene(start, end);
      start++;
      end--;
    
    }
      
  } 
  
  return mutated; 

}


void evolution::crossover(double p){

  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd(1, _p[0].get_length()-1);
  std::uniform_int_distribution<> rnd1(1, _Nindividual);
  unsigned int how_many=rnd1(gen);
  unsigned int index_mother=0, index_father=0, index_cut=0;
  unsigned int N_mutations=0, mutation=0, bond=0;
  chromosome* temp;
  
  
  //selection of the Genetic heritage
  bond=_p[0].get_length();
  restore_fitness();
  for(int times=0; times<how_many; times++){
    index_mother=selection(p)-1;
    index_father=selection(p)-1;
    temp=new chromosome(_p[index_mother]);
  
    //Recombine the Genetic material
    index_cut=rnd(gen);  //cut the path of mother and father random in the same point
    N_mutations=bond-index_cut;  //number of total mutations
 
    /*
      Complete the cut paths of the mother (son1) with the missing
      cities adding them in the order in which they appear in the
      consort selected father
    */
    mutation=0;
    for(int j=0;  j<bond; j++){   
  
      if(_p[index_mother].is_the_gene_inside(0, index_cut-1, _p[index_father].get_gene(j))==false){

        _p[index_mother].set_gene(index_cut+mutation, _p[index_father].get_gene(j));
        mutation++;
        
      }
  
    }
  
    /*
      Complete the cut paths of the father (son2) with the missing
      cities adding them in the order in which they appear in the
      consort selected mother
    */
    mutation=0;
    for(int j=0;  j<bond; j++){   
      
      if(_p[index_father].is_the_gene_inside(0, index_cut-1, temp -> get_gene(j)) ==  false){
    
        _p[index_father].set_gene(index_cut+mutation, temp -> get_gene(j));
        mutation++;
        
      }
      
    }
    delete temp;
    
  }
  
  _p[index_mother].recompute_L();
  _p[index_father].recompute_L();
  reset_best_worst();
    
}


void evolution::CR(unsigned int mother, unsigned int father){

  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd(1, _p[0].get_length()-1);
  unsigned int index_cut=0, N_mutations=0, mutation=0, bond=0;
  chromosome* temp;
  

  bond=_p[0].get_length();
  temp=new chromosome(_p[mother]);
  
  //Recombine the Genetic material
  index_cut=rnd(gen);  //cut the path of mother and father random in the same point
  N_mutations=bond-index_cut;  //number of total mutations
 
  /*
    Complete the cut paths of the mother (son1) with the missing
    cities adding them in the order in which they appear in the
    consort selected father
  */
  mutation=0;
  for(int j=0;  j<bond; j++){   
  
    if(_p[mother].is_the_gene_inside(0, index_cut-1, _p[father].get_gene(j))==false){

        _p[mother].set_gene(index_cut+mutation, _p[father].get_gene(j));
        mutation++;
        
    }
  
  }
  
  /*
    Complete the cut paths of the father (son2) with the missing
    cities adding them in the order in which they appear in the
    consort selected mother
  */
  mutation=0;
  for(int j=0;  j<bond; j++){   
      
    if(_p[father].is_the_gene_inside(0, index_cut-1, temp -> get_gene(j)) ==  false){
    
      _p[father].set_gene(index_cut+mutation, temp -> get_gene(j));
      mutation++;
        
    }
      
  }
  
  _p[mother].recompute_L();
  _p[father].recompute_L();
  delete temp;  
    
}


std::vector <chromosome> evolution::crossover(unsigned int mother, unsigned int father){

  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd(1, _p[0].get_length()-1);
  unsigned int index_cut=0, N_mutations=0, mutation=0, bond=0;
  chromosome* temp;
  std::vector <chromosome> offspring={};
  

  bond=_p[mother].get_length();
  temp=new chromosome(_p[mother]);
  offspring.push_back(_p[mother]);
  offspring.push_back(_p[father]);
  
  //Recombine the Genetic material
  index_cut=rnd(gen);  //cut the path of mother and father random in the same point
  N_mutations=bond-index_cut;  //number of total mutations
 
  /*
    Complete the cut paths of the mother (son1) with the missing
    cities adding them in the order in which they appear in the
    consort selected father
  */
  mutation=0;
  for(int j=0; j<bond; j++){   
  
    if(offspring[0].is_the_gene_inside(0, index_cut-1, offspring[1].get_gene(j))==false){

        offspring[0].set_gene(index_cut+mutation, offspring[1].get_gene(j));
        mutation++;
        
    }
  
  }
  
  /*
    Complete the cut paths of the father (son2) with the missing
    cities adding them in the order in which they appear in the
    consort selected mother
  */
  mutation=0;
  for(int j=0;  j<bond; j++){   
      
    if(offspring[1].is_the_gene_inside(0, index_cut-1, temp -> get_gene(j)) ==  false){
    
      offspring[1].set_gene(index_cut+mutation, temp -> get_gene(j));
      mutation++;
        
    }
      
  }
  delete temp;
  
  offspring[0].recompute_L();
  offspring[1].recompute_L();
  return offspring;  
    
}


/*
  This method implements what happens
  during the evolution of a single generation
*/
void evolution::new_generation(double p, double prob_PP, double prob_CP, double prob_IN, double prob_CR){

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> rnd(0, 1);
  std::vector <chromosome> new_population={};
  std::vector <chromosome> offspring={};
  double PP=0, CP=0, IN=0, CR=0; 
  unsigned int parent1, parent2;  //more politically correct than
  				  //mother and father
  
  
  restore_fitness();  //Order the population based
  		      //on chromosomes fitness
  new_population.clear();
  for(int evlt=0; evlt<int(_Nindividual/2); evlt++){
  
    //I want to select two chromosomes each time 
    //to get a new _Nindividual offspring paths
    offspring.clear();
    parent1=selection(p)-1;
    parent2=selection(p)-1;
    
    //I have children with a certain probability
    if(rnd(gen)<prob_CR){
      
      offspring=crossover(parent1, parent2);
      CR++;
    
    }
    //or not
    else{
    
      offspring.push_back(_p[parent1]);
      offspring.push_back(_p[parent2]);
    
    }
    
    //I apply certain mutations to children
    //always with a certain probability
    if(rnd(gen)<prob_PP){
    
      pair_permutation(offspring[0]);
      pair_permutation(offspring[1]);
      PP++;
    
    }
    if(rnd(gen)<prob_CP){
    
      contiguous_permutation(offspring[0]);
      contiguous_permutation(offspring[1]);
      CP++;
    
    }
    if(rnd(gen)<prob_IN){
    
      inversion(offspring[0]);
      inversion(offspring[1]);
      IN++;
    
    }
    
    new_population.push_back(offspring[0]);
    new_population.push_back(offspring[1]);
  
  }
  
  if(new_population.size()==_Nindividual){
    //replace the new population
    _p=new_population;
    reset_best_worst();
    //restore_fitness();
    //check_population();
  }
  else std::cerr << "!! Error during the evolution phase !!" << std::endl;

  /*
  std::cout << "Crossover p_c = " << (CR*100.0)/(_Nindividual/2) << "%" << std::endl;
  std::cout << "Pair Permutation p_m = " << (PP*100.0)/(_Nindividual/2) << "%" << std::endl;
  std::cout << "Contiguous Permutation p_m = " << (CP*100.0)/(_Nindividual/2) << "%" << std::endl;
  std::cout << "Inversion p_m = " << (IN*100.0)/(_Nindividual/2) << "%" << std::endl << std::endl;
  */

}
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


/***********************************************************************************************/
/***********************************************************************************************/
/****************************************CLASS ANNEALER*****************************************/
/***********************************************************************************************/
/***********************************************************************************************/
annealer::annealer(chromosome path, std::map <double, unsigned int> s){

  //Copy Assignment Operator of the
  //Class chromosome called
  _p=path;
  
  //Copy Assignment Operator of the
  //built-in Class map
  _schedule=s;
  _acceptance=0;

}


annealer::annealer(chromosome path, std::vector <double> beta, std::vector <unsigned int> steps){

  //Copy Assignment Operator of the
  //Class chromosome called
  _p=path;
  
  if(beta.size()==steps.size()){
  
    for(int j=0; j<beta.size(); j++)
      _schedule[beta[j]] = steps[j];
  
  }
  else{ 
  
    std::cerr << "!! Error in the creation of the Annealing Schedule !!" << std::endl;
    std::cout << "Problem with dimensions which do not match" << std::endl;
  
  }
  
  _acceptance=0;

}


void annealer::set_schedule(const std::vector <double>& beta, const std::vector <unsigned int>& steps){

  _schedule.clear();
  if(beta.size()==steps.size()){
  
    for(int j=0; j<beta.size(); j++)
      _schedule[beta[j]] = steps[j];
  
  }
  else{ 
  
    std::cerr << "!! Error in the assignment of the Annealing Schedule !!" << std::endl;
    std::cout << "Problem with dimensions which do not match" << std::endl;
  
  }

}


double annealer::get_temperature(unsigned int index_schedule) const{

  if(index_schedule<_schedule.size()){
    auto it=std::next(_schedule.cbegin(), index_schedule);
    double beta=it->first;
    return 1.0/beta;
    
  }
  else{
  
    std::cerr << "!! Index out of range !!\nRequired temperature not present in the Annealing Schedule" << std::endl;
    return -1;
  
  }
  
}


void annealer::print_path() const{

  _p.print_chromosome();

}


void annealer::print_schedule() const{

  unsigned int j=1;
  std::cout << "--------------------" << std::endl;
  std::cout << "Annealing Schedule" << std::endl;
  std::cout << "--------------------" << std::endl;
  
  if(_schedule.size()<=10){
  
    //for(auto it=_schedule.cbegin(); it!=_schedule.cend(); ++it){  //1st way
    for(const auto& elem : _schedule){  //2nd way : C++ for in loop
  
      std::cout << "\t(β_" << j << ", n_" << j << ") = (" << elem.first << ", " << elem.second << ")\n";
      j++;
  
    }
    
  }
  else{
  
    auto it1=_schedule.cbegin();
    auto it2=std::next(it1, _schedule.size()-2);
    for(int k=0; k<4; k++){
    
      std::cout << "\t(β_" << k+1 << ", n_" << k+1 << ") = (" << it1->first << ", " << it1->second << ")\n";
      it1++;
    
    }
    for(int w=0; w<4; w++)
      std::cout << "\t..." << std::endl;
    for(int z=0; z<2; z++){
    
      std::cout << "\t(β_" << _schedule.size()-1+z << ", n_" << _schedule.size()-1+z << ") = (" << it2->first << ", " << it2->second << ")\n";
      it2++;
    
    }
  
  
  }

}


void annealer::check() const{

  if(_p.check())
    std::cout <<"The path to be optimized fulfils all the bonds imposed by the TSP" << std::endl;
  else
    std::cout <<"The path to be optimized does not fulfill all the bonds imposed by the TSP" << std::endl;

}


void annealer::pair_permutation(){

  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd(1, _p.get_length()-1); 
  unsigned int g1, g2;
  
    
  g1=rnd(gen);  //choose the two 
  g2=rnd(gen);  //genes to swap random
    
  while(g1==g2)
    g2=rnd(gen);
    
  _p.swap_gene(g1, g2);  //this already recalculates the fitness 
  			 
}


void annealer::pair_permutation(chromosome& target){

  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd(1, target.get_length()-1); 
  unsigned int g1, g2;
  
    
  g1=rnd(gen);  //choose the two 
  g2=rnd(gen);  //genes to swap random
    
  while(g1==g2)
    g2=rnd(gen);
    
  target.swap_gene(g1, g2);  //this already recalculates the fitness 
  			 
}


void annealer::contiguous_permutation(){

  std::random_device rd;  
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd1(1, _p.get_length()/2);  //the contiguous cities in the chosen chromosome
  std::uniform_int_distribution<> rnd2(1, _p.get_length()-1);				       
  unsigned int start_gene1=0, start_gene2=0, m=0;
  
  
  start_gene1=rnd2(gen);  //the 1st-first of the 1st m-contiguous cities
  start_gene2=rnd2(gen);  //the 2nd-first of the 2nd m-contiguous cities
  m=rnd1(gen);  //the number of contiguous cities
    
  while(start_gene1+m >= _p.get_length() || start_gene2+m >= _p.get_length() || start_gene2 <= start_gene1+m){
      
    m=rnd1(gen);
    start_gene1=rnd2(gen);
    start_gene2=rnd2(gen);
    
  }
      
  for(int count=0; count<m; count++)
    _p.swap_gene(start_gene1+count, start_gene2+count);

}


void annealer::contiguous_permutation(chromosome& target){

  std::random_device rd;  
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd1(1, target.get_length()/2);  //the contiguous cities in the chosen chromosome
  std::uniform_int_distribution<> rnd2(1, target.get_length()-1);				       
  unsigned int start_gene1=0, start_gene2=0, m=0;
  
  
  start_gene1=rnd2(gen);  //the 1st-first of the 1st m-contiguous cities
  start_gene2=rnd2(gen);  //the 2nd-first of the 2nd m-contiguous cities
  m=rnd1(gen);  //the number of contiguous cities
    
  while(start_gene1+m >= target.get_length() || start_gene2+m >= target.get_length() || start_gene2 <= start_gene1+m){
      
    m=rnd1(gen);
    start_gene1=rnd2(gen);
    start_gene2=rnd2(gen);
    
  }
      
  for(int count=0; count<m; count++)
    target.swap_gene(start_gene1+count, start_gene2+count);

}


void annealer::inversion(){

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd(1, _p.get_length()-1);  //the contiguous cities in the chosen chromosome					       
  unsigned int start_gene=0, m=0, start=0, end=0;
  
  
  start_gene=rnd(gen);  //the first of the m-contiguous cities
  m=rnd(gen);  //the number of contiguous cities
    
  while(start_gene+m >= _p.get_length()){
    
    m=rnd(gen);
    start_gene=rnd(gen);
      
  }
      		          
  start=start_gene;
  end=start_gene+m-1;
    
  while(start<end){
    
    _p.swap_gene(start, end);
    start++;
    end--;
    
  }
  
}


void annealer::inversion(chromosome& target){

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rnd(1, target.get_length()-1);  //the contiguous cities in the chosen chromosome					       
  unsigned int start_gene=0, m=0, start=0, end=0;
  
  
  start_gene=rnd(gen);  //the first of the m-contiguous cities
  m=rnd(gen);  //the number of contiguous cities
    
  while(start_gene+m >= target.get_length()){
    
    m=rnd(gen);
    start_gene=rnd(gen);
      
  }
      		          
  start=start_gene;
  end=start_gene+m-1;
    
  while(start<end){
    
    target.swap_gene(start, end);
    start++;
    end--;
    
  }
  
}


void annealer::annealing(unsigned int index_schedule, double prob_PP, double prob_CP, double prob_IN){

  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> rnd(0, 1);
  chromosome* new_path;
  double energy=0, beta=0, n=0;
  unsigned int accept=0;
  
  //Define the current temperature
  if(index_schedule<_schedule.size()){
    auto it=std::next(_schedule.cbegin(), index_schedule);
    beta=it->first;  //actual inverse temperature and
    n=it->second;    //Metropolis steps in the annealing
  }
  else
    std::cerr << "!! Index out of range !!\nRequired temperature not present in the Annealing Schedule" << std::endl;
  
  //Sample new path configuration via 
  //Metropolis Algorithm
  //std::cout << std::fixed << std::showpoint << std::setprecision(2) << "T = " << 1.0/beta << " \t";
  for(int M=0; M<n; M++){
  
    //if(M%50==0)
      //std::cout << "*"; 
    new_path=new chromosome(_p);
    /*
      Propose a new path 
      by applying genetic mutation
    */
    if(rnd(gen)<prob_PP)
      pair_permutation(*new_path);
    if(rnd(gen)<prob_CP)
      contiguous_permutation(*new_path);
    if(rnd(gen)<prob_IN)
      inversion(*new_path);
      
    /*
      Accept or reject the new path
      using Metropolis acceptance technique
    */
    //H=L
    energy=new_path->get_L()-_p.get_L();
    if(rnd(gen) < std::min(1., exp(-beta*energy))){
    
        set_path(*new_path);
        accept++;
    
    }
      
  }
  
  set_acceptance(((1.0*accept)/(1.0*n))*100.0);
  //std::cout << "\t acceptance ratio = " << std::fixed << std::showpoint << std::setprecision(1) << get_acceptance() << "\%" << std::endl;
  delete new_path;
  
  /*
  std::cout << "prob_PP = " << (PP*1.0)/n << std::endl;
  std::cout << "prob_CP = " << (CP*1.0)/n << std::endl;
  std::cout << "prob_IN = " << (IN*1.0)/n << std::endl;
  */

}
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
