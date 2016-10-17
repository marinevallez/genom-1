
#ifndef GENOME_HPP
#define GENOME_HPP

#include "Protein.hpp"

#include <stdio.h>
#include <vector>

using namespace std;

class Genome
{
public:
    Genome();
    ~Genome();
    
private:
    //Sequence sequence;
    vector<Protein> proteins;
    
};




#endif /* Genome_hpp */
