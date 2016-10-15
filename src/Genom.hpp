
#ifndef GENOM_HPP
#define GENOM_HPP

#include "Protein.hpp"

#include <stdio.h>
#include <vector>

using namespace std;

class Genom
{
public:
    Genom();
    ~Genom();
    
private:
    //Sequence sequence;
    vector<Protein> proteins;
    
};




#endif /* Genom_hpp */
