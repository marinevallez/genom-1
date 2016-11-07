
#include "Letters.hpp"
#include <math.h>




void test(double& L) // assure that the size of a lettre isn't equal to 2 if the probability is 0 (for aesthtic consideration)
{
    if (L == 0) {
        L = pow(10, -100);
    }
}



Letters::Letters(vector<Motif> motif_)
{
    double compteurA(0);
    double tailleA(0);
    
    double compteurT(0);
    double tailleT(0);
    
    double compteurG(0);
    double tailleG(0);
    
    double compteurC(0);
    double tailleC(0);
    
    double nombreTot(0);
    for (size_t j(0) ; j < motif_[0].lettres.size() ; ++j)
    {
        for ( size_t i(0) ; i < motif_.size() ; ++i)
        {
            
            nombreTot = nombreTot + motif_[i].binding_score ;
            if ( motif_[i].lettres[j] == 'A')
            {
                compteurA = compteurA + motif_[i].binding_score*log2(motif_[i].binding_score) ; // calculates the Shannon entropy Hi=-∑ƒA,i*log2(ƒA,i) with ƒa,i the relative frequency of base A at position i
                tailleA += (2 - compteurA)*motif_[i].binding_score; // calculates information content given by log2(4)-Hi, the 4 is given by the 4 base possible and gets the height (Ri*ƒa,i)
                
            }
            if ( motif_[i].lettres[j] == 'T')
            {
                compteurT = compteurT + motif_[i].binding_score*log2(motif_[i].binding_score) ; // calculates the Shannon entropy Hi=-∑ƒT,i*log2(ƒT,i) with ƒT,i the relative frequency of base T at position i
                tailleT += (2 - compteurT)*motif_[i].binding_score;// calculates information content given by log2(4)-Hi, the 4 is given by the 4 base possible and gets the height (Ri*ƒa,i)
                
            }
            if ( motif_[i].lettres[j] == 'G')
            {
                compteurG = compteurG + motif_[i].binding_score*log2(motif_[i].binding_score) ;// calculates the Shannon entropy Hi=-∑ƒA,i*log2(ƒG,i) with ƒG,i the relative frequency of base G at position i
                tailleG += (2 - compteurG)*motif_[i].binding_score;// calculates information content given by log2(4)-Hi, the 4 is given by the 4 base possible and gets the height (Ri*ƒa,i)
                
            }
            if ( motif_[i].lettres[j] == 'C')
            {
                compteurC = compteurC + motif_[i].binding_score*log2(motif_[i].binding_score) ; // calculates the Shannon entropy Hi=-∑ƒC,i*log2(ƒC,i) with ƒC,i the relative frequency of base C at position i
                tailleC += (2 - compteurC)*motif_[i].binding_score;// calculates information content given by log2(4)-Hi, the 4 is given by the 4 base possible, and gets the height (Ri*ƒa,i)
                
            }
            
        }
        test(tailleA); // verifes that the size isn't 0, if it is, makes very close to 0.
        test(tailleT);
        test(tailleG);
        test(tailleC);
        
        Letters_size.push_back({tailleA ,tailleC, tailleG ,tailleT});
        compteurA = 0;
        tailleA = 0;
        compteurT = 0;
        tailleT = 0;
        compteurG = 0;
        tailleG = 0;
        compteurC = 0;
        tailleC = 0;
        nombreTot = 0;
        
    }
    
}

vector<vector<double>>  Letters::get_lettre_size() // getter for the size of each lettre
{
    return Letters_size;
}


