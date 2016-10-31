#include "Protein.hpp"

using namespace std;

// ==============================================================================================================CONSTRUCTEUR ET DESTRUCTEUR

Protein::Protein() 
{
	
}

Protein::~Protein()
{
	
}

// ==============================================================================================================METHODES

void Protein::fillPattern(double const& bScore, vector<char> const& site)
{
	pattern.bScore = bScore;
	pattern.site = site;
}


void Protein::fillVectorPatterns(Pattern pattern)					//à revoir en fonction de si on pré remplit le tableau ou pas
{
    /*if (matrix.empty()) {
		matrix.push_back(pattern);
	} else {
		for(size_t i(0); i < matrix.size() ; ++i) {
			if (matrix[i].empty()) { //comment vérifier empty() pour une seule case?
				matrix[i] = pattern;
			} else {
				++i;
			}
		}
    }	*/
}


double probas(double n, double tot)
{
    return n/tot;
}



Matrice create_PWM()
{
    // typedef std::vector< Motif>Tableau;
    typedef std::vector<std::vector<double>> Matrice;
    // Tableau sequences;
    double compteurA(0);
    double compteurT(0);
    double compteurG(0);
    double compteurC(0);
    double nombreTot(0);
    std::vector <double> tabA;
    std::vector <double> tabT;
    std::vector <double> tabG;
    std::vector <double> tabC;
    Matrice finale;
    
    /* Motif m1 = {{'A', 'T', 'G', 'C', 'T', 'G', 'T'}, 0.8};
     sequences.push_back( m1) ;
     Motif m2 = {{'T', 'C', 'A', 'C', 'T','C','C'}, 0.2};
     sequences.push_back( m2) ;
     Motif m3 = {{'T', 'T', 'G', 'C' , 'C','A', 'C'}, 0.4};
     sequences.push_back( m3) ;
     Motif m4 = {{ 'A' , 'A', 'C', 'G', 'G','T', 'C'}, 0.7};
     sequences.push_back( m4) ;
     */
    
    
    for (size_t j(0) ; j < 7 ; ++j)
    {
        
        for ( size_t i(0) ; i < patterns.size() ; ++i)
        {
            
            nombreTot = nombreTot + patterns[i].Bscore ;
            if ( sequences[i].site[j] == 'A')
            {
                compteurA = compteurA + patterns[i].Bscore ;
            }
            if ( sequences[i].site[j] == 'T')
            {
                compteurT = compteurT + patterns[i].Bscore ;
            }
            if ( sequences[i].site[j] == 'G')
            {
                compteurG = compteurG + patterns[i].Bscore ;
            }
            if ( sequences[i].site[j] == 'C')
            {
                compteurC = compteurC + patterns[i].Bscore ;
            }
        }
        
        tabA.push_back(probas(compteurA , nombreTot));
        tabT.push_back(probas(compteurT , nombreTot));
        tabG.push_back(probas(compteurG , nombreTot));
        tabC.push_back(probas(compteurC , nombreTot));
        compteurA = 0;
        compteurT = 0;
        compteurG = 0;
        compteurC = 0;
        nombreTot = 0;
    }
    
    finale.push_back(tabA);
    finale.push_back(tabT);
    finale.push_back(tabG);
    finale.push_back(tabC);
    
    return finale;
    
}


void display_PWM( Matrice finale)
{
    
    for ( int i(0); i < finale.size() ; ++i)
    {
        
        for (int j(0) ; j < finale[i].size() ; ++j)
        {
            
            std::cout  << std::setprecision(5) << std::setw(5) << finale[i][j] << " | " ; 
        } 
        
        std::cout << std::endl; 
    }
    
}