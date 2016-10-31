#ifndef PROTEIN_HPP
#define PROTEIN_HPP

//#include "Matrix.hpp"

#include <vector>

using namespace std;

struct Pattern {
			double bScore;
			vector<char> lisOfSites;
        };



class Protein 
{
	private:
		//Attributs
        //Matrix matrix;
		Pattern pattern;
		vector<Pattern> patterns;
	
	public:		
		//Constructeur et destructeur
		Protein();  						//voir ce qu'on veut initialiser et à quelles valeurs?
		~Protein();
	
		//Méthodes
        void fillPattern(double const& bScore, vector<char> const& site);   		//à utiliser éventuellement dans le constructeur? const ou pas?
        void fillVectorPatterns(Pattern pattern);								//const ou pas? Comment faire pour mettre directement plusieurs patterns?
    
        double probas(double n, double tot);
        Matrice create_PWM();
        void display_PWM(Matrice finale);
	
};


#endif
