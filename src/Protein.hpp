#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <vector>
#include "matrix.hpp"
using namespace std;

struct Pattern {
			double bScore;
			vector<char> lisOfSites;
        };



class Protein 
{
	private:
		//Attributs
		Pattern pattern;
		vector<Pattern> patterns;
		Matrix mtrx;
	
	public:		
		//Constructeur et destructeur
		Protein();  						//voir ce qu'on veut initialiser et à quelles valeurs?
		~Protein();
	
		//Méthodes
        void fillPattern(double const& bScore, vector<char> const& site);   		//à utiliser éventuellement dans le constructeur? const ou pas?
        void fillVectorPatterns(Pattern pattern);								//const ou pas? Comment faire pour mettre directement plusieurs patterns?
    
        double probas(double n, double tot);
		matrix loadmatrix_fromscore();
        void display_PWM(matrix finale);
        //Load matrix from a file
        matrix loadmatrix_fromfile(std::string Data);
        
	
};


#endif
