#ifndef PROTEIN_HPP
#define PROTEIN_HPP

//#include "Matrix.hpp"

#include <vector>

using namespace std;

struct Pattern {
			double bScore;
			vector<char> site;
		}

class Protein 
{
	private:
		//Attributs
		Matrix matrix;
		Pattern pattern;
		vector<Pattern> patterns;
	
	public:		
		//Constructeur et destructeur
		Protein();  						//voir ce qu'on veut initialiser et à quelles valeurs?
		~Protein();
	
		//Méthodes
		fillPattern(double const& bScore, vector<char> const& site);   		//à utiliser éventuellement dans le constructeur? const ou pas?
		fillVectorPatterns(Pattern pattern);								//const ou pas? Comment faire pour mettre directement plusieurs patterns?
	
}


#endif
