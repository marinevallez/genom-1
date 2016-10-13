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

