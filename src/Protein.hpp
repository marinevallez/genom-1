#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <vector>

using namespace std;

struct Pattern 
{
		double bScore;
		vector<char> listOfSites;
};



class Protein 
{
	private:
		//Attributs
		vector<Pattern> patterns;
		Matrix mtrx;
	
	public:		
	
	
		//Constructeur et destructeur
		Protein();  						//voir ce qu'on veut initialiser et à quelles valeurs?
		~Protein();
	
		//Méthodes
        void fillPattern(double const& bScore, vector<char> const& site);   		//à utiliser éventuellement dans le constructeur? const ou pas?
        void fillVectorPatterns(vector<vector<char>> sites, double threshold = 0 ) // fills Patterns from a list of binding sites, if the user doesn't give a threshold it is calculated as the avarage binding score from random sequence (NB this threshold is also used if the user enter 0 since this only filter the random binding)
    
        double probas(double n, double tot);
		matrix loadmatrix_fromscore();
        void display_PWM(matrix finale);
        //Load matrix from a file
        matrix loadmatrix_fromfile(std::string Data);
        vector<Pattern> getPatterns();
        
	
};


#endif
