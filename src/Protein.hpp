#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <vector>

using namespace std;

typedef std::vector<std::vector<double> > matrix;

struct Pattern
{
    double bScore;
    vector<char> site;
};



class Protein
{
<<<<<<< Updated upstream
private:
    //Attributs
    vector<Pattern> patterns;
    matrix pssm_abs;
    matrix pwm_abs;
    matrix pssm_rel;
    matrix pwm_rel;
    
    
public:
    
    
    //Constructeur et destructeur
    Protein();  						//voir ce qu'on veut initialiser et à quelles valeurs?
    ~Protein();
    
    //Méthodes
    
    //Patterns
    void fillPattern(double const& bScore, vector<char> const& site);   		//à utiliser éventuellement dans le constructeur? const ou pas?
    void fillVectorPatterns(vector<vector<char>> sites, double threshold = 0 ); // fills Patterns from a list of binding sites, if the user doesn't give a threshold it is calculated as the avarage binding score from random sequence (NB this threshold is also used if the user enter 0 since this only filter the random binding)
    
    double probas(double n, double tot);
    
    
    //Matrix
    double get_afinity_score_from_matrix(vector<vector<double>> matrix,vector<char> sequence);
    void loadmatrix_fromscore();
    void display_PWM(matrix finale);
    void setpwm_abs(matrix Matrix);
    
    //Load matrix from a file
    matrix loadmatrix_fromfile(std::string Data);
    vector<Pattern> getPatterns();
    
    //Calculation of default threshold
    vector<char> seq_generator(double length);
    char genRandomChar();
    double set_average(matrix Matrice, double size);
    
=======
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
        vector<Pattern> getPatterns();
        void check_if_empty(string Data);
        
	
>>>>>>> Stashed changes
};


#endif

