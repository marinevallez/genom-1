#ifndef MATRIXPROTEIN_HPP
#define MATRIXPROTEIN_HPP

#include <iostream>
#include <vector>
using namespace std;


typedef vector<vector<double> > matrix;

struct Pattern
{
    double bScore;
    string site;
};


class MatrixProtein
{
	
	private:
    
    //-----Attributes-----
    vector<Pattern> patterns;
    //Temporary matrix with which we generate all other types of matrix
    matrix mx;
    matrix pssm_abs;
    matrix pwm_abs;
    matrix pssm_rel;  
    matrix pwm_rel;
    
	public:
    
    //-----Constructor & Destructor-----
    MatrixProtein();
    ~MatrixProtein();
    
    //-----Methods-----
    
   /* //Conversions 
    double To_double(const string& string);
    int toInt(const string& str); 
	vector<char> toVector(string str);
	string toString(vector<char> vec);*/


  
    // Patterns (i.e Motifs)
    void fillVectorPatterns(vector<vector<char> > sites, double threshold = 0 ); // takes a list of sites (even size) and set the Patterns attribute with the list of all the site having an affinity score above a certain threshold (can be given or taken by default(the default number isn't 0 this is just a way to assure that if the User gives 0 or nothing the same default threshold will be used (see calculation of the default threshold)))
    vector<Pattern> getPatterns() const;
    void setPatterns(vector<Pattern>);
    void get_relevent_site(vector<vector<char>> Input, int set, int = 0);//this function takes a list of motif (all the motif don't need to have the same lenght) and return the list of all the site of size 'set' having a affinity score above a certain threshold (this threshold can either be given (last argument) or one will be calculated by default (see calculation of default threshold))
    /*!
     * The findPatterns method finds all theoretical motifs from a .mat file. With a given threshold, it calculates all(depending on the size of the .mat, yet assumed seven columns long) motifs and keeps those with a score above the threshold.
     */
    vector<Pattern> findPatterns(string matFile, double threshold); //finds motifs from a matrix
    
    //General
    double probas(double n, double tot);
    
    //Matrix
    //void display_PWM(matrix finale);
    void loadmatrix_fromfile(const string& Data);
    void setrw(int value);
    void swaptopssm(matrix& mtx);
    void swaptopwm(matrix& mtx);
    void swaptoabsolute(matrix& mtx);
    void swaptorelative(matrix& mtx);
    bool absolute(matrix matrice);
    void PWM_to_PSSM(matrix& matrice);
    void PWM_to_PSSM_2(matrix& matrice);
    bool which_PWM_to_PSSM(matrix matrice);
    bool check_if_pmworpssm(matrix matrice);
    bool possible(matrix matrice);
    vector <bool> matrix_status(matrix matrice);
    void matrix_generation();
    void readjust_values(matrix& mtx);
    double get_affinity_score_from_matrix(vector<vector<double> > matrix,vector<char> sequence);
    void display_PWM_rel();
    matrix getpwm_abs();
    matrix getpssm_abs();
    matrix getmx()
    void calcul( vector<char> tab_, vector<double>& tabA_, vector<double>& tabT_, vector<double>& tabG_, vector<double>& tabC_);

	// fonctions utilisée pour l'algorithme EM; 
	void RemetZero( vector<double>& vec_);
	void calculScore(vector<vector<double>> finale_, vector<double>& score_, int j_, vector<char> tab2_); 
	void FindMotif(vector<vector<double>> finale_,vector<string> FromFasta_, int longueur_motif_, vector<vector<vector<char>>>& best_seqs_);
	void fillPattern(vector<vector<vector<char>>> best_seqs_);
	double calculScoreFinal(string seq);
	
	
	// Algorithme EM complet --> Met la matrice PWM dans l'attribut mx + change le pattern
	void EMalgorithm(int longueur_motif, vector<string> FromFasta);  
    
    //Calculation of the default threshold
    char genRandomChar();
    vector<char> seq_generator(double length);
    double set_average(matrix Matrice, double size);
    
};

#endif
