#ifndef MATRIXPROTEIN_HPP
#define MATRIXPROTEIN_HPP
#include <iostream>
#include <vector>
using namespace std;
typedef vector<vector<double> > matrix;


struct Pattern 
{
		double bScore;
		vector<char> listOfSites;
};



class MatrixProtein 
{
	private:
	
		//-----Attributs-----
		Pattern pattern;
		vector<Pattern> patterns;
		//Matrice provisoire grace à laquelle on génère toutes les autres.
		matrix mx;
		matrix pssm_abs;
		matrix pwm_abs;
		matrix pssm_rel; 
		matrix pwm_rel;
	
	public:		
	
	
		//-----Constructeur&Destructeur------
		MatrixProtein();  						
		~MatrixProtein();
	
		//-----Méthodes-----
        void fillPattern(double const& bScore, vector<char> const& site);   		//à utiliser éventuellement dans le constructeur? const ou pas?
        void fillVectorPatterns(Pattern pattern);								//const ou pas? Comment faire pour mettre directement plusieurs patterns?
    
        double probas(double n, double tot);
        //void display_PWM(matrix finale);
        void loadmatrix_fromscore();
        void loadmatrix_fromfile(string Data);
        vector<Pattern> getPatterns();
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
		double To_double (const string& string );
        
	
};


#endif
