#ifndef MATRIXPROTEIN_HPP
#define MATRIXPROTEIN_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <cmath>
#include <map>
#include <cassert>
#include <stdexcept>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "utilities.hpp"

using namespace std;


typedef vector<vector<double> > matrix;

 /*!
 *The Pattern struct groups together all the informations related to one pattern.  Pattern contains motif itself (site),  
 *the bScore, the chromosome number, the position of this motif, and the direction (forward '+' or reverse '-'). 
 * */
 struct Pattern
{
    double bScore;
    string site;
    int pos;
    string chrNb;
    char dir;
    
};

struct SeqPos {
    vector<char> sequence;
    size_t position;
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
    
 /*!
  * Constructor.
  * */
	MatrixProtein();
	
 /*!
  * Destructor.
  * */
	~MatrixProtein();
    
    //----- Methods -----
    
    // Patterns (i.e Motifs)
 /*!
 * Methods that takes a list of sites (even size) and set the Patterns attribute with the list of 
 * all the site having an affinity score above a certain threshold (can be given or 
 * taken by default(the default number isn't 0 this is just a way to assure that if the 
 * user gives 0 or nothing the same default threshold will be used (see calculation of the default threshold)))
 * */
    void fillVectorPatterns(vector<vector<char> > sites, double threshold = 0 );
    
 /*!
 * Methods takes a list of motif (all the motif don't need to have the same lenght) and returns the 
 * list of all the site of size 'set' having an affinity score above a certain threshold (this threshold 
 * can either be given (last argument) or one will be calculated by default. (see calculation of default threshold))
 * */
	void get_relevent_site(vector<vector<char> > Input, int set, int = 0);
	
 /*!
 * The findPatterns method finds all theoretical motifs from a .mat file. With a given threshold, it calculates all(depending on the size of the .mat, yet assumed seven columns long) motifs and keeps those with a score above the threshold.
 * */
    vector<Pattern> findPatterns(string matFile, double threshold); 
    
 /*!
 * Return n/tot value.
 * */
    double probas(double n, double tot);
    
 /*!
 * Loadmatrix_fromfile method allows us to recover data from a file and stores it in the mx Matrix.
 * */
    void loadmatrix_fromfile(const string& Data);
    
 /*!
 * Swap a PWM rel/abs matrix to its PSSM form.
 * */
    void swaptopssm(matrix mtx);
    
 /*!
 * Swap a PSSM rel/abs matrix to its PWM form.
 * */
    void swaptopwm(matrix mtx);
    
 /*!
 * Swap a PSSM/PWM relative matrix to its absolute form.
 * */
    void swaptoabsolute(matrix& mtx);
    
 /*!
 * Swap a PSSM/PWM absolute matrix to its relative form.
 * */
    void swaptorelative(matrix& mtx);
        
 /*!
 * A bool equal to 1 if the matrix is absolute, 0 otherwise. This function checks every line of the matrix.
 * */
    bool absolute(matrix matrice);
 /*!
 * Utilitary method used in possible to compute in matrix lines and columns.
 * */
    void PWM_to_PSSM(matrix& matrice);
 /*!
 * Utilitary method used in possible to compute in matrix lines and columns.
 * */
    void PWM_to_PSSM_2(matrix& matrice);
    
 /*!
  * Utilitary method used to check if we multiplied by 0.25 or not matrix slots.
  * */
    bool which_PWM_to_PSSM(matrix matrice);
 
 /*!
  * Utilitary method that returns 1 if the matrix is a PWM and checks if the matrix is in a proper format.
  * */
    bool check_if_pmworpssm(matrix matrice);
    
 /*!
  * Method that checks if a matrix is possible, testing twice.
  * */
    bool possible(matrix matrice);
    
 /*!
  * Method that return the matrix status (PWM/PSSM & abs/rel) in a vec <bool>.
  * */
    vector <bool> matrix_status(matrix matrice);
    
 /*!
  * Method that allows us to load all matrices from any single one that has the proper format.
  * */
    void matrix_generation();
    
 /*!
  * Method that defines a minimum value for matrices slots.
  * */
    void readjust_values(matrix& mtx);
    
 /*!
  * Method that calculates the binding score (double) of a certain sequence based on a PMW matrix.
  * */
    double get_affinity_score_from_matrix(matrix mx,vector<char> sequence);
    
 /*!
  * Method that displays the PWM rel.
  * */
    void display_PWM_rel();
    
 /*!
  * Getters.
  * */
    matrix getpssm_abs();
    matrix getpssm(matrix);
    matrix getpssm_rel();
    matrix getpwm_rel();
    matrix getpwm_abs();
    matrix getmx();
	vector<Pattern> getPatterns() const {return patterns;};
 /*!
  * Setters.
  * */
    void setpssm_rel(matrix mtx);
    void setpwm_rel(matrix mtx);
    void set_mx(matrix mtx);
    void setPatterns(vector<Pattern>);
   
 /*!
 * Methods that computes the score of a sequence.
 * */
    void calcul( vector<char> tab_, vector<double>& tabA_, vector<double>& tabT_, vector<double>& tabG_, vector<double>& tabC_);
    
    // ----- Fonctions utilisées pour l'algorithme EM -----
    
 /*! 
 * The "RemetZero" method allows to reset a vector by putting each box to Zero without changing the lenght of the vector.  
 * */
	void RemetZero( vector<double>& vec_);
  
 /*!
 * The "calculScore" method calculate the score of a motif tab2_ with the matrix finale and put the score in 
 * the box j_ of the table  score_. 
 * */
	void calculScore(vector<vector<double> > finale_, vector<double>& score_, int j_, vector<char> tab2_);
    
 /*!
 * The method "FindMotif" puts the motifs with the best scores in the vector best_seqs, using the matrix finale
 * and the vector of string FromFasta_. Each motif in best_seqs is represnted by a SeqPos, which contains also 
 * its position. 
 * */
    void FindMotif(vector<vector<double> > finale_,vector<string> FromFasta_, int longueur_motif_, vector<vector<SeqPos> >& best_seqs_);
    
 /*!
 * The "fillPattern" method fill the attribute "patterns" with the motifs with the higher scores coming from "best_seqs_". 
 * In this pattern, the score, the position on the sequence, the direction, the chromosome number are specified. 
 * This fills the Pattern with the relevant information : The position given for a reverse seq is the position of 
 * it's foward seq.
 * */
	void fillPattern(vector<vector<SeqPos> > best_seqs_, int sizeint, vector<int>);
    
  /*!
 * The "calculScoreFinal" method calculates the score of a given sequence "seq". 
 * */
	double calculScoreFinal(string seq);
    
 /*!
 * The "enleveZero" method removes all the zeros of the matrix mx and add a very small number to each box of the matrix, 
 * then divide again to make sure that the sum of line is still 1. 
 * */
	void enleveZero(vector<vector<double> >& mx, double somme_);
    
 /*!
 * The "EMalgorithm" method creates the PWM matrix from a Fasta file, and put the result in the attribute
 * mx. This method also fills the attribute patterns with the the motifs (length = longueur_motif) with the higher scores. 
 * */
    // Algorithme EM complet --> Met la matrice PWM dans l'attribut mx + change le pattern
    void EMalgorithm(int longueur_motif, vector<string> FromFasta, vector<int>, int);
    
    //Calculation of the default threshold
 /*!
  * Generates a random nucleotide.
  * */
	char genRandomChar();
	
 /*!
  * Generates a random sequence.
  * */
    vector<char> seq_generator(double length);
    
 /*!
  * Set average score for defined length sequences.
  * */
    double set_average(matrix Matrice, double size);

};

#endif
