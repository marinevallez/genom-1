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
 * 
 * The Pattern struct groups together all the informations related to one pattern.  Pattern contains motif itself (site),  
 * the bScore, the chromosome number, the position of this motif, and the direction (forward '+' or reverse '-'). 
 * 
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
    
    //-----Constructor & Destructor-----
    MatrixProtein();
    ~MatrixProtein();
    
    //-----Methods-----
    // Patterns (i.e Motifs)
    void fillVectorPatterns(vector<vector<char> > sites, double threshold = 0 ); // takes a list of sites (even size) and set the Patterns attribute with the list of all the site having an affinity score above a certain threshold (can be given or taken by default(the default number isn't 0 this is just a way to assure that if the User gives 0 or nothing the same default threshold will be used (see calculation of the default threshold)))
    
    vector<Pattern> getPatterns() const {return patterns;};
    
    void setPatterns(vector<Pattern>);
    
    void get_relevent_site(vector<vector<char> > Input, int set, int = 0);//this function takes a list of motif (all the motif don't need to have the same lenght) and return the list of all the site of size 'set' having a affinity score above a certain threshold (this threshold can either be given (last argument) or one will be calculated by default (see calculation of default threshold))
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
    void swaptopssm(matrix mtx);
    void swaptopwm(matrix mtx);
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
    double get_affinity_score_from_matrix(matrix mx,vector<char> sequence);
    void display_PWM_rel();
    matrix getpssm_abs();
    matrix getpssm(matrix);
    matrix getpssm_rel();
    matrix getpwm_rel();
    matrix getpwm_abs();
    matrix getmx();
    void setpssm_rel(matrix mtx);
    void setpwm_rel(matrix mtx);
    void set_mx(matrix mtx);
    void calcul( vector<char> tab_, vector<double>& tabA_, vector<double>& tabT_, vector<double>& tabG_, vector<double>& tabC_);
    
    // fonctions utilisée pour l'algorithme EM;
    
 /*!
 * 
 * The "RemetZero" method allows to reset a vector by putting each box to Zero without changing the lenght of the vector. 
 * 
* */

    void RemetZero( vector<double>& vec_);
  
/*!
 * The "calculScore" method calculate the score of a motif tab2_ with the matrix finale and put the score in 
 * the box j_ of the table  score_. 
 * 
 * */
 
    void calculScore(vector<vector<double> > finale_, vector<double>& score_, int j_, vector<char> tab2_);
    
/*!
 * 
 *  The method "FindMotif" puts the motifs with the best scores in the vector best_seqs, using the matrix finale
 *  and the vector of string FromFasta_. Each motif in best_seqs is represnted by a SeqPos, which contains also 
 * its position. 
 * 
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
 * 
 * */
 
    void enleveZero(vector<vector<double> >& mx, double somme_);
    
 /*!
 * The "EMalgorithm" method creates the PWM matrix from a Fasta file, and put the result in the attribute
 * mx. This method also fills the attribute patterns with the the motifs (length = longueur_motif) with the higher scores. 
 * 
 * */
    // Algorithme EM complet --> Met la matrice PWM dans l'attribut mx + change le pattern
    void EMalgorithm(int longueur_motif, vector<string> FromFasta, vector<int>, int);
    
    //Calculation of the default threshold
    char genRandomChar();
    vector<char> seq_generator(double length);
    double set_average(matrix Matrice, double size);

};

#endif
