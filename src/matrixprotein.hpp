#ifndef MATRIXPROTEIN_HPP
#define MATRIXPROTEIN_HPP
#include <iostream>
#include <vector>
using namespace std;
typedef vector<vector<double> > matrix;


struct Pattern
{
    double bScore;
    vector<char> site;
};



class MatrixProtein
{
private:
    
    //-----Attributs-----
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
    // Patterns
    void fillVectorPatterns(vector<vector<char> > sites, double threshold );
    vector<Pattern> getPatterns() const;
    void setPatterns(vector<Pattern>);
    
    //general
    double probas(double n, double tot);
    double To_double (const string& string );
    
    //Matrix
    //void display_PWM(matrix finale);
    void loadmatrix_fromscore();
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
    
    //calculation of the default threshold
    char genRandomChar();
    vector<char> seq_generator(double length);
    double set_average(matrix Matrice, double size);
    
};


#endif
