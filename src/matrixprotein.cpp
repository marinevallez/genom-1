#include "matrixprotein.hpp"
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
#include "utilities.cpp"

using namespace std;

/*//-----Converions-----

double MatrixProtein::To_double( const string& string ) // allows a convertion from string to double
{
    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}


int MatrixProtein::toInt(const string& str) // allows a convertion from string to int
{
    istringstream stream(str);
    int a;
    if (!(stream >> a))
        return 0;
    return a;
}

vector<char> MatrixProtein::toVector(string str) 
{
	vector<char> vec;
	for (size_t i(0); i <= str.size(); ++i)
	{
		vec.push_back(str[i]);
	}
	return vec;
}

string MatrixProtein::toString(vector<char> vec)
{
	string str;
	for (size_t i(0); i <= vec.size(); ++i)
	{
		str += vec[i];
	}
	return str;
}*/
// ==============================================================================================================CONSTRUCTEUR ET DESTRUCTEUR

MatrixProtein::MatrixProtein()
{
    /*Il nous faut une méthode qui choisi loadmatrix_fromscore ou loadmatrix_fromfile
     * pour génerer les différentes matrices*/
}

MatrixProtein::~MatrixProtein()
{
    
}

// ==============================================================================================================METHODES


void MatrixProtein::fillVectorPatterns(vector<vector<char> > sites, double threshold )// Should we ask if the user rather have a binding score calculated from a relative matrix ? woukd that make a difference (NB this function calculates the BS from a PWM absolute) this takes all the binding site for a given prot(in the fasta file)

{
    double binding_score(0);
    if (threshold == 0) {
        threshold = set_average(pwm_abs, sites[0].size());
    }
    
    for (size_t i(0); i < sites.size(); ++i) {
        
        try {
            binding_score = get_affinity_score_from_matrix(pwm_abs, sites[i]);
        } catch (const runtime_error& error) {
            cout << error.what();
        }
        
        
        if(binding_score > threshold)
        {
            Pattern pattern;
            pattern.bScore = binding_score; 
            pattern.site = toString(sites[i]);
            patterns.push_back(pattern);
            
        }
        
        
    }
}




double MatrixProtein::probas(double n, double tot)
{
    return n/tot;
}


/*
void MatrixProtein::loadmatrix_fromscore()
{
    
    double compteurA(0);
    double compteurT(0);
    double compteurG(0);
    double compteurC(0);
    double nombreTot(0);
    std::vector <double> tabA;
    std::vector <double> tabT;
    std::vector <double> tabG;
    std::vector <double> tabC;
    matrix finale;
    
    
    for (size_t j(0) ; j < 7 ; ++j)
    {
        
        for ( size_t i(0) ; i < patterns.size() ; ++i)
        {
            
            nombreTot = nombreTot + patterns[i].bScore ;
            if ( patterns[i].site[j] == 'A')
            {
                compteurA = compteurA + patterns[i].bScore ;
            }
            if ( patterns[i].site[j] == 'T')
            {
                compteurT = compteurT + patterns[i].bScore ;
            }
            if ( patterns[i].site[j] == 'G')
            {
                compteurG = compteurG + patterns[i].bScore ;
            }
            if ( patterns[i].site[j] == 'C')
            {
                compteurC = compteurC + patterns[i].bScore ;
            }
        }
        
        tabA.push_back(probas(compteurA , nombreTot));
        tabT.push_back(probas(compteurT , nombreTot));
        tabG.push_back(probas(compteurG , nombreTot));
        tabC.push_back(probas(compteurC , nombreTot));
        compteurA = 0;
        compteurT = 0;
        compteurG = 0;
        compteurC = 0;
        nombreTot = 0;
    }
    
    finale.push_back(tabA);
    finale.push_back(tabT);
    finale.push_back(tabG);
    finale.push_back(tabC);
    
    mx = finale;
    matrix_generation();
    
}
*/

void MatrixProtein::display_PWM_rel()
{
    
    for (size_t i(0); i < pwm_rel.size() ; ++i)
    {
        
        for (size_t j(0) ; j < pwm_rel[i].size() ; ++j)
        {
            
            cout  << setprecision(5) << setw(5) << pwm_rel[i][j] << " | " ;
        }
        
        std::cout << std::endl;
    }
    
}

void MatrixProtein::loadmatrix_fromfile(const string& Data){ // the function stores data from a file containing a PWM or a PSSM assuming that the bases are stored in column and in the ACGT order
    
    int colmn;
    int row;
    vector<double> temp; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open("../test/" + Data);
    
    
    if (file.fail())
    {
        throw runtime_error("Erreur de lecture du fichier de donnée");
    }
    
    string var;
    while (!file.eof())
    {
        while (!file.eof())			//why are you using two while loops that do the same thing?
        {
            file >> var >> ws;
            double dbl = to_Double(var); // allows to get the data in double
            temp.push_back(dbl);
        }
    }
    file.close();
    
    int z(0);
    colmn = 4;
    row = (temp.size())/4;  // here we make the assumption that the base are stored in column
    vector<double> tmp (colmn,0.0 );
    matrix matrix_ (row,tmp); // all the case are initialize at 0.0
    for (int i(0); i < row; ++i) {
        for (int j(0); j < colmn; ++j) {
            matrix_[i][j]=temp[z];
            ++z;
        }
    }
    mx = matrix_;
    matrix_generation();
}

vector<Pattern> MatrixProtein::getPatterns() const
{
    return patterns;
}


void MatrixProtein::matrix_generation()
{
    std::vector<bool> check(2);
    assert (possible(mx)==1);
    check=matrix_status(mx);
    if( (check[0] == 0) & (check[1] == 0))
    {
        pssm_abs=mx;
        swaptorelative(mx);
        pssm_rel = mx;
        swaptopwm(mx);
        pwm_rel = mx;
        mx=pssm_abs;
        swaptopwm(mx);
        pwm_abs=mx;
        
    }
    else if ( (check[0] == 0) & (check[1] == 1))
    {
        pssm_rel=mx;
        swaptoabsolute(mx);
        pssm_abs = mx;
        swaptopwm(mx);
        pwm_abs = mx;
        mx = pssm_rel;
        swaptopwm(mx);
        pwm_rel=mx;
        
    }
    else if ( (check[0] == 1) & (check[1] == 0))
    {
        pwm_abs=mx;
        swaptopssm(mx);
        pssm_abs = mx;
        swaptorelative(mx);
        pssm_rel = mx;
        swaptopwm(mx);
        pwm_rel=mx;
    }
    else
    {
        pwm_rel=mx;
        swaptopssm(mx);
        pssm_rel = mx;
        swaptoabsolute(mx);
        pssm_abs=mx;
        swaptopwm(mx);
        pwm_abs=mx;
    }
}

//Since the computer doesn't like to compute things with -inf, here's a
//function that replaces these values by -100, 2^-100 will be approximated
//by 0.

void MatrixProtein::readjust_values(matrix& mtx)
{
    for (size_t i(0); i < mtx.size(); ++i)
    {
        for (size_t j(0); j < mtx[0].size(); ++j)
        {
            if ( mtx[i][j] < -100 )
            {
                mtx[i][j] = -100;
            }
        }
    }
}

//A function to determine if the transformation from PSSM to PWM was done
//with or without the use of the 0.25. 1 is with, 0 is without.

bool MatrixProtein::which_PWM_to_PSSM(matrix matrice)
{
    assert (possible(matrice));
    
    PWM_to_PSSM(matrice);
    
    if (possible(matrice))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void MatrixProtein::swaptorelative(matrix& mtx)
{
    
    assert (absolute(mtx));
    assert (possible(mtx));
    for(size_t i(0); i < mtx.size(); ++i)
    {
        double max(0.0);
        
        for(size_t j(0); j<4; ++j)
        {
            if ( mtx[i][j] > max )
            {
                max = mtx[i][j];
            }
        }
        
        for(size_t k(0); k<4; ++k)
        {
            mtx[i][k] = mtx[i][k]/max;
        }
    }
}

void MatrixProtein::swaptoabsolute(matrix& mtx)
{
    
    assert (absolute(mtx)==0);
    assert (possible(mtx));
    for(size_t i(0); i < mtx.size(); ++i)
    {
        double total(0.0);
        
        for(size_t j(0); j<4; ++j)
        {
            total += mtx[i][j];
        }
        double x(1/total);
        for(size_t k(0); k<4; ++k)
        {
            mtx[i][k] *= x;
        }
    }
    
}

void MatrixProtein::swaptopssm(matrix& mtx){
    assert (check_if_pmworpssm(mtx));
    for (unsigned int i(0); i < mtx.size() ; ++i)
    {
        for (unsigned int j(0); j < mtx[i].size(); ++j)
        {
            if ( mtx[i][j] == 0 )
            {
                mtx[i][j] = -100;
            }
            else
            {
                mtx[i][j] = log2(mx[i][j]/0.25);
            }
        }
    }
    
    
}
void MatrixProtein::swaptopwm(matrix& mtx){
    assert (check_if_pmworpssm(mtx) == 0);
    for (unsigned int i(0); i < mx.size() ; ++i)
    {
        for (unsigned int j(0); j < mx[i].size(); ++j)
        {
            mx[i][j] = exp2(mx[i][j])*0.25;
        }   // we choose 0.25 as a backgroud because each aa has the same probability to appear randomly
    }
}


vector <bool> MatrixProtein::matrix_status(matrix matrice)
{
    
    vector <bool> a(2);
    
    assert(possible(matrice));
    if (check_if_pmworpssm(matrice))
    {
        a[0] = 1;
        PWM_to_PSSM(matrice);
    }
    else
    {
        a[0] = 0;
    }
    if (absolute(matrice))
    {
        a[1] = 0;
    }
    else
    {
        a[1] = 1;
    }
    
    
    return a;
}

//1 if it's a possible matrix. The idea is simple; we test if the matrix is a PSSM,
//if it's a PWM it may not pass the test. After that we convert the matrix to a PSSM
//and we do the test again. If it doesn't pass the test twice, the matrix is impossible.
//If it's a PSSM it would 100% succeed the first test and probably not the second. If it's a
//PWM it would probably not pass the first test but it would 100% succeed the second.

bool MatrixProtein::possible(matrix matrice){
    
    bool a1(1);
    bool a2(1);
    bool a3(1);
    
    for(size_t i(0); i<matrice.size(); ++i)
    {
        double a(matrice[i][0]);
        double b(matrice[i][1]);
        double c(matrice[i][2]);
        double d(matrice[i][3]);
        
        double sum(a + b + c + d);
        if ( sum < 0.999 or sum > 4.0 )
        {
            a1 = 0;
        }
        if (a != 1 and b != 1 and c != 1 and d != 1 and ( sum < 0.999 or sum > 1.001))
        {
            a1 = 0;
        }
        if (a < 0 or b < 0 or c < 0 or d < 0)
        {
            a1 = 0;
        }
        if (a > 1 or b > 1 or c > 1 or d > 1)
        {
            a1 = 0;
        }
    }
    
    matrix matrice_2;
    matrice_2 = matrice;
    PWM_to_PSSM(matrice_2);
    
    for(size_t j(0); j<matrice_2.size(); ++j)
    {
        double a(matrice_2[j][0]);
        double b(matrice_2[j][1]);
        double c(matrice_2[j][2]);
        double d(matrice_2[j][3]);
        
        double sum(a + b + c + d);
        
        if ( sum < 0.99 or sum > 4.01 )
        {
            a2 = 0;
        }
        if ((a != 1.0) and (b != 1.0) and (c != 1.0) and (d != 1.0) and (sum < 0.999 or sum > 1.001))
        {
            a2 = 0;
        }
    }
    
    PWM_to_PSSM_2(matrice);
    
    for(size_t k(0); k<matrice.size(); ++k)
    {
        double a(matrice[k][0]);
        double b(matrice[k][1]);
        double c(matrice[k][2]);
        double d(matrice[k][3]);
        
        double sum(a + b + c + d);
        
        if ( sum < 0.999 or sum > 4 )
        {
            a3 = 0;
        }
        if ((a != 1.0) and (b != 1.0) and (c != 1.0) and (d != 1.0) and (sum < 0.999 or sum > 1.001))
        {
            a3 = 0;
        }
    }
    
    return (a1 + a2 + a3);
}

void MatrixProtein::PWM_to_PSSM_2(matrix& matrice)
{
    for (size_t i(0); i < matrice.size() ; ++i)
    {
        for (size_t j(0); j < matrice[i].size(); ++j)
        {
            matrice[i][j] = (pow(2, matrice[i][j]));
        }
    }
}

//Used in the possible function. Sometimes, they don't multiply by 0.25 so
//We added a function that does it without to do a third test in the possible function.

void MatrixProtein::PWM_to_PSSM(matrix& matrice)
{
    for (size_t i(0); i < matrice.size() ; ++i)
    {
        for (size_t j(0); j < matrice[i].size(); ++j)
        {
            matrice[i][j] = (pow(2, matrice[i][j]))*0.25;
        }
    }
}

//1 if the matrix is absolute, 0 otherwise. This function checks every
//line of the matrix. A matrix with relative lines and absolute lines
//would be categorised as relative. The input has to be a PSSM. The reason
//why I check every line is because the line ( 1.0 0.0 0.0 0.0 ) where
//1.0 can be at any place, is the same when it's converted to relative.

bool MatrixProtein::absolute(matrix matrice)
{
    assert(possible(matrice));
    
    size_t b(0);
    for(size_t i(0); i<matrice.size(); ++i)
    {
        double a(0);
        for( size_t j(0); j<matrice[0].size(); ++j)
        {
            a += matrice[i][j];
        }
        if (a > 0.999 and a < 1.001)
        {//Since we always compute a limited amount of digits, the sum of
            //the values of each column may be slightly different from 1
            //that's why I didn't put a == 1.0
            ++b;
        }
    }
    if ( b == matrice.size() )
    {
        return 1;
    }
    
    return 0;
}

/*
 Matrix::Matrix(matrix matrice)
 : mx(matrice) {
	readjust_values(mx);
	matrix_generation();}
 */


//1 if it's a PWM. T~he function doesn't launch if the matrix is impossible.
//There's 3 tests, a PSSM can't have a sum of its columns higher than 4
//because the maximum is (1.0 1.0 1.0 1.0) which would be (0.25 0.25 0.25 0.25)
//converted to a relative matrix. If it's an absolute matrix a.k.a
//a matrix that have no 1.0 values, then the sum of its values must be 1.
//The last test is obvious, it can't have any negative values

bool MatrixProtein::check_if_pmworpssm(matrix matrice)
{
    assert (possible(matrice));
    
    for(size_t i(0); i<matrice.size(); ++i)
    {
        double a(matrice[i][0]);
        double b(matrice[i][1]);
        double c(matrice[i][2]);
        double d(matrice[i][3]);
        
        double sum(a + b + c + d);
        if ( sum < 0.999 or sum > 4.0 )
        {
            return 1;
        }
        if (a != 1 and b != 1 and c != 1 and d != 1 and ( sum < 0.999 or sum > 1.001))
        {
            return 1;
        }
        if (a < 0 or b < 0 or c < 0 or d < 0)
        {
            return 1;
        }
    }
    
    return 0;
}

static const char nucleotides[] =		// list of nucleotides
"ATCG"
;

int stringLength = sizeof(nucleotides) - 1;

char MatrixProtein::genRandomChar()  					// Generates a random nucleotide
{
    
    return nucleotides[rand() % stringLength];
}

vector<char> MatrixProtein::seq_generator(double length)
{
    vector <char> sequence;
    for(int i=0; i < length; i++)      		// creates the sequence
    {
        sequence.push_back(genRandomChar());
        
    }
    
    return sequence;
}

double MatrixProtein::set_average(matrix Matrice, double size)
{
    double total(0);
    double average(0);
    
    for (int i(0); i < 7; ++i) {
        vector<char> seq = seq_generator(size);
        try {
            total += get_affinity_score_from_matrix(Matrice, seq);
        } catch (const string& errorMessage) {
            cout << errorMessage;
        }
        
        
    }
    
    
    average = total/7;
    return average;
    
}

double MatrixProtein::get_affinity_score_from_matrix(vector<vector<double> > matrix,vector<char> sequence) // calculates the binding score (double) of a certain sequence based on a Matrix of type PWM
// the function throws an exception which should be catched !
{
    if (matrix.size() < sequence.size()) // checks that the sequence is entirly contained in the matrix
    {
        throw runtime_error("Error : The Matrix doesn't fit the sequence"); //throw exeption
    }
    if (matrix.size()==0) { // checks if the Matrix is configurated
        throw runtime_error("Error : The Matrix isn't configurated "); //throw exeption
    }
    
    double score(1);
    for (unsigned int i(0); i <= sequence.size(); ++i) {
        switch (sequence[i]) {
            case 'A':
                score *= matrix[i][0];
                break;
            case 'C':
                score *= matrix[i][1];
                break;
            case 'G':
                score *= matrix[i][2];
                break;
            case 'T':
                score *= matrix[i][3];
                break;
                
                
            default:
                break;
        }
    }
    return score;
}

void MatrixProtein::setPatterns(vector<Pattern> template_)
{
    patterns = template_;
}

matrix MatrixProtein::getpwm_abs()
{
    return pwm_abs;
}

matrix MatrixProtein::getpssm_abs()
{
    return pssm_abs;
}

void MatrixProtein::get_relevent_site(vector<vector<char>> Input, int set, int threshold)
{
    vector<vector<char>> before_threshold;
    for (size_t i(0); i < Input.size(); ++i)
    {
        for (size_t j(0); j < Input[i].size() - set; ++j) {
            vector<char> tmp;
            for (size_t z(0); z < set; ++z)
            {
                tmp.push_back(Input[i][j+z]); // construct small vector of size n (given by the User)
            }
            before_threshold.push_back(tmp);
        }
    }
    fillVectorPatterns(before_threshold, threshold); //gives all the possible combination of size n to test the affinity score
}

//COMMENTED AS DOES NOT WORK YET + IS IT USEFUL ?
/*vector<Pattern> MatrixProtein::findPatterns(string matFile, double threshold) //need to be tested
{
	
					//we map matrix colums (in the assumed ACGT order) to a character so that we can concatenate a motif
		map<char, size_t> column2nuc{{0,'A'}, {0, 'a'},
									{1,'C'}, {1, 'c'},
									{2,'G'}, {2, 'g'},
									{3,'T'}, {3, 't'},
									};

	loadmatrix_fromfile(matFile);
	vector<Pattern> listOfMotifs;
	Pattern motif;
	
	double temporaryScore;
	double score(1.0);
	vector<char> buildingMotif;
	double compteur(0.0); 
	
	
	do {
		//startagain and keeplooking are two labels
		startagain:
		
		for(size_t row(0); row < pwm_abs.size(); ++row)				//for each row we look for a maximum value...
		{
			keeplooking:
			
			temporaryScore = pwm_abs[row][0];						//...to find score that is first the first value encountered
			for(size_t column(0); column < (pwm_abs[row]).size(); ++column)  // careful here : here column is an int as we used it in the map
			{
				if (temporaryScore < pwm_abs[row][column])
				{
					temporaryScore = pwm_abs[row][column];	
				} 
				
				score *=temporaryScore;
				buildingMotif.push_back(column2nuc[column]);
				
				if(score < threshold) 
				{
					score = 1.0;
					goto startagain;					//if we reached a score that is below the threshold, we break the loop to start again
				} else
				{
					goto keeplooking;					//question : will it increment row with ++row ? do I have to add before goto ?
				}
			
			}
		}
		motif.bScore = score;							//we set a newfound motif
		motif.site = buildingMotif;
		listOfMotifs.push_back(motif);					//we had this motif to the list
		++compteur;
	
	} while (compteur != pow(4.0, 7.0));  //while we haven't reached all motif possibilites(4^7), we keep trying to find some
	
	return listOfMotifs;
}*/



void calcul( vector<char> tab_, vector<double>& tabA_, vector<double>& tabT_, vector<double>& tabG_, vector<double>& tabC_)
{
	
	for ( size_t w(0) ; w < tab_.size() ; ++w) // la je compte le nombre de A, T, G, C à chaque position
			{
					
				if ( tab_[w] == 'A')
				{
					++ tabA_[w] ;
				}
				
				if ( tab_[w] == 'T')
				{
					++ tabT_[w] ;
				}
				
				if ( tab_[w] == 'G')
				{	
					++ tabG_[w];
				}	
				
				if ( tab_[w]  == 'C')
				{
					++ tabC_[w] ;
				}
			}
}

void RemetZero( vector<double>& vec_)
{
	for ( size_t i(0); i < vec_.size(); ++i)
	{
		vec_[i] = 0; 
	}
}

void calculScore(vector<vector<double>> finale_, vector<double>& score_, int j_, vector<char> tab2_)
{
		int indice(0); 
		for (size_t l(0) ; l < finale_[0].size() ; ++l) // on traduit la lettre du tableau temporaire en indice du tableau finale
			{ 
				if (tab2_[l] == 'A')
				{
					indice = 0; 
				} 
		
				else if (tab2_[l] == 'T')
				{
					indice = 1; 
				}
						
				else if (tab2_[l] == 'G')
				{
					indice = 2; 
				}
							
				else if (tab2_[l] == 'C')
				{
					indice = 3; 
				}
							
							
				score_[j_] += finale_[indice][l]; // on calcul le "score" pour chaque combinaison de taille longueur_motif possible
			}
}


void FindMotif(vector<vector<double>> finale_ , vector<Pattern> patterns_, int longueur_motif_, vector<vector<vector<char>>>& best_seqs_)
{
	 int index(0);
	 int b(0);
	 double nf(0); 
	 vector<char> tab2(longueur_motif_);
	 
	for( size_t i(0) ; i < patterns_.size() ; ++i) // on va dans chaque motif de la liste
		{
			b = patterns_[i].site.size() - longueur_motif_ + 1; // nombre de combinaison dans la séquence de taille  "longueur_motif". 
			vector<double> score(b, 0);
			vector<int> indices(0); 
			vector<vector<char>> vec1(0);
			for (size_t j(0) ; j < b ; ++j) 
				{
					if (j == 0) // soit je remplit le tableau temporaire pour la premiere fois 
					{
						for ( size_t k(0) ; k < longueur_motif_ ; ++k) 
						{
							tab2[k] = patterns_[i].site[k]; 
						}
					}
					
					else { // soit j'enleve la premiere case du tableu temporaire et je rajoute la derniere. On décale le motif  
							tab2.erase(tab2.begin()); 
							tab2.push_back(patterns_[i].site[longueur_motif_ + j - 1]); 
						 }
						 
						 
						 
						calculScore(finale_, score, j, tab2); 
						
					
					 if( score[j] == nf) // pour rajouter un autre indice si y en a dejà un avec le même score
						{  
							indices.push_back(j); 
						}
				 
					 if( score[j] > nf) // si le score est plus grand, enlever les autres dans le vecteur d'indice et mettre celui la à la place
						{
							nf = score[j]; 
							if (j > 0)
							{
								indices.clear(); 
							}
							indices.push_back(j);
						}
						
				}
				
			for (size_t t(0); t < indices.size() ; ++t)
			{ 
				index = indices[t]; 
				vector<char>vec2(0); 
				for ( size_t d(index) ; d < index + longueur_motif_ ; ++d)
				{
					vec2.push_back(patterns_[i].site[d]);  // on met dans un tableau les bases du motif de la position retrnue par indice de longueur "longueur_mtif"
				}
				vec1.push_back(vec2); // on le met dans le tableau pour chaque motif séparement 
				  
			}
			
			best_seqs_.push_back(vec1); // tous les motifs ensemble
			nf = 0;  
		}
}



void LoadMatrixFromListOfSites(int longueur_motif)
{
//patterns est l'attribut de matrixprotein. 
//--------------------------------------------------------------------------------------------------------------------- 
    // LISTE D'INITIALISATION
    
    vector<char> tab(longueur_motif); // les tableaux temporaires
    int somme(0); 
    std::vector <double> tabA(longueur_motif,0);
    std::vector <double> tabT(longueur_motif,0);
    std::vector <double> tabG(longueur_motif,0);
    std::vector <double> tabC(longueur_motif,0);
    vector < vector<double> > finale; // premiere matrice de l'algorithme
    vector <vector<double> > PWM; // La matrice à renvoyé au final
    int a(0); 
	vector<vector<vector<char>>> best_seqs(0); 
	
//---------------------------------------------------------------------------------------------------------------------
	//REMPLIR MATRICE FINALE
    
    for( size_t i(0) ; i < patterns.size() ; ++i) // on va dans chaque motif de la liste
    {
		a = patterns[i].site.size() - longueur_motif + 1; // nombre de combinaison dans la séquence de taille  "longueur_motif". 
		
		somme = somme + a; // pour calculer le nombre total (dans toute la liste de site) de combinaison de taille  "longueur_motif".
		
		for (size_t j(0) ; j < a ; ++j) 
		{
			if (j == 0) // soit je remplit le tableau temporaire pour la premiere fois 
			{
				for ( size_t k(0) ; k < longueur_motif ; ++k) 
				{
					tab[k] = patterns[i].site[k]; 
				}
			}
			
			else { // soit j'enleve la premiere case du tableu temporaire et je rajoute la derniere. On décale le motif  
					tab.erase(tab.begin()); 
					tab.push_back(patterns[i].site[longueur_motif + j - 1]); 
				 }
					
			calcul(tab, tabA, tabT, tabG ,tabC); 
		}
        
    }
    
    // on ajoute dans la matrice finale le nombre de chaque lettre a chaque position du motif
    finale.push_back(tabA);//pour A
    finale.push_back(tabT);//pour T
    finale.push_back(tabG);//...
    finale.push_back(tabC);
    
    
//---------------------------------------------------------------------------------------------------------------------    
    // DIVISION PAR LA SOMME   
    for ( size_t i(0) ; i < finale.size() ; ++i) 
    {
		for (size_t j(0); j < longueur_motif; ++j)
		{
			finale[i][j] = finale[i][j]/somme; 
		}
	}
	
//---------------------------------------------------------------------------------------------------------------------
	//TROUVER LE MEILLEUR MOTIF DE CHAQUE SITE
	
	FindMotif(finale, patterns, longueur_motif, best_seqs); 
	
//----------------------------------

	somme = 0; 
	vector<vector<double>> copiePWM (4, vector<double>(longueur_motif)); 
	
	RemetZero(tabA); 
	RemetZero(tabT);
	RemetZero(tabG);
	RemetZero(tabC);
	
	for ( size_t i(0); i < best_seqs.size() ; ++i)
	{
		for ( size_t j(0); j < best_seqs[i].size() ; ++j)
		{
			calcul(best_seqs[i][j], tabA, tabT, tabG, tabC); 
			++somme; 
		}
	}
	PWM.clear(); 
	PWM.push_back(tabA);//pour A
    PWM.push_back(tabT);//pour T
    PWM.push_back(tabG);//...
    PWM.push_back(tabC);  
	
    
	for ( size_t i(0) ; i < PWM.size() ; ++i) 
    {
		for (size_t j(0); j < longueur_motif; ++j)
		{
			PWM[i][j] = PWM[i][j]/somme; 
			// on divise par la somme pour avoir la probabilité 
		}
	} 

	do{
		
		copiePWM = PWM;  
		best_seqs.clear(); 
		FindMotif(PWM, patterns, longueur_motif, best_seqs);   // recalcul du meilleur motif 
		somme = 0; 
		RemetZero(tabA); 
		RemetZero(tabT);
		RemetZero(tabG);
		RemetZero(tabC);
		
		for ( size_t i(0); i < best_seqs.size() ; ++i)
		{
			for ( size_t j(0); j < best_seqs[i].size() ; ++j)
			{
				calcul(best_seqs[i][j], tabA, tabT, tabG, tabC); // remodification de tabA, tabT, tabG, tabC; 
				++ somme; 
			}
		}
		
		PWM.clear(); 
		PWM.push_back(tabA);//pour A
		PWM.push_back(tabT);//pour T
		PWM.push_back(tabG);//...
		PWM.push_back(tabC); 
		
		 for ( size_t i(0) ; i < 4 ; ++i) 
		{
			
			for (size_t j(0); j < longueur_motif; ++j)
			{
				PWM[i][j] = PWM[i][j]/somme; 
			}
		}
	} while (PWM != copiePWM); 
	
	// Ma matrice est pour l'instant sous la forme: 1 2 3 ....
												//A
												//T
												//G
												//C
												
	// je dois encore remettre PWM dans l'affichage que je veux: A C G T
															//1
															//2
															//3
															//...

	vector <vector<double> > PWMfinale(0); 
	for ( size_t i(0); i < longueur_motif ; ++i)
	{
		vector<double> tabpos; 
		tabpos.push_back(PWM[0][i]);
		tabpos.push_back(PWM[3][i]);
		tabpos.push_back(PWM[2][i]);
		tabpos.push_back(PWM[1][i]);
		PWMfinale.push_back(tabpos); 
		tabpos.clear(); 
	}
	
	mx = PWMfinale; 
	 	
}


