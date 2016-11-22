#include "Protein.hpp"
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;

double To_double( const string& string ) // allows a convertion from string to double
{
    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}

// ==============================================================================================================CONSTRUCTEUR ET DESTRUCTEUR

Protein::Protein()
{
    
}

Protein::~Protein()
{
    
}

// ==============================================================================================================METHODES

static const char nucleotides[] =		// list of nucleotides
"ATCG"
;

int stringLength = sizeof(nucleotides) - 1;

char Protein::genRandomChar()  					// Generates a random nucleotide
{
    
    return nucleotides[rand() % stringLength];
}

vector<char> Protein::seq_generator(double length)
{
    vector <char> sequence;
    for(int i=0; i < length; i++)      		// creates the sequence
    {
        sequence.push_back(genRandomChar());
        
    }
    
    return sequence;
}

double Protein::set_average(matrix Matrice, double size)
{
    double total(0);
    double average(0);
    
    for (int i(0); i < 7; ++i) {
        vector<char> seq = seq_generator(size);
        try {
            total += get_afinity_score_from_matrix(Matrice, seq);
        } catch (runtime_error message) {
            cout << message.what();
        }
        
        
    }
    
    
    average = total/7;
    return average;
    
}



void Protein::fillVectorPatterns(vector<vector<char>> sites, double threshold )// Should we ask if the user rather have a binding score calculated from a relative matrix ? woukd that make a difference (NB this function calculates the BS from a PWM absolute) this takes all the binding site for a given prot(in the fasta file)

{
    double binding_score(0);
    if (threshold == 0) {
        threshold = set_average(pwm_abs, sites[0].size());
    }
    
    for (size_t i(0); i < sites.size(); ++i) {
        
        try {
            binding_score = get_afinity_score_from_matrix(pwm_abs, sites[i]);
        } catch (runtime_error message) {
            cout << message.what();
        }
        
        
        if(binding_score > threshold)
        {
            Pattern pattern ({binding_score, sites[i]});
            patterns.push_back(pattern);
            
        }
        
        
    }
}


double Protein::probas(double n, double tot)
{
    return n/tot;
}



void Protein::loadmatrix_fromscore()
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
    
    mtrx.mx = finale; // du coup ca remplit quoi ?
    mtrx.matrix_generation();
    
}

void Protein::check_if_empty(string Data)
{	
	bool a;
	
	matrix b;
	
	b = mtrx.mx;
	
	if ( b.empty() == 0 )
	{
		cout<<"There's already a matrix loaded, do you want to replace it ?"<<endl;
		cout<<"Press 1 for yes or 0 for no"<<endl;
		cin>>a;
		if (a)
		{
			loadmatrix(Data);
		}
	}
	else
	{
		loadmatrix(Data);
	}
}
			

void Protein::display_PWM()
{
    
    for ( int i(0); i < mtrx.size() ; ++i)
    {
        
        for (int j(0) ; j < mtrx[i].size() ; ++j)
        {
            
            cout  << setprecision(5) << setw(5) << mtrx[i][j] << " | " ;
        }
        
        std::cout << std::endl;
    }
    
}


matrix Protein::loadmatrix_fromfile(std::string Data){ // the function stores data from a file containing a PWM or a PSSM assuming that the bases are stored in column and in the ACGT order
    
    int colmn;
    int row;
    string fichier (Data);
    vector<double> temp; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open(fichier);
    
    
    if (file.fail())
    {
        throw std::runtime_error("Erreur de lecture du fichier de donnée ");
    }
    
    string var;
    while (!file.eof())
    {
        while (!file.eof())
        {
            file >> var >> ws;
            double dbl = To_double(var); // allows to get the data in double
            temp.push_back(dbl);
        }
    }
    file.close();
    
    int z(0);
    colmn = 4;
    row = (temp.size())/4;  // here we make the assumption that the base are stored in column
    vector<double> tmp (colmn,0.0 );
    matrix mtrix (row,tmp); // all the case are initialize at 0.0
    for (int i(0); i < row; ++i) {
        for (int j(0); j < colmn; ++j) {
            mtrix[i][j]=temp[z];
            ++z;
        }
    }
    //setrw(row);
    return mtrix;
}

vector<Pattern> Protein::getPatterns()
{
    return patterns;
    
}

double Protein::get_afinity_score_from_matrix(vector<vector<double>> matrix,vector<char> sequence) // calculates the binding score (double) of a certain sequence based on a Matrix of type PWM
// the function throws an exception which should be catched !
{
    if (matrix.size() < sequence.size()) // checks that the sequence is entirly contained in the matrix
    {
        throw std::runtime_error("Error : The Matrix doesn't fit the sequence"); //throw exeption
    }
    if (matrix.size()==0) { // checks if the Matrix is configurated
        throw std::runtime_error("Error : The Matrix isn't configurated "); //throw exeption
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

