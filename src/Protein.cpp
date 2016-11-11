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
	pattern.listOfSites = site;
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
    
    /* Motif m1 = {{'A', 'T', 'G', 'C', 'T', 'G', 'T'}, 0.8};
     sequences.push_back( m1) ;
     Motif m2 = {{'T', 'C', 'A', 'C', 'T','C','C'}, 0.2};
     sequences.push_back( m2) ;
     Motif m3 = {{'T', 'T', 'G', 'C' , 'C','A', 'C'}, 0.4};
     sequences.push_back( m3) ;
     Motif m4 = {{ 'A' , 'A', 'C', 'G', 'G','T', 'C'}, 0.7};
     sequences.push_back( m4) ;
     */
    
    
    for (size_t j(0) ; j < 7 ; ++j)
    {
        
        for ( size_t i(0) ; i < patterns.size() ; ++i)
        {
            
            nombreTot = nombreTot + patterns[i].Bscore ;
            if ( patterns[i].site[j] == 'A')
            {
                compteurA = compteurA + patterns[i].Bscore ;
            }
            if ( patterns[i].listOfSites[j] == 'T')
            {
                compteurT = compteurT + patterns[i].Bscore ;
            }
            if ( patterns[i].listOfSites[j] == 'G')
            {
                compteurG = compteurG + patterns[i].Bscore ;
            }
            if ( patterns[i].listOfSites[j] == 'C')
            {
                compteurC = compteurC + patterns[i].Bscore ;
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
    
    mtrx.mx = finale; 
    mtrx.matrix_generation(); 
    
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

matrix Protein::loadmatrix(string Data){ // the function stores data from a file containing a PWM or a PSSM assuming that the bases are stored in column and in the ACGT order

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
    setrw(row);
    return mtrix;
}

vector<Pattern> Protein::getPatterns()
{
	return patterns; 
}
