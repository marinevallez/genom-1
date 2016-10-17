#include <vector>
#include <iostream>
#include <iomanip>

double probas(double n, double tot) 
{
	return n/tot; 
}

struct Motif 
{
	std::vector<char> lettres; 
	double binding_score; 
};

int main()
{
	typedef std::vector< Motif>Tableau; 
	typedef std::vector<std::vector<double>> Matrice;
	Tableau sequences;
	double compteurA(0); 
	double compteurT(0);
	double compteurG(0);
	double compteurC(0);
	double nombreTot(0); 
	std::vector <double> tabA;    
	std::vector <double> tabT;
	std::vector <double> tabG;
	std::vector <double> tabC;
	Matrice finale; 
	
	Motif m1 = {{'A', 'T', 'G', 'C', 'T', 'G', 'T'}, 0.8};
	sequences.push_back( m1) ;
	Motif m2 = {{'T', 'C', 'A', 'C', 'T','C','C'}, 0.2};
	sequences.push_back( m2) ;
	Motif m3 = {{'T', 'T', 'G', 'C' , 'C','A', 'C'}, 0.4};
	sequences.push_back( m3) ;
	Motif m4 = {{ 'A' , 'A', 'C', 'G', 'G','T', 'C'}, 0.7};
	sequences.push_back( m4) ;
	

	for (size_t j(0) ; j < 7 ; ++j) 
	{

		for ( size_t i(0) ; i < sequences.size() ; ++i) 
		{
			
			nombreTot = nombreTot + sequences[i].binding_score ;
		  if ( sequences[i].lettres[j] == 'A')
			{ 
				compteurA = compteurA + sequences[i].binding_score ; 
			}
		  if ( sequences[i].lettres[j] == 'T')
			{ 
				compteurT = compteurT + sequences[i].binding_score ; 
			}
		  if ( sequences[i].lettres[j] == 'G')
			{ 
				compteurG = compteurG + sequences[i].binding_score ; 
			}
		  if ( sequences[i].lettres[j] == 'C')
			{ 
				compteurC = compteurC + sequences[i].binding_score ;  
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
	
		for ( int i(0); i < finale.size() ; ++i)
		{
			
			for (int j(0) ; j < finale[i].size() ; ++j)
			{
				
				std::cout  << std::setprecision(5) << std::setw(5) << finale[i][j] << " | " ; 
			} 
			 
			std::cout << std::endl; 
		}
	return 0;  

}
