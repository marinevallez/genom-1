#include <iostream>
#include <vector>

//This code gives two functions, one that converts a matrix to its absolute matrix and the other that converts
//a matrix to its relative matrix. The two other function, afficher_matrix and creer_matrice are just here
//to show if the two others are working correctly.

using namespace std;

void relative_matrix(vector<vector<double> >& matrix);
void absolute_matrix(vector<vector<double> >& matrix);
void afficher_matrice(vector<vector<double> > matrice);
vector<vector<double> > creer_matrice();

int main()
{	
	vector <vector<double> > M;
		
	M = creer_matrice();
	
	cout<<"La matrice relative correspondante : "<<endl;
	
	relative_matrix(M);
	afficher_matrice(M);
	
	cout<<"La matrice remise en absolue : "<<endl;
	//ça ne marche que si la somme des lignes est égale à 1
	absolute_matrix(M);
	afficher_matrice(M);
	
	return 0;
}

//----------------------------------------------------------------------

void relative_matrix(vector<vector<double> >& matrix)
{
	size_t nb_lignes(matrix.size());
	size_t nb_colonnes(matrix[0].size());
	
	for(size_t i(0); i<nb_lignes; ++i)
	{
		double max(0.0);
		
		for(size_t j(0); j<nb_colonnes; ++j) 
		{                                    
			if ( matrix[i][j] > max )
			{
				max = matrix[i][j];
			}
		}
		
		for(size_t k(0); k<nb_colonnes; ++k)
		{
			matrix[i][k] = matrix[i][k]/max;
		}
	}
}

//----------------------------------------------------------------------

void absolute_matrix(vector<vector<double> >& matrix)
{
	size_t nb_lignes(matrix.size());
	size_t nb_colonnes(matrix[0].size());
	
	for(size_t i(0); i<nb_lignes; ++i)
	{	
		double total(0.0);
		
		for(size_t j(0); j<nb_colonnes; ++j)
		{
			total += matrix[i][j];
		}
		double x(1/total);
		for(size_t k(0); k<nb_colonnes; ++k)
		{
			matrix[i][k] *= x;
		}
	}
}
	
//----------------------------------------------------------------------	
			
void afficher_matrice(vector<vector<double> > matrice)
{
	for(size_t i(0); i<matrice.size(); ++i)
	{
		for( size_t j(0); j<matrice[0].size(); ++j)
		{
			cout<<matrice[i][j]<<" ";
		}
		cout<<endl;
	}
} 

//----------------------------------------------------------------------

vector<vector<double> > creer_matrice()
{
	size_t nb_lignes (0);
	
	do
	{
	  cout<<"Nombre de lignes : ";
	  cin>>nb_lignes;
    }
    while (nb_lignes<1);
    
	size_t nb_colonnes (0);
	
	do
	{
	cout<<"Nombre de colonnes : ";
	cin>>nb_colonnes;
    }
    while (nb_colonnes<1);
    
	vector <vector<double> > M (nb_lignes, vector<double> (nb_colonnes));
	for ( size_t i(0); i<nb_lignes ; ++i)
	{
		for ( size_t j(0); j<nb_colonnes; ++j)
		{
			cout<<"M["<<i+1<<","<<j+1<<"]=";
			cin>>M[i][j];
		}
	}
	return M;
}	
	


	
	


