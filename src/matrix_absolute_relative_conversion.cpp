#include <iostream>
#include <vector>
#include <string>

using namespace std;

void conversion_of_matrix(vector<vector<double> >& matrix);
void afficher_matrice(vector<vector<double> > matrice);
void reverse_conversion_of_matrix(vector<vector<double> >& matrix);

int main()
{
	vector<vector<double> > a(2, vector<double> (4));
	a[0][0] = 0.184211;
	a[0][1] = 0.000000;
	a[0][2] = 0.105263;
	a[0][3] = 0.710526;
	a[1][0] = 0.25;
	a[1][1] = 0.25;
	a[1][2] = 0.25;
	a[1][3] = 0.25;
	
	vector<vector<double> > b(2, vector<double> (4));
	b[0][0] = 0.259260;
	b[0][1] = 0.000000;
	b[0][2] = 0.148148;
	b[0][3] = 1.000000;
	b[1][0] = 1.0;
	b[1][1] = 1.0;
	b[1][2] = 1.0;
	b[1][3] = 1.0;
	
	conversion_of_matrix(a);
	afficher_matrice(a);
	
	reverse_conversion_of_matrix(b);
	afficher_matrice(b);
	
	return 0;
}

void conversion_of_matrix(vector<vector<double> >& matrix)
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

void reverse_conversion_of_matrix(vector<vector<double> >& matrix)
{
	size_t nb_lignes(matrix.size());
	size_t nb_colonnes(matrix[0].size());
	
	for(size_t i(0); i<nb_lignes; ++i)
	{	
		double total(0);
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
	
	


	
	

