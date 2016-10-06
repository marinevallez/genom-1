#include <iostream>
#include <vector>
#include <string>

//A function that takes a sequence s and an int n and gives how many time on average
//s would appear in a n long sequence.

using namespace std;

double appearance_probability(string sequence, int longueur);

const double at (0.2); //We know from a research on internet that this is the
const double gc (0.3); //probabilty of appearance for these nucleic acid

int main()
{
	double a;
	
	a = appearance_probability("gattaca", 100);
	
	cout<<a<<endl;
	
	return 0;
}

double appearance_probability(string sequence, int longueur)
{
	double proba (1.0);
	int longueur_2(sequence.length());
	
	for (size_t i(0); i<sequence.length(); ++i)
	{
		if (sequence[i] != 'a' or 't' or 'c' or 'g')
		{
			--longueur_2; //In case the input contains other letters, 
			              //it won't take them in count.
		}
	}
			
	double position_possible(longueur - longueur_2 + 1);
	
	if (position_possible < 0)
	{
		position_possible = 0;
	}
	
	for(size_t i(0); i<sequence.length(); ++i)
	{
		if ((sequence[i] == 'a') or (sequence[i] == 't'))
		{
			proba*=at;
		}
		if ((sequence[i] == 'g') or (sequence[i] == 'c'))
		{
			proba*=gc;
		}
	} //We just multiply the probabilty of appearance of the sequence by the
	  //number of places it can appear.
	
	proba = proba*position_possible; 
	
	return proba;
}
