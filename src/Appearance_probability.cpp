#include <iostream>
#include <vector>
#include <string>

using namespace std;

double appearance_probability(string sequence, int longueur);

const double at (0.2);
const double gc (0.3);

int main()
{
	double a;
	
	a = appearance_probability("gtagatgatga", 9999);
	
	cout<<a<<endl;
	
	return 0;
}

double appearance_probability(string sequence, int longueur)
{
	double proba (1.0);
	double position_possible(longueur-sequence.length()+1);
	
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
	}
	
	proba = proba*position_possible;
	
	return proba;
}
