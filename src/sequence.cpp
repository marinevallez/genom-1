#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

class Sequence {
	
private:
	
	vector<char> seq1;	//the two sequences we are working with
	vector<char> seq2;
	
public:

	void loadSeq();					// read two sequences from a file
	size_t searchSeq(string) const;	//find the position of subsequence within the bigger sequence
	void display();					//display the two sequences
	
};


void Sequence::loadSeq()
{
	const string fileName("promoters.fasta.txt"); //name of the file
	string seq1_;	//first and second sequences
	string seq2_;
	
	ifstream file(fileName);	//the file we are reading
	
	if(file.fail())		// if it didnt open -> show an error 
	{
		cerr << "File could not be opened!" << endl;
	}
	else
	{
		string skip; 		//we wanna skip the first line
		file >> skip;		
		file >> seq1_;		//get the first sequence
		file >> skip;		//skip the second useless part
		file >> seq2_;		//get the second sequence
			
		
		file.close();
	}
	
	for(char& c : seq1_) 	//now we convert strings into tables of characters to facilitate the work
	{
		seq1.push_back(c);
	}
	
	for(char& c : seq2_)
	{
		seq2.push_back(c);
	}
	
	
};

void Sequence::display()
{
	for(auto c : seq1)	//we display the sequence character by chacter
	{
		cout << c;
	}
	
	cout << endl << endl;
	
	for(auto c : seq2)
	{
		cout << c;
	}
	
	cout << endl << endl;
};

size_t Sequence::searchSeq(string subStr) const
{
	vector<char> subStr_;
	
	for(char& c : subStr)	//we convert the substring into a table of characters as well
	{
		subStr_.push_back(c);
	}
	
	
	for(size_t i(0); i < seq1.size(); i++)
	{
		if(seq1[i] == subStr_[i] and seq1[i+1] == subStr_[i+1]
			and seq1[i+2] == subStr_[i+2] and seq1[i+3] == subStr_[i+3])
		{
			return i+1;
		}
	}
	
	return 456;
			
};

int main()
{
	Sequence seq_;
	seq_.loadSeq();
	seq_.display();
	cout << seq_.searchSeq("CCCA");
	return 0;
}
