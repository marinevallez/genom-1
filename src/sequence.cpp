#include <iostream>
#include <fstream>
#include "sequence.hpp"
using namespace std;

void Sequence::loadFile(string fileName) 
{
	
	string seq1_, seq2_;						//first and second sequences;
	
	ifstream file;								//the file we are going to read
	file.open("../test/" + fileName); 			//since out text files are in the test folder, we need to include a path to it
	
	if(file.fail())								// if it didnt open -> show an error 
	{
		cerr << "File could not be opened!" << endl;
	}
	else
	{
		string skip; 		//we wanna skip the first line
		
		file >> skip;		//skip the first line
		file >> seq1_;		//get the first sequence
		file >> skip;		//skip the second useless part
		file >> seq2_;		//get the second sequence
			
		file.close();
	}
	
	for(char& c : seq1_) 	//now we convert strings into tables of characters to facilitate the work later
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
	for(auto c : seq1)	//we display the sequence character by character
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


vector<int> Sequence::searchMotif(string subStr, int seqNb) const
{
	vector<int> output; 	//the vector of all matching positions
	
	vector<char> subStr_;	//used to convert a string into a table of char
	
	for(char& c : subStr)	//we convert the substring into a table of characters as well
	{
		subStr_.push_back(c);
	}
		
	vector<char> seqCopy; //used to make the method more versatile
	
	switch(seqNb)		//we set seqCopy to be the sequence we want to search (either seq1 or seq2) depending on the parameter entered
	{
		case 1: 
			seqCopy = seq1;
			break;
		case 2:			
			seqCopy = seq2;
			break;
		default:
			break;
	}
	
	for(size_t i(0); i < seqCopy.size() - 6; ++i)	//we check the whole sequence 
			{
					bool matchingCondition(seqCopy[i] == subStr_[0] and seqCopy[i+1]== subStr_[1] and seqCopy[i+2]== subStr_[2]
						and seqCopy[i+3]== subStr_[3] and seqCopy[i+4] == subStr_[4] and seqCopy[i+5] == subStr_[5]
						and seqCopy[i+6] == subStr_[6]); //we know that 7 nucleotides in a row have to match
					
					if(matchingCondition)		//if 7 nucleotides in a row do match then
					{
						output.push_back(i+1);	//we add the position OF THE FIRST MATCHING NUCLEOTIDE to the table 
					}		
			}
	
	return output;
};

/*int main()
{
	Sequence seq_;
	seq_.loadFile("promoters.fasta");
	seq_.display();
	for(auto c : seq_.searchMotif("TTCCCCA", 2))
	{
		cout << c << " ";
	}
	return 0;
}
*/
