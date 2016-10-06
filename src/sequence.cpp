#include <iostream>
#include <fstream>
#include "sequence.hpp"
using namespace std;


void Sequence::loadSeq()
{
	const string fileName("promoters.fasta.txt"); //name of the file
	string seq1_;								  //first and second sequences
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

size_t Sequence::searchSeq(string subStr) const		//so far this method only searches the first sequence and only 4 letters
{
	vector<char> subStr_;
	
	for(char& c : subStr)	//we convert the substring into a table of characters as well
	{
		subStr_.push_back(c);
	}
		
	for(size_t i(0); i < seq1.size() - 3; ++i)	//we check the whole sequence 
	{
			if(seq1[i] == subStr_[0] and seq1[i+1]== subStr_[1] and seq1[i+2]== subStr_[2]
				and seq1[i+3]== subStr_[3]) {return i+1;}		//if the we find that the position where ther sequences match then we return it
	}
	
};

vector<int> Sequence::searchSeq_(string subStr) const
{
	vector<int> output; 	//the vector of all matching positions
	
	vector<char> subStr_;	//used to convert a string into a table of char
	
	for(char& c : subStr)	//we convert the substring into a table of characters as well
	{
		subStr_.push_back(c);
	}
		
	for(size_t i(0); i < seq1.size() - 3; ++i)	//we check the whole sequence 
	{
			if(seq1[i] == subStr_[0] and seq1[i+1]== subStr_[1] and seq1[i+2]== subStr_[2]
				and seq1[i+3]== subStr_[3]) {output.push_back(i+1);}		//if the we find that the position where ther sequences match then we add its value
																			//to the table
	}
	
	return output;
};

int main()
{
	Sequence seq_;
	seq_.loadSeq();
	seq_.display();
	for(auto c : seq_.searchSeq_("CCCC"))
	{
		cout << c << " ";
	}
	return 0;
}
