#include <iostream>
#include <fstream>
#include "sequence.hpp"
using namespace std;

//Sequence::Sequence() {}

//Sequence::~Sequence() {}

/*!
 * The chromosomeNb method extracts the chromosome number (e.g. from 1 to 23 in humans) from .fasta file
 */
 
string chromosomeNb(const string& str)  //a function that extracts the chromosome number from .fasta
{
	string chr;
	
	if(str[5] == '|')
	{
		chr = str.substr(1,3) + " " + str[4];
	}
	else
	{
		chr = str.substr(1,3) + " " + str.substr(4,2);
	}
	
	return chr;	
} 


/*!
 * The outputSites method outputs a file containing sequence numbers (from a .fasta file) and the position on which a motif is found, 
 * as well as the direction of the reading frame (+ means forward, - means reverse).
 */
 
void Sequence::outputSites(string motif)	//a method that outputs a file with seq nb and a position where the motif was found
{
	ofstream file;
	file.open("../test/sites.bed");	//we create a file where we will output results
	if(file.fail())
	{
		cerr << "The file couldn't be created!";
	}
	else
	{
		vector<vector<char>> seq = quickRead("promoters.fasta");	//we get read both sequences from .fasta
		vector<vector<char>> seqComp = {giveComplementarySeq(seq[0]), giveComplementarySeq(seq[1])};	//then we get their compliments
		for(size_t i(0); i < 2; ++i)
		{
			vector<int> pos = searchMotif(seq[i], motif);	//we find where the motif appears in the provided sequences
			if(!pos.empty())	//if it appears then we save it to file
			{
				file << "seq" << i+1 << " " << pos[0] << " to " << pos[0]+7 <<  " + " << motif << endl;
			} 
			pos = searchMotif(seqComp[i], motif);	//then we check the compliment of those sequences
			if(!pos.empty())	//if it appears there then we also save it
			{
				file << "seq" << i+1 << " " << 400 - pos[0]<< " to " << 400 - pos[0] - 7 << " - "<< motif << endl;
			}
		}
		file.close();		//we close the file
	}
}

/*!
 * The quickRead method reads a given .fasta file and outputs on the terminal both sequences of the file, withough saving them somewhere for efficiency.
 */
 
vector<vector<char>> Sequence::quickRead(string fileName) const	//a method that reads .fasta and outputs both sequences instead of saving them
{
	vector<vector<char>> output;
	string seq1_(""), seq2_("");
	
	ifstream file;								//the file we are going to read
	file.open("../test/" + fileName); 			//since out text files are in the test folder, we need to include a path to it
	
	if(file.fail())								// if it didnt open -> show an error 
	{
		cerr << "File could not be opened!" << endl;
	}
	else
	{
		
		
		string text;
		string chr(""); 		
		
		file >> text;		//skip the first line
			
		//cout << chromosomeNb(text) << " ";	
		
		file >> seq1_;		//get the first sequence
		file >> text;		//skip the second useless part
		
		//cout << chromosomeNb(text) << " ";
		
		file >> seq2_;		//get the second sequence
			
		file.close();
	}
	
	vector<char> v1, v2;
	
	for(char& c : seq1_) 	//now we convert strings into tables of characters to facilitate the work later
	{
		v1.push_back(c);
	}
	
	output.push_back(v1);
	
	for(char& c : seq2_)
	{
		v2.push_back(c);
	}
	
	output.push_back(v2);
	
	return output;
}

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
		
		string text;
		string chr(""); 		
		
		file >> text;		//skip the first line
			
		cout << chromosomeNb(text) << " ";	
		
		file >> seq1_;		//get the first sequence
		file >> text;		//skip the second useless part
		
		cout << chromosomeNb(text) << " ";
		
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

/*!
 * The display method displays the two sequences from a .fasta file  on the terminal.
 */

void Sequence::display()
{
	for(char& c : seq1)	//we display the sequence character by character
	{
		cout << c;
	}
	
	cout << endl;
	cout << endl;
	
	for(char& c : seq2)
	{
		cout << c;
	}
	
	cout << endl;
	cout << endl;
};

/*!
 * The searchMotif method looks for a particular motif within the two sequences.
 */
 
vector<int> Sequence::searchMotif(const vector<char>& seq, const string& subStr) const
{
	vector<int> output; 	//the vector of all matching positions
	
	vector<char> subStr_;	//used to convert a string into a table of char
	
	for(const char& c : subStr)	//we convert the substring into a table of characters as well
	{
		subStr_.push_back(c);
	}
		
	for(size_t i(0); i < seq.size() - 6; ++i)	//we check the whole sequence 
			{
					bool matchingCondition(seq[i] == subStr_[0] and seq[i+1]== subStr_[1] and seq[i+2]== subStr_[2]
						and seq[i+3]== subStr_[3] and seq[i+4] == subStr_[4] and seq[i+5] == subStr_[5]
						and seq[i+6] == subStr_[6]); //we know that 7 nucleotides in a row have to match
					
					if(matchingCondition)		//if 7 nucleotides in a row do match then
					{
						output.push_back(i+1);	//we add the position OF THE FIRST MATCHING NUCLEOTIDE to the table 
					}		
			}
	
	return output;
};

/*!
 * The giveComplementarySeq method gives the reverse complementary sequence of any nucleotidic sequences, whether a motif or a genomic sequence.
 */
 
//A method giving the REVERSE complementary sequence from one of the two DNA strand
vector<char> Sequence::giveComplementarySeq(vector<char> seq)
{
		vector<char> complementarySequence;    							//the reverse comp. sequence we get 
		
		for (int position(seq.size()); position !=-1; --position)  	//we start from the end (seq.size()) then go upward in the vector (--position) until top is reached
			{
				if(seq[position] == 'C')								//conversion of nucleotides
					{
						complementarySequence.push_back('G');
					}
				else if(seq[position] == 'G')
					{
						complementarySequence.push_back('C');
					}
				else if(seq[position] == 'A')
					{
						complementarySequence.push_back('T');
					}
				else if(seq[position] == 'T')
					{
						complementarySequence.push_back('A');
					}					
			}
			
			//TESTING NEW METHOD : displaying but not supposed to ?
			for (auto c : complementarySequence)
			{
				cout << c;
			}
			//  cout << complementarySequence.size();    -> we see that the vectors are the same size
			return complementarySequence;
};

//Getters for giveComplementarySequence :
/*!
 * This method makes accessible the first DNA sequence from a .fasta file.
 */
vector <char>  Sequence::getSequence1()
{
	return seq1;
};
/*!
 * This method makes accessible the second DNA sequence from a .fasta file.
 */
vector <char>  Sequence::getSequence2()
{
	return seq2;
};




int main()
{
	Sequence seq_;
	seq_.outputSites("GGATTGG");
	
	//TESTING NEW METHOD GIVECOMPLEMENTARYSEQ
	//cout << "Reverse complementary sequence of first sequence in fasta file :" << endl;
	//seq_.giveComplementarySeq(seq_.getSequence1());
	//END OF TEST
	
	return 0;
}

