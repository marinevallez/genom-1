#include <iostream>
#include <fstream>
#include <stdexcept>
#include <list>
#include "sequence.hpp"
using namespace std;

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

Sequence::Sequence() {}

Sequence::~Sequence() {}

void Sequence::outputSites(const string& motif, const string& fileName)	//a method that outputs a file with seq nb and a position where the motif was found
{
	ofstream file;
	file.open("../test/sites.txt");	//we create a file where we will output results
	if(file.fail())
	{
		cerr << "File wasn't found!" << endl;
	}
	else
	{
		vector<vector<char>> seq = quickRead(fileName);	//we get read both sequences from .fasta
		vector<vector<char>> seqComp = {giveComplementarySeq(seq[0]), giveComplementarySeq(seq[1])};	//then we get their compliments
		for(size_t i(0); i < 2; ++i)
		{
			vector<int> pos = searchMotif(seq[i], motif);	//we find where the motif appears in the provided sequences
			if(!pos.empty())	//if it appears then we save it to file
			{
				for(auto& p : pos)
				{
					file << "The motif was found on seq" << i+1 << " starts at: " << p << " ends at: " << p+7 <<  " + " << motif << endl;
				}
			} 
			pos = searchMotif(seqComp[i], motif);	//then we check the compliment of those sequences
			if(!pos.empty())						//if it appears there then we also save it
			{
				for(auto& p : pos)
				{
					file << "The motif was found on seq" << i+1 << " stars at: " << 400 - p << " ends at: " << 400 - p - 7 << " - "<< motif << endl;
				}
			}
		}
		file.close();		//we close the file
	}
}

bool compare(const vector<char>& v1, const vector<char>& v2)
{
	if(v1.size() == v2.size())
	{
		for(size_t i(0); i < v1.size(); ++i)
		{
			if(v1[i] != v2[i]) {return false;}
		}
		
		return true;
	}
	else
	{
		cerr << "Can't compare these vectors!" << endl;		//maybe throw an error
		return false;
	}
}

vector<PosDir> Sequence::motifRecognition(const string& motif) const 
{
	ifstream file;								//the file we are going to read
	char c1,c2,c3,c4,c5,c6,c7, nucl; 
	string line("");
	vector<PosDir> positions;
	size_t compteur(1);
	
	file.open("../test/promoters.fasta"); 			//since out text files are in the test folder, we need to include a path to it
	
	vector<char> motif_;							//used to convert a string into a table of char
	
	for(const char& c : motif)	//we convert the substring into a table of characters as well
	{
		motif_.push_back(c);
	}
	
	if(file.fail())								// if it didnt open -> show an error 
	{
		cerr << "File could not be opened!" << endl;
	}
	else
	{
		file >> line >> ws;
		file >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7;
		
		list<char> l = {c1,c2,c3,c4,c5,c6,c7};
		vector<char> seq = {l.begin(), l.end()};
		
		bool matchingCondition(compare(seq, motif_));
		
		if(matchingCondition)
		{
			positions.push_back({compteur,'+'});
		}
		
		while(file >> nucl)
		{
			if(nucl == '>') 
			{
				break;
			}
			else
			{
				++compteur;
				
				l.pop_front();
				l.push_back(nucl);
				
				seq = {begin(l), end(l)};
				matchingCondition = compare(seq, motif_);
				
				if(matchingCondition)
				{
					positions.push_back({compteur, '+'});
				}
			}
		}
	}
	
	file.close();
	for(const PosDir& c : positions) {cout << c.pos << " ";}
	return positions;
}



char revComp(const char& nucl) 
{
	if(nucl == 'A' or nucl == 'T' or nucl == 'C' or nucl == 'G')
	{
		if(nucl == 'A') return 'T';
		else if(nucl == 'T') return 'A';
		else if(nucl == 'C') return 'G';
		else if(nucl == 'G') return 'C';
	}
	else 
	{
		throw runtime_error("Error! Wrong nucleotide entered!");
	}
}

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
		//cout << c;
	}
	
	output.push_back(v1);
	
	//cout << endl << endl;
	
	for(char& c : seq2_)
	{
		v2.push_back(c);
		//cout << c;
	}
	
	//cout << endl << endl;
	
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
			
			cout << endl << endl;
			//  cout << complementarySequence.size();    -> we see that the vectors are the same size
			return complementarySequence;
};

//Getters for giveComplementarySequence :
vector <char>  Sequence::getSequence1()
{
	return seq1;
};

vector <char>  Sequence::getSequence2()
{
	return seq2;
};




int main()
{
	Sequence seq_;
	//seq_.outputSites("GGAGAGT", "promoters.fasta");
	seq_.motifRecognition("TGACTAT");
	//TESTING NEW METHOD GIVECOMPLEMENTARYSEQ
	//cout << "Reverse complementary sequence of first sequence in fasta file :" << endl;
	//seq_.giveComplementarySeq(seq_.getSequence1());
	//END OF TEST
	
	return 0;
}

