#include <iostream>
#include <fstream>
#include <stdexcept>
#include <list>
#include <algorithm>		//we want to use std::reverse method from this library 
#include <sstream>
#include "sequence.hpp"
using namespace std;


Sequence::Sequence() {}

Sequence::~Sequence() {}

string chromosomeNb(const string& str)  //a function that extracts the chromosome number from .fasta
{
	string chr;
	 
	if(str[4] == '|')
	{
		chr = str.substr(0,3) + str[3];
	}
	else
	{
		chr = str.substr(0,3) + str.substr(3,2);
	}
	
	return chr;	
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

char giveComplementaryBase(const char& nucl)				//it's a function because it has nothing to do with the sequence class, it's a function to convert letters
{
	if(nucl == 'A' or nucl == 'T' or nucl == 'C' or nucl == 'G'
		or nucl == 'a' or nucl == 't' or nucl == 'c' or nucl == 'g'
			or nucl == '-' or nucl == '.' or nucl == 'N')
	{
		if(nucl == 'A' or nucl == 'a') 	 return 'T';
		else if(nucl == 'T' or nucl == 't') return 'A';
		else if(nucl == 'C' or nucl == 'c') return 'G';
		else if(nucl == 'G' or nucl == 'g') return 'C';
		else if(nucl == '-') return '-';
		else if(nucl == '.') return '.';
		else if(nucl == 'N') return 'N';
	}
	else 
	{
		cerr << endl;
		throw runtime_error("Error: Wrong nucleotide detected in the .fasta: ");
	}
}

void Sequence::outputSites(const vector<PosDir>& info) const	//a method that outputs a file with seq nb and a position where the motif was found
{
	ofstream file;
	file.open("../test/ListofSites.txt");	//we create a file where we will output results
	if(file.fail())
	{
		cerr << "File wasn't found!" << endl;
	}
	else
	{
		for(const PosDir& c : info)
		{
			file << "seq" << c.seqNb << "/" << c.chrNb << " " << c.pos << " " << c.dir << " " << endl;
		}
		
		file.close();		//we close the file
	}
}



vector<PosDir> Sequence::motifRecognition(const string& motif, const string& fileName) const 
{
	ifstream file;								//the file we are going to read
	char c1,c2,c3,c4,c5,c6, c7, nucl, extra; 
	string line("");
	vector<char> seq;
	list<char> l;
	vector<PosDir> positions;
	size_t compteur(0);
	size_t compteurSeq(0);
	string chrNb_;
	
	//file.open("../test/promoters.fasta"); 			//since out text files are in the test folder, we need to include a path to it
	
	file.open("../test/" + fileName); //how to include the right path to get it ? have to do that ?

	vector<char> motif_;							//used to convert a string into a table of char
													//used to convert a string into a table of char, the motif we are looking for
													//instead of checking the complementary strand we check for the reverse complement of the motif which should appear 
	for(const char& c : motif)						//we convert the substring into a table of characters as well
	{
		motif_.push_back(c);
	}
		
	vector<char> motifTemp = motif_;
	
	reverse(motifTemp.begin(), motifTemp.end());
		
	vector<char> motifComplementary_ = motifTemp;
	
	if(file.fail())								// if it didnt open -> show an error 
	{
		cerr << "The file could not be opened!" << endl;
	}
	else
	{	
		while(!file.eof())
		{
			file >> nucl;
			if(nucl != '>') 
			{
				cout << endl;
				throw runtime_error("Error: missing header in .fasta!");
			}
			else
			{
				file >> line;
				chrNb_ = chromosomeNb(line);
				
				++compteurSeq; 

				file >> c1 >> noskipws >> c2 >> noskipws >> c3 >> noskipws >>
					noskipws >> c4 >> noskipws>> c5 >> c6 >> noskipws >> c7;
				
				cout << c1 << c2 << c3 << c4 << c5 << c6 << c7;
				
				l = {c1,c2,c3,c4,c5,c6,c7};

				while(file >> noskipws >> nucl)
				{
					cout << nucl;
					file.get(extra);
					if(nucl == '\n' and extra != '>') 
					{
						if(!file.eof()) 
						{
							cout << endl;
							throw runtime_error("Error: new lines inside are not allowed!");
						}
						if(file.eof()) {file.close();}
					}
					else if(nucl == '\n') 
					{
						file.putback(extra);
						break;
					}
					else if(isspace(nucl)) 
					{
						cout << endl;
						throw runtime_error("Error: no spaces allowed!");
					}
					else
					{
						file.putback(extra);
						l.pop_front();
						l.push_back(nucl);
					}
					
					++compteur;
					seq = {begin(l), end(l)};
					
					/*for(auto& c : seq) {cout << c;}
					cout << " ";*/

					if(compare(seq, motif_))
					{
						positions.push_back({compteur,compteurSeq,chrNb_,'+'});
					}
					
					vector<char> cDNA;				//we get the second strand of DNA from .fasta 
					
					for(const char& c : seq)
					{
						cDNA.push_back(giveComplementaryBase(c));
					}
					
					
					if(compare(cDNA, motifComplementary_))
					{
						positions.push_back({400 - compteur - 5, compteurSeq,chrNb_, '-'});
					}
				}
			}
		}
	}
	
	file.close();
	for(const PosDir& c : positions) {cout << c.pos << " " << c.dir << " ";}
	return positions;
}


//The motifRecognition method is overloaded
/*vector<PosDir> Sequence::motifRecognition(const MatrixProtein& protein, const string& fileName) const
{
    ifstream file;										//the file we are going to read
    char c1,c2,c3,c4,c5,c6, nucl;
    string line("");
    vector<char> seq;
    list<char> l;
    vector<PosDir> positions;
    size_t compteur(0);
    size_t compteurSeq(0);
    string chrNb_;
    
    //file.open("../test/promoters.fasta"); 			//since out text files are in the test folder, we need to include a path to it
    
    file.open(fileName); 								//how to include the right path to get it ? have to do that ?
    
    vector<char> motif_;								//used to convert a string into a table of char
    //used to convert a string into a table of char, the motif we are looking for
    //instead of checking the complementary strand we check for the reverse complement of the motif which should appear
    
    if(file.fail())										// if it didnt open -> show an error
    {
        cerr << "The file could not be opened!" << endl;
    }
    
    
    while(!file.eof())
    {
        while(file >> nucl)
        {
            if(nucl == '>')
            {
                file >> line;
                chrNb_ = chromosomeNb(line);
                break;
            }
            else
            {
                if(compteur == 394 or compteur == 0) 			//reinitialisaiton of the compteur for a new sequence of the fasta file
                {
                    compteur = 0;
                    ++compteurSeq;
                    file >> c1 >> c2 >> c3 >> c4 >> c5 >> c6;
                    l = {nucl,c1,c2,c3,c4,c5,c6};
                    
                }
                else
                {
                    l.pop_front();
                    l.push_back(nucl);
                }
                
                ++compteur;
                seq = {begin(l), end(l)};
                
                for (size_t motifNb(0); motifNb < protein.getPatterns().size(); ++motifNb)
                {
                    
                    vector<char> motif(protein.getPatterns()[motifNb].site);
                    for(const char& c : motif)									//we convert the substring into a table of characters as well
                    {
                        motif_.push_back(c);
                    }
                    
                    if(compare(seq, motif_))
                    {
                        positions.push_back({compteur,compteurSeq,chrNb_,'+', protein.getPatterns()[motifNb]});
                    }
                }
                
                for (size_t motifNb(0); motifNb < protein.getPatterns().size(); ++motifNb)
                {
                    vector<char> motif(protein.getPatterns()[motifNb].site);
                    for(const char& c : motif)									//we convert the substring into a table of characters as well
                    {
                        motif_.push_back(c);
                    }
                    vector<char> motifTemp = motif_;
                    reverse(motifTemp.begin(), motifTemp.end());	
                    vector<char> motifComplementary_ = motifTemp;
                    
                    vector<char> complementaryDNA;				//we get the second strand of DNA from .fasta 				
                    for(const char& c : seq)
                    {
                        complementaryDNA.push_back(giveComplementaryBase(c));
                    }	
                    
                    if(compare(complementaryDNA, motifComplementary_))
                    {
                        positions.push_back({400 - compteur - 5, compteurSeq,chrNb_, '-', protein.getPatterns()[motifNb]});
                    }
                }
            }
        }
    }
    
    file.close();
    for(const PosDir& c : positions) {cout << c.pos << " " << c.dir << " ";}
    return positions;
}*/



//A method giving the REVERSE complementary sequence from one of the two DNA strand
vector<char> Sequence::giveReverseComplementarySeq(const vector<char>& seq) const
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

vector<Coordinate> Sequence::readBedGraph(const string& fileName) // the function stores data from a file containing a chromosome n°, 2 positions and a score
{
    
    string fichier(fileName);
    
    vector<Coordinate> Coordinates;
    
    vector<string> temp; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open(fichier);
    
    
    if (file.fail()) {
        cerr << "This file can't be openned" <<endl;
    }
    
    else {
        string var;
        
        while (!file.eof()) {
            
            while (!file.eof()) {
                file >> var >> ws;
                temp.push_back(var);
            }
        }
        
        file.close();
        
        
        unsigned int z(0);
        
        while (z < temp.size()) {
            Coordinate c;
            c.chromosome = temp[z];
            ++z;
            c.start = To_int(temp[z]);
            ++z;
            c.end = To_int(temp[z]);
            ++z;
            c.score = To_double(temp[z]);
            ++z;
            Coordinates.push_back(c);
            
        }
    }
    return Coordinates;
    
};


//Conversions

double To_double(const string& str) // allows a convertion from string to double
{
    istringstream stream(str);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
};

double To_int(const string& str) // allows a convertion from string to int
{
    istringstream stream(str);
    int a;
    if (!(stream >> a))
        return 0;
    return a;
};

int main()
{
	Sequence seq_;
	try
	{
		vector<PosDir> info = seq_.motifRecognition("ACTGTCA", "SeqFail2.fasta");
		seq_.outputSites(info);
	}
	catch(const runtime_error& e) {cout << e.what() << endl;}
	
	return 0;
}

/*ACTGTCA

5' TGACAGT 3'
3' ACTGTCA 5'

GAGCTC
CTCGAG*/

