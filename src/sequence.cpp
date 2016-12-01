#include <iostream>
#include <fstream>
#include <stdexcept>
#include <list>
#include <algorithm>													//use std::reverse 
#include <sstream>
#include <iomanip>
#include "sequence.hpp"
#include "utilities.hpp"
using namespace std;

//-----CONSTRUCTOR AND DESTRUCTOR-----

Sequence::Sequence() {}

Sequence::~Sequence() {}

//-----

//-----LOCAL FUNCTIONS-----
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
			if(v1[i] != v2[i] and v1[i] != 'N'){return false;}
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
	
	return -1; //added to avoid a warning message
}

//-----SEQUENCE METHODS-----
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
	char base, nucl, extra; 
	string line("");
	vector<char> seq;
	list<char> l;
	vector<PosDir> positions;
	size_t compteur(0);
	size_t compteurSeq(0);
	string chrNb_;
	
	file.open("../test/" + fileName);

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

				for(size_t i(0); i < motif.size(); ++i)
				{
					if(i == 0)
					{	
						l.clear();
						file >> noskipws >> base;
						
						if(isspace(base))
						{
							cout << endl;
							throw runtime_error("Error: no spaces allowed!");
						}
						
						cout << base;
					}
					else
					{
						file >> noskipws >> base;
						
						if(isspace(base))
						{
							cout << endl;
							throw runtime_error("Error: no spaces allowed!");
						}
						cout << base;
					}
					
					l.push_back(base);
				}

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
						positions.push_back({compteur +1,compteurSeq,chrNb_,'+'});
					}
					
					vector<char> cDNA;				//we get the second strand of DNA from .fasta 
					
					for(const char& c : seq)
					{
						cDNA.push_back(giveComplementaryBase(c));
					}
					
					
					if(compare(cDNA, motifComplementary_))
					{
						positions.push_back({compteur +1, compteurSeq,chrNb_, '-'});
					}
				}
			}
		}
	}
	
	file.close();
	for(const PosDir& c : positions) {cout << c.pos << " " << c.dir << " ";}
	return positions;
}

void Sequence::find(MatrixProtein& mtrx, const string& fasta)
{
	ifstream file;								
	char base, nucl, extra; 
	string line("");
	vector<char> seq, secondStrand, reverseMotifTemp_;
	vector<vector<char>> reverseMotifs;	//a table of reverse of ALL motifs in the matrix
	list<char> l;
	vector<PosDir> positions;
	size_t compteur(0);
	size_t compteurSeq(0);
	string chrNb_;
	
	const vector<Pattern> patterns_(mtrx.getPatterns());
	
	for(const Pattern& ptrn : patterns_)
	{
		reverseMotifTemp_.clear(); //we make sure it's clear
		reverseMotifTemp_ = toVector(ptrn.site);	//we take a motif from the patterns vector
		reverse(reverseMotifTemp_.begin(), reverseMotifTemp_.end());	//we get its reverse
		reverseMotifs.push_back(reverseMotifTemp_);		//we store the reverse
	}
	
	file.open("../test/" + fasta);
	
	if(file.fail())
	{
		cout << endl;
		throw runtime_error("The file couldn't be opened!");
	}
	else
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
				file >> line >> ws;
				chrNb_ = chromosomeNb(line);
				
				++compteurSeq;
				compteur = 0;
				
				for(size_t i(0); i < reverseMotifs[0].size(); ++i)
				{
					
					if(i == 0)
					{	
						l.clear();
						file >> noskipws >> base;
						
						if(isspace(base))
						{
							cout << endl;
							throw runtime_error("Error: no spaces allowed!");
						}
						
						cout << base;
					}
					else
					{
						file >> noskipws >> base;
						
						if(isspace(base))
						{
							cout << endl;
							throw runtime_error("Error: no spaces allowed!");
						}
						cout << base;
					}
					
					l.push_back(base);
				}
				
				seq.clear();
				seq = {begin(l), end(l)};
				
				for(const Pattern& ptrn : patterns_)
				{
					if(compare(seq, toVector(ptrn.site)))
					{
						motifs4output.push_back({compteur + 1, compteurSeq, chrNb_, '+', ptrn});
					}
				}
				
				secondStrand.clear();
				
				for(const char& base_ : seq)
				{
					secondStrand.push_back(giveComplementaryBase(base_));
				}
				
				for(size_t i(0); i < reverseMotifs.size(); ++i)
				{
					if(compare(secondStrand, reverseMotifs[i]))
					{
						motifs4output.push_back({compteur +1, compteurSeq, chrNb_, '-', patterns_[i]});
					}
				}
				
				while(file >> noskipws >> nucl)
				{
					cout << nucl;
					file.get(extra);
					if(nucl == '\n' and extra != '>')
					{
						if(!file.eof())
						{
							cout << endl;
							throw runtime_error("Error: new lines are not allowed!");
						}
					}
					else if(nucl == '\n')
					{
						file.putback(extra);
						cout << endl;
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
					
					seq.clear();
					seq = {begin(l), end(l)};
					
					for(const Pattern& ptrn : patterns_)
				{
					if(compare(seq, toVector(ptrn.site)))
					{
						motifs4output.push_back({compteur + 1, compteurSeq, chrNb_, '+', ptrn});
					}
				}
				
				secondStrand.clear();
				
				for(const char& base_ : seq)
				{
					secondStrand.push_back(giveComplementaryBase(base_));
				}
				
				for(size_t i(0); i < reverseMotifs.size(); ++i)
				{
					if(compare(secondStrand, reverseMotifs[i]))
					{
						motifs4output.push_back({compteur +1, compteurSeq, chrNb_, '-', patterns_[i]});
					}
				}
			}
		}
	}
	file.close();
}

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
			
			//TESTING NEW METHOD : displaying but not supposed to !
			for (auto c : complementarySequence)
			{
				cout <<"Generation of reverse complementary sequence :" << endl;
				cout << c;
			}
			
			cout << endl << endl;
			//  cout << complementarySequence.size();    -> we see that the vectors are the same size
			return complementarySequence;
};

vector<Coordinate> Sequence::readBedGraph(const string& fileName) // the function stores data from a file containing a chromosome n°, 2 positions and a score
{
    
    vector<Coordinate> coordinates;
    
    vector<string> temporarySites; // stock in a 1x1 Matrix to calculate the number of data 
    ifstream file;
    file.open(fileName);
    
    
    if (file.fail()) 
		{
			cerr << "This file could not be opened !" <<endl;
		}
    
    else 
    {
        string var;
        
        while (!file.eof()) 
        {
            
            while (!file.eof()) 
            {
                file >> var >> ws;					//we read while we encounter a white space
                temporarySites.push_back(var);
            }
		}
        
        file.close(); 
        
        size_t z(0);									//we now store all informations read in an organized vector (in structures Coordinate)
        while (z < temporarySites.size()) 
        {
            Coordinate site;
            site.chromosome = temporarySites[z];		//we add the chromosome number
            ++z;
            site.start = toInt(temporarySites[z]); 	//then the starting position, etc.
            ++z;
            site.end = toInt(temporarySites[z]);
            ++z;
            site.score = to_Double(temporarySites[z]);
            ++z;
            coordinates.push_back(site);
    
        }
    }
    return coordinates;
    
}

vector<BedCoordinate> Sequence::ReadBed(const string& fileName) // the function stores data from a file containing a chromosome n°and 2 positions
{
    
    vector<BedCoordinate> BedCoordinates;
    
    vector<string> temporarySites; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open(fileName);
    
    
    if (file.fail()) {
        cerr << "This file could not be opened !" <<endl;
    }
    
    else {
        string var;
        
        while (!file.eof()) {
            
            while (!file.eof()) {
                file >> var >> ws;
                temporarySites.push_back(var);
            }
        }
        
        file.close();
        
        
        size_t z(0);     							//we now store all informations read in an organized vector (in strucures BedCoordinate)
        while (z < temporarySites.size()) 
        {
            BedCoordinate c;
            c.chromosome = temporarySites[z];
            ++z;
			c.start = toInt(temporarySites[z]); 	//does not compile here but don't know why ?
            ++z;
			c.end = toInt(temporarySites[z]);
            ++z;
            BedCoordinates.push_back(c);
            
        }
    }
    return BedCoordinates;
    
}

double Sequence::interval_addition(int pos, vector<Coordinate> Coordinates)
{
    
    double a(0);
    size_t i(0);
    
    while (Coordinates[i].end < (pos-50))
    {
        ++i;
    }
    
    size_t j(i);
    
    while (Coordinates[j].end < (pos+50))
    {
        ++j;
    }
    
    for (size_t k(i); k <= j; ++k)
    {
        int start(Coordinates[k].start);
        int end(Coordinates[k].end);
        
        if (start < (pos-50))
        {
            start = (pos-50);
        }
        
        if (end > (pos+50))
        {
            end = (pos+50);
        }
        
        a += (end-start+1)*Coordinates[k].score;
    }
    
    return a;
}




vector<string> Sequence::scanFasta(vector<Coordinate>& coordinates, const string& fileName, int nbrOfSeq)
{
		
	vector<string> listOfMotifs;
	string motif;
	char read;
	string readMotif;
	long long int start, end;
	ifstream file;
	
	file.open("../test/" + fileName);
	if(file.fail())
	{
		cerr << "This file could not be opened !" << endl;
	}
	else 
	{
		file >> read;
		if(read != '>') 
		{
			cout << endl;						
			throw runtime_error("Error: missing header in .fasta!");  //will need to catch
		}
		else
		{ 	
			//trying for all :
			for(int i(0); i <= nbrOfSeq - 1 ; ++i)
			{
				start = (coordinates[i]).start;
				end = (coordinates[i]).end;	
						
				file.seekg(start + 6);		//we go directly to that position, beware the first line is the header !
				long long int sizeOfMotif = (end - start);
				//cout << sizeOfMotif << endl;
				for(long long int i(0); i <= sizeOfMotif; ++i)				 //init at 0 or 1 ?
				{
					file >> read;
					motif += read; 		//we discover the motif string by string
				}
				listOfMotifs.push_back(motif);
			
			}
		}
		}
	//TEST
	file.close();
	return listOfMotifs;
}

vector<vector<char>> Sequence::delete_vectors_too_small(size_t n, vector<vector<char>> target)
{
    
    for(size_t i(0); i <target.size(); ++i) {
        if (target[i].size() < n) {
            target.erase(target.begin() + i);
        }
        
    }
    
    return target;
}

void loadResultsOnFile(const string& fileName, const vector<PosDir>& posdir, double sommeScores)    //fonction that loads on a file (fileName) all the information of a/several sequence(s)
{
    ofstream sortie;
    //sortie.open("../test/" + fileName); //mode écrasement
    sortie.open("../test/" + fileName, ios::out|ios::app);	//mode append (ajout)
    
    if (sortie.fail()) {
        cerr << "coulnd't open the file" << endl;
    } else {
        for(const PosDir& entry : posdir)
        {
			sortie << "seq" << entry.seqNb << " "  << entry.pos << " " << entry.dir << " " ;
			sortie << entry.pattern.site << " "  << entry.pattern.bScore << " "  << sommeScores << endl;
			sortie.close();
		}
    }
}
 
/*int main()
{
	Sequence seq_;
	MatrixProtein mtrx;
	try
	{
		seq_.find(mtrx, "promoters.fasta");
		loadResultsOnFile("List of Sites", seq_.getMotifs4Output() , 0.0);
	}
	catch(const runtime_error& e) {cout << e.what() << endl;}
	
	return 0;
} */




void loadMatrixOnFile(const string& fileName, matrix matrice)   //fonction that loads a matrix on a txt file
{
    ofstream sortie;
    //sortie.open("../test/" + fileName); //mode écrasement
    sortie.open("../test/" + fileName, ios::out|ios::app);	//mode append (ajout)
    
    if (sortie.fail()) {
        cerr << "coulnd't open the file" << endl;
    } else {
        for (size_t i(0); i < matrice.size() ; ++i) {
            for (size_t j(0) ; j < matrice[i].size() ; ++j) {
                sortie << setprecision(5) << setw(5) << matrice[i][j] << " | " ;
            }
            sortie << endl;
        }
        sortie.close();
    }
}
