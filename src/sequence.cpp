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
			if(v1[i] != v2[i] and v1[i] != 'N') {return false;}
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
       or nucl == '-' or nucl == '.' or nucl == 'N' or nucl == 'n' or nucl == '\n')
    {
        if(nucl == 'A' or nucl == 'a') 	 return 'T';
        else if(nucl == 'T' or nucl == 't') return 'A';
        else if(nucl == 'C' or nucl == 'c') return 'G';
        else if(nucl == 'G' or nucl == 'g') return 'C';
        else if(nucl == '-') return '-';
        else if(nucl == '.') return '.';
        else if(nucl == 'N') return 'N';
        else if (nucl == 'n') return 'n';
        else if (nucl == '\n') return '\n';
        
    }
    else
    {
        throw runtime_error("Error: Wrong nucleotide detected in the .fasta!");
    }
    
    return -1; //added to avoid a warning message

}

//-----SEQUENCE METHODS-----



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
	
	file.open("../Resources/" + fileName);

	vector<char> motif_;							//used to convert a string into a table of char
													//used to convert a string into a table of char, the motif we are looking for
													//instead of checking the complementary strand we check for the reverse complement of the motif which should appear 
	for(const char& c : motif)						//we convert the substring into a table of characters as well
	{
		motif_.push_back(c);
	}
		
	vector<char> motifTemp = motif_;
	
	reverse(motifTemp.begin(), motifTemp.end());
		
	vector<char> motifReverse_ = motifTemp;
	
	if(file.fail())								// if it didnt open -> show an error 
	{
		cerr << "The file could not be opened!" << endl;
	}
	else
	{	
		cout << "Scanning " + fileName + "...";
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
				compteur = 0;

				for(size_t i(0); i < motif_.size(); ++i)
				{
					if(i == 0) {l.clear();}
				
					file >> ws >> base;
					
					if(isspace(base))
					{
						cout << endl;
						throw runtime_error("Error: no spaces allowed!");
					}
								
					l.push_back(base);			
				}

				while(file >> noskipws >> nucl)
				{
					file.get(extra);
					if(nucl == '\n' and extra != '>') 
					{
						if(!file.eof()) 
						{
							
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

					if(compare(seq, motif_))
					{
						positions.push_back({compteur +1,compteurSeq,chrNb_,'+', toString(motif_), 0.0});
					}
					
					vector<char> secondStrand;				//we get the second strand of DNA from .fasta 
					
					if(compare(seq, giveReverseComplementarySeq(motif_)))
					{
						positions.push_back({compteur+1, compteurSeq,chrNb_, '-', toString(motif_), 0.0});
					}
				}
			}
		}
	}
	
	file.close();
	cout << endl << "Done reading. Here are the motifs found in " + fileName + ":" << endl;
	for(const PosDir& c : positions) {cout << "pos: " <<c.pos << ", seqNb: " << c.seqNb << ", chrNb: "
			<< c.chrNb << ", dir: " << c.dir << ", bScore: " << c.bindingscore << ", site: " 
			<< c.sequence << endl;}	
	return positions; 
}

void Sequence::fastaPlusMatrix(const string& fastaFile, matrix& pwm, const double& seuil)
{
	ifstream file;								//the file we are going to read
	char base, nucl, extra; 
	string line("");
	vector<char> seq, secondStrand;
	list<char> l;
	size_t compteur(0);
	size_t compteurSeq(0);
	double score;
	string chrNb_;
	
	file.open("../Resources/" + fastaFile);
	
	if(file.fail())								// if it didnt open -> show an error 
	{
		cerr << "The file could not be opened!" << endl;
	}
	else
	{	
		cout << "Scanning " + fastaFile + "..." << endl;
		while(!file.eof())
		{
			file >> nucl;
			if(nucl != '>' and nucl != ' ') 
			{
				cout << endl;
				throw runtime_error("Error: missing header in .fasta!");
			}
			else
			{
				file >> line;
				chrNb_ = chromosomeNb(line);
				compteur = 0;
				++compteurSeq; 

				for(size_t i(0); i < pwm.size(); ++i)
				{
					if(i == 0) {l.clear();}
				
					file >> ws >> base;
						
					l.push_back(base);			
				}
				
				seq = {l.begin(), l.end()};
				
				score = calculateScore(pwm, seq);

				if(score >= seuil)
				{
					motifs4output.push_back({compteur+1,compteurSeq,chrNb_,'+', toString(seq),score});
				}
	
				secondStrand = giveReverseComplementarySeq(seq);
				
				score = calculateScore(pwm, secondStrand);
				
				if(score >= seuil)
				{
					motifs4output.push_back({compteur+1,compteurSeq,chrNb_,'-', toString(secondStrand),score});
				}

				while(!file.eof())
				{
					file >> noskipws >> nucl;
					cout << nucl;
					file.get(extra);
					if(nucl == '\n' and extra != '>') 
					{
						if(extra == '\n') break;
						else
						{
							cout << endl;
							throw runtime_error("Error: new lines inside are not allowed!");
						}
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

					score = calculateScore(pwm, seq);
				
					if(score >= seuil)
					{
						motifs4output.push_back({compteur+1,compteurSeq,chrNb_,'+', toString(seq),score});
					}
						
					secondStrand = giveReverseComplementarySeq(seq);
					
					score = calculateScore(pwm, secondStrand);
					
					if(score >= seuil)
					{
						motifs4output.push_back({compteur+1,compteurSeq,chrNb_,'-', toString(secondStrand),score});
					}
				}
			}
		}
	}
	
	file.close();
	cout << endl << "Done reading " + fastaFile + "." << endl;
}


//A method giving the REVERSE complementary sequence from one of the two DNA strand
vector<char> Sequence::giveReverseComplementarySeq(const vector<char>& seq) const
{

    vector<char> complementarySequence;    							//the reverse comp. sequence we get
    
    for (int position(seq.size() - 1); position > -1; --position)  	//we start from the end (seq.size()) then go upward in the vector (--position) until top is reached
    {
        if((seq[position] == 'C') or (seq[position] == 'c'))								//conversion of nucleotides
        {
            complementarySequence.push_back('G');
        }
        else if((seq[position] == 'G') or (seq[position] == 'g'))
        {
            complementarySequence.push_back('C');
        }
        else if((seq[position] == 'A') or (seq[position] == 'a'))
        {
            complementarySequence.push_back('T');
        }
        else if((seq[position] == 'T') or (seq[position] == 't'))
        {
            complementarySequence.push_back('A');
        }
        
        else if ((seq[position] != 'T') or (seq[position] != 'A') or (seq[position] != 'C') or (seq[position] != 'G') or (seq[position] != 'N'))
        {
			if(seq[position] == ' ') 
			{
				throw runtime_error("Error: No spaces allowed in .fasta files!");
			}
			throw runtime_error("Error: Wrong nucleotide detected in the .fasta!");	
		}
    }
    
    //  cout << complementarySequence.size();    -> we see that the vectors are the same size
    return complementarySequence;
}


vector<Coordinate> Sequence::readBedGraph(const string& fileName, string chrsought) // the function stores data from a file containing a chromosome n°, 2 positions and a score
{
    
    vector<Coordinate> coordinates;
    
    vector<string> temporarySites; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open("../Resources/" + fileName);
    
    

    if (file.fail())
    {
        throw runtime_error ( "This file could not be opened !");
    }
    
    else
    {
        string var;
        
        while (!file.eof())
        {
            Coordinate site;
            for (int f(0); f < 4; ++f) {
                file >> var >> ws;
                if (f == 0 && var == chrsought ) {
                    site.chromosome = chrsought;
                    file >> var >> ws;
                    site.start = toInt(var);
                    file >> var >> ws;
                    site.end = toInt(var);
                    file >> var >> ws;
                    site.score = to_Double(var);
                    coordinates.push_back(site);
                } else {
                    string trash;
                    file >> trash;
                    file >> trash;
                    file >> trash;
                }
                
            }
            
        }
        
        file.close();
        
    }
    return coordinates;
    
}


vector<Coordinate> Sequence::readBed(const string& fileName, string chrsought) // the function stores data from a file containing a chromosome n°and 2 positions
{
    
    vector<Coordinate> Coordinates;

    
    vector<string> temporarySites; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open("../Resources/" + fileName);
    
    
    if (file.fail()) {
        throw runtime_error("This file could not be opened!");
        
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

            Coordinate c;
            if (temporarySites[z] == chrsought) {
                c.chromosome = temporarySites[z];
                ++z;
                c.start = toInt(temporarySites[z]); 
                ++z;
                c.end = toInt(temporarySites[z]);
                ++z;
                Coordinates.push_back(c);
            } else {
                z += 3;
            }
            
            
        }
    }
    return Coordinates;  
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
        
        a += Coordinates[k].score;
    }
    
    return a;
}

vector<string> Sequence::scanFasta(vector<Coordinate>& coordinates, const string& fileName) //version without motif
{

    
    vector<string> listOfMotifs;
    string motif;
    char read;
    string readMotif;
    long long int start, end;
    ifstream file;
    
    file.open("../Resources/" + fileName);
    if(file.fail())
    {
        cerr << "This file could not be opened ! Scan fasta" << endl;
    }
    else
    {
        file >> read;
        if(read != '>')
        {
            
            throw runtime_error("Error: missing header in .fasta! scan");  //will need to catch
        }
        
        else
        {
            //trying for all :
            for(size_t i(0); i < coordinates.size() ; ++i)
            {
                
                start = (coordinates[i]).start;
                end = (coordinates[i]).end;
                motif = "";
                file.seekg(start + 6);		//we go directly to that position, beware the first line is the header !
                long long int sizeOfMotif = (end - start);
                //cout << sizeOfMotif << endl;
                for(long long int j(0); j < sizeOfMotif; ++j)				 //init at 0 or 1 ?
                {
                    int score (0);
                    file >> read;
                    if ((read != 'A') and (read != 'C') and (read != 'T') and (read != 'G') and (read != 'a') and (read != 't') and (read != 'g') and (read != 'c') and (read != 'N') and (read != '>')) {
                        throw runtime_error("Error: the .fasta file contains some unknown characters!");
                    }
                    if ((read == 'A') or (read == 'C') or (read == 'T') or (read == 'G') or (read == 'a') or (read == 't') or (read == 'g') or (read == 'c')) {
                        ++ score; // checks that we don't have a entire line of n
                    }
                    if (score != 0) {
                        motif += read; //we discover the motif string by string
                    }
                    
                }
                if (motif != "") {
                    {
                        listOfMotifs.push_back(motif);
                    }
                    
                }
                
                motif.clear();
                
            }
        }
    }
    //TEST
    file.close();
    makeFasta(listOfMotifs, coordinates);
    return listOfMotifs;
}


vector<string> Sequence::scanFasta(vector<Coordinate>& coordinates, const string& fileName, int size_motif) //version with motif
{

    
    vector<string> listOfMotifs;
    string motif;
    char read;
    string readMotif;
    long long int start, end;
    ifstream file;
    
    file.open("../Resources/" + fileName);
    if(file.fail())
    {
        cerr << "This file could not be opened ! Scan fasta" << endl;
    }
    else
    {
        file >> read;
        if(read != '>')
        {
            
            throw runtime_error("Error: missing header in .fasta! scan");  //will need to catch
        }
        
        else
        {
            //trying for all :
            for(size_t i(0); i < coordinates.size() ; ++i)
            {
                
                start = (coordinates[i]).start;
                end = (coordinates[i]).end;
                motif = "";
                file.seekg(start + 6);		//we go directly to that position, beware the first line is the header !
                long long int sizeOfMotif = (end - start);
                //cout << sizeOfMotif << endl;
                for(long long int j(0); j < sizeOfMotif; ++j)				 //init at 0 or 1 ?
                {
                    int score (0);
                    file >> read;
                    if ((read != 'A') and (read != 'C') and (read != 'T') and (read != 'G') and (read != 'a') and (read != 't') and (read != 'g') and (read != 'c') and (read != 'N') and (read != '>')) {
                        throw runtime_error("Error: the .fasta file contains some unknown characters!");
                    }
                    if ((read == 'A') or (read == 'C') or (read == 'T') or (read == 'G') or (read == 'a') or (read == 't') or (read == 'g') or (read == 'c')) {
                        ++ score; // checks that we don't have a entire line of n
                    }
                    if (score != 0) {
                        motif += read; //we discover the motif string by string
                    }
                    
                }
                if (motif != "") {
                    if ((motif.size() >= size_motif)) {
                        listOfMotifs.push_back(motif);
                    }
                    
                }
                
                motif.clear();
                
            }
        }
    }
    //TEST
    file.close();
    makeFasta(listOfMotifs, coordinates);
    return listOfMotifs;
}

void Sequence::makeFasta(const vector<string>& regions, const vector<Coordinate>& coordinates) const
{
	unsigned int fileNb(1);
	string newFastaName("../Output/ExtractedSequences" + std::to_string(fileNb) + ".fasta"), header;
	ofstream newFasta;
	while(ifstream(newFastaName))
	{
		++fileNb;
		newFastaName = "../Output/ExtractedSequences" + std::to_string(fileNb) + ".fasta";
	}
	
	
	newFasta.open(newFastaName);
	
	if(newFasta.fail()) {throw runtime_error("File could not be created (makeFasta method)!");}
	else
	{
		for(size_t i(0); i < regions.size(); ++i)
		{
			header = ">" + coordinates[i].chromosome + "|" + coordinates[i].chromosome + ":" 
				+ std::to_string(coordinates[i].start) + "-" + std::to_string(coordinates[i].end);
			string addition(regions[i]);
			newFasta << header << endl << addition << endl;
		}
		
	}
}


void Sequence::loadResultsOnFile(const string& fileName, const vector<PosDir>& posdir, vector <double> sommeScores )    //fonction that loads on a file (fileName) all the information of a/several sequence(s)

{
    ofstream sortie;
    sortie.open("../Output/" + fileName, ios::out|ios::app);	//mode append (ajout)
    
    if (sortie.fail()) {
        cerr << "coulnd't open the file" << endl;
    } else {

        int i(0);
        for(const PosDir& entry : posdir)
        {
            sortie << "seq" << i << " " <<entry.chrNb <<" "  << entry.pos << " " << entry.dir << " " ;
            sortie << entry.sequence << " "  << entry.bindingscore << " " <<sommeScores[i] ;
            sortie << "\n";
            ++i;
        }

    }
    sortie.close();
}

void Sequence::loadResultsOnFile(const string& fileName, const vector<PosDir>& posdir )    //fonction that loads on a file (fileName) all the information of a/several sequence(s)

{
    ofstream sortie;
    sortie.open("../Output/" + fileName, ios::out|ios::app);	//mode append (ajout)
    
    if (sortie.fail()) {
        cerr << "coulnd't open the file" << endl;
    } else {
        int i(0);
        for(const PosDir& entry : posdir)
        {
            sortie << "seq" << i << setfill(' ') << setw(4) << entry.chrNb << entry.pos << entry.dir;
            sortie << entry.sequence << entry.bindingscore ;
            sortie << "\n";
            ++i;
        }
    }
    sortie.close();
}

void Sequence::loadResultsOnFile(const string& fileName)    //fonction that loads on a file (fileName) all the information of a/several sequence(s)
{
    ofstream sortie;
    sortie.open("../Output/" + fileName);
    
    if (sortie.fail()) 
    {
        throw runtime_error("Coulnd't open the file");
    } 
    else 
    {
		int const col1(8);
		int const col2(7);
		int const col3(13);
		int const col4(10);

		sortie << "Seq#/chr#  Position Strand     Motif     Score" << endl;
        for(const PosDir& entry : motifs4output)
        {
            sortie << "seq" << entry.seqNb << "/" << entry.chrNb << setw(col1)  << entry.pos << setw(col2) << entry.dir << setw(col3);
            sortie << entry.sequence << setw(col4) << entry.bindingscore;
            sortie << "\n";
        }
    }
    sortie.close();
}


/*int main()
{
	Sequence seq_;
	MatrixProtein mtrx;
	try
	{
		//auto testing = seq_.readBed("BMAL1_sites.bed", "chr7");
		//seq_.scanFasta(testing,"chr7.fa", );
		seq_.motifRecognition("AAATCG", "sites1.fasta");
	}
	catch(const runtime_error& e) {cout << e.what() << endl;}
	
	return 0;
} */

// A TESTER A PLACE DE L'AUTRE, en commentaire pour ne pas faire bug
/*void Sequence::loadMatrixOnFile(const string& fileName, matrix matrice)   //fonction that loads a matrix on a txt file
{
    ofstream sortie;
    sortie.open("../Output/" + fileName); //mode écrasement
    MatrixProtein Matrix;
    
    Matrix.setmx(matrice);
    
    Matrix.matrix_generation();
        
    if (sortie.fail()) {
        cerr << "coulnd't open the file" << endl;
    } else {
		
		sortie << "La matrice PWM relative : "<<endl;
		
        for (size_t i(0); i < pwm_rel.size() ; ++i) {
            for (size_t j(0) ; j < pwm_rel[i].size() ; ++j) {
                sortie << setprecision(5) << setw(5) << pwm_rel[i][j] << " | " ;
            }
            sortie << endl;
        }
        
        sortie << "La matrice PWM absolue : "<<endl;
        
        for (size_t i(0); i < pwm_abs.size() ; ++i) {
            for (size_t j(0) ; j < pwm_abs[i].size() ; ++j) {
                sortie << setprecision(5) << setw(5) << pwm_abs[i][j] << " | " ;
            }
            sortie << endl;
        }
        
        sortie << "La matrice PSSM absolue : "<<endl;
        
        for (size_t i(0); i < pssm_abs.size() ; ++i) {
            for (size_t j(0) ; j < pssm_abs[i].size() ; ++j) {
                sortie << setprecision(5) << setw(5) << pssm_abs[i][j] << " | " ;
            }
            sortie << endl;
        }
        
        sortie << "La matrice PSSM relative : "<<endl;
        
        for (size_t i(0); i < pssm_rel.size() ; ++i) {
            for (size_t j(0) ; j < pssm_rel[i].size() ; ++j) {
                sortie << setprecision(5) << setw(5) << pssm_rel[i][j] << " | " ;
            }
            sortie << endl;
        }
        
        sortie.close();
    }
}*/

void Sequence::loadMatrixOnFile(const string& fileName, matrix matrice)   //fonction that loads a matrix on a txt file
{
    ofstream sortie;
    sortie.open("../Output/" + fileName); //mode écrasement
    
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

void Sequence::Clean_Motif_Output(const string& fileName)
{
    ofstream sortie;
    sortie.open("../test/" + fileName);
    if (sortie.fail()) {
        cerr << "coulnd't open the file" << endl;
    }
}
vector<string> Sequence::loadSeq( string fileName )
{
    vector<string> sequences;
    ifstream file(fileName);    //the file we are reading
    file.open("../Resources/" + fileName);
    
    if(file.fail())        // if it didnt open -> show an error
    {
        throw runtime_error("The Fasta file name is invalid.");
    }
    else
    {
        string seq;
        while(!file.eof())  // EOF is false here
        {
            file >> seq;
            if (seq[0] != '>') {
                sequences.push_back(seq);
            }
        }
        
        file.close();
    }
    return sequences;
}


void Sequence::fillPosDir(MatrixProtein& MX, string chrn)
{
    vector<Pattern> patterns_ = MX.getPatterns();
    for(size_t i(0); i < patterns_.size(); ++i) 
    {
		motifs4output.push_back({static_cast<size_t>(patterns_[i].pos),2, chrn, patterns_[i].dir, patterns_[i].site, patterns_[i].bScore}); 
    }
    
}

vector<Coordinate> Sequence::Convert_BedCToC(vector<BedCoordinate> BedC)
{
    vector<Coordinate> C;
    for (size_t i(0); i < BedC.size(); ++i) {
        Coordinate X( {BedC[i].chromosome, BedC[i].start, BedC[i].end, 0});
        C.push_back(X);
    }
    return C;
}

vector<PosDir> Sequence::getMotifs4Output() const
{
    return motifs4output;
}

