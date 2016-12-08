#include <iostream>
#include "sequence.hpp"
#include "utilities.hpp"
#include "Logo.hpp"


using namespace std;

int main()
{
    char answer;
    int nbr(0);
    
    cout << "Hello ! Welcome to the Genom-1 DNA binding site analysis package ! \n" << endl;
    cout << "This program provides a few fonctionalities that will help you analyse and work on genomic sequences. \n" << endl;
    //info to documentation
    
    do {
        cout << " What files would you like to work with : \n " <<endl;
        cout << " - Fasta files (write '1') : this will give you a Position Weight Matrix \n" <<endl;
        cout << " - Genomic fasta files and Bed files (write '2'): this could give you a Position Weight Matrix and a list of motif \n" <<endl;
        cout << " -  Genomic fasta files and BedGraph files (write '3') : this could give you a Position Weight Matrix and a list of motif completed with a bedgraph score \n" <<endl;
        cout << " - Non-genomic fasta file, .mat and a threshold for the score (enter '4') to produce a list of sites \n" << endl;
        cin >> answer;
        ++nbr;
    } while ( (answer != '1') and (answer != '2') and (answer != '3') and (answer != '4') and  (nbr <= 10));
    
    if(nbr > 10) { cout << " A non-sense answer has been entered too many times ! " <<endl; }
    
    else if (answer == '1') {
        cout << "Enter the Fasta files name :" ;
        string fasta;
        cin >> fasta;
        int nbr;
        do {
            cout << "How long are the motifs you want ? (between 1 and 15)";
            cin >> nbr;
        } while ((nbr < 0) and (nbr > 10));
        
        MatrixProtein Protein;
        Sequence sequence;
        vector<string> sequences_;
        
        try {
            sequences_ = sequence.loadSeq(fasta);
        } catch (runtime_error message) {
            cout << message.what();
        }
        
        /* for (size_t i(0); i < sequences_.size(); ++i) {
         vector<char> sequences_char;
         sequences_char = sequence.toVector(sequences_[i]);
         vector<char> sequences_char_inverse;
         sequences_char_inverse = sequence.giveReverseComplementarySeq(sequences_char);
         sequences_.push_back(sequence.toString(sequences_char_inverse));
         
         }*/
        
        
        
        Protein.EMalgorithm(nbr,sequences_, {0,0}, 0);
        
        sequence.loadMatrixOnFile("Matrix_Output", Protein.getmx());
        cout << "The Matrix has been saved on the Output file on the test folder ";
        cout << "Do you want to display the Logo of the PWM ? (type 1)";
        int response;
        cin >> response;
        if (response == 1) {
            Logo logo;
            logo.afficher_logo(Protein.getPatterns(), Protein.getmx());
        }
        
        
        
    }
    
    else if (answer == '2') {
        MatrixProtein Protein;
        Sequence sequence;
        cout << "Enter a Genomic Fasta file \n";
        string Genom;
        cin >> Genom;
        
        
     
        string chr(Genom);
        chr.pop_back();
        chr.pop_back();
        chr.pop_back();
        
        cout << "Enter a Bed file \n";
        string Bed;
        cin >> Bed;
        
        do {
            cout << "How long are the motifs you want ? (between 1 and 15)";
            cin >> nbr;
        } while ((nbr < 0) and (nbr > 10));
        
        
        vector<Coordinate> coordinates;
        coordinates = sequence.readBed(Bed, chr);
        //vector<Coordinate> Coordinate_ (sequence.Convert_BedCToC(coordinates));
        vector<string> fasta_seq (sequence.scanFasta(coordinates, Genom, nbr));
        size_t sizein(fasta_seq.size());
        for (size_t k(0); k < sizein; ++k) {
            vector<char> vect (toVector(fasta_seq[k]));
            fasta_seq.push_back(toString(sequence.giveReverseComplementarySeq(vect)));
        }
        vector<int> debut;
        for (size_t z(0); z < coordinates.size(); ++z) {
            debut.push_back(coordinates[z].start);
        }
        
        Protein.EMalgorithm(nbr,fasta_seq, debut, sizein );
        sequence.loadMatrixOnFile("Matrix_Output", Protein.getmx());
        cout << "The Matrix has been saved on the Matrix_Output file on the test folder \n The List of site has been saved on the Motif_Output on the test folder";
        
        sequence.fillPosDir(Protein, chr);
        sequence.Clean_Motif_Output("Motif_Output");
        
        sequence.loadResultsOnFile("Motif_Output", sequence.getMotifs4Output() );
        cout << "Do you want to display the Logo of this PWM ? (type 1)";
        int response;
        cin >> response;
       if (response == 1) {
            Logo logo;
            logo.afficher_logo(Protein.getPatterns(), Protein.getmx());
        }
    }
    
    
    
    else if (answer == '3') {
        
        MatrixProtein Protein;
        Sequence sequence;
        cout << "Enter a Genomic Fasta file \n";
        string Genom;
        cin >> Genom;
        
        

        string chr(Genom);
        chr.pop_back();
        chr.pop_back();
        chr.pop_back();
        
        cout << "Enter a GraphBed file \n";
        string Bed;
        cin >> Bed;
        
        
        
        
        do {
            cout << "How long are the motifs you want ? (between 1 and 15)";
            cin >> nbr;
        } while ((nbr < 0) and (nbr > 10));
        
        
        vector<Coordinate> Coordinate_ ;
        try {
            Coordinate_ = sequence.readBedGraph(Bed, chr);
        } catch (runtime_error message) {
            cout << message.what();
        }
        vector<Coordinate> relevantCoordinate;
        for (auto& n : Coordinate_) {
            if (n.score >= 2) {
                relevantCoordinate.push_back(n);
            }
        }
        vector<string> fasta_seq (sequence.scanFasta(relevantCoordinate, Genom, nbr));
        size_t sizein(fasta_seq.size());
        for (size_t k(0); k < sizein; ++k) {
            vector<char> vect (toVector(fasta_seq[k]));
            fasta_seq.push_back(toString(sequence.giveReverseComplementarySeq(vect)));
        }
        vector<int> debut;
        for (size_t z(0); z < relevantCoordinate.size(); ++z) {
            debut.push_back(relevantCoordinate[z].start);
        }
        
        Protein.EMalgorithm(nbr,fasta_seq, debut, sizein );
        sequence.loadMatrixOnFile("Matrix_Output", Protein.getmx());
        cout << "The Matrix has been saved on the Matrix_Output file on the test folder \n The List of site has been saved on the Motif_Output on the test folder";
        sequence.fillPosDir(Protein, chr);
        vector<double> Sommes;
        for (size_t x(0); x < sequence.getMotifs4Output().size(); ++x) {
            double tmp (sequence.interval_addition(sequence.getMotifs4Output()[x].pos, Coordinate_));
            Sommes.push_back(tmp);
        }
        sequence.Clean_Motif_Output("Motif_Output");
        sequence.loadResultsOnFile("Motif_Output", sequence.getMotifs4Output(), Sommes );
        cout << "Do you want to display the Logo of this PWM ? (type 1)";
        int response;
        cin >> response;
        if (response == 1) {
            Logo logo;
            logo.afficher_logo(Protein.getPatterns(), Protein.getmx());
        }
    }
    
    else if(answer == '4')
    {
		Sequence seq;
		MatrixProtein matrix_;
		cout << "Enter .mat file name: " << endl;
		string matName;
		cin >> matName;							//need to check what the user enters (to do later)
		matrix_.loadmatrix_fromfile(matName);
		matrix pssm_ = matrix_.getpssm_rel();			// we assume it's a PWM
		
		cout << "Please enter a threshold: " << endl;
		double seuil_;
		cin >> seuil_;
		
		cout << "Please enter .fasta name: " << endl;
		string fastaName_;
		cin >> fastaName_;
		try
		{
		seq.fastaPlusMatrix(fastaName_, pssm_, seuil_);
		seq.loadResultsOnFile("List of sites (fasta + mat).txt");
		}
		catch(runtime_error& e)
		{
			cout << e.what() << endl;
		}
	}
    
    return 0;
}
