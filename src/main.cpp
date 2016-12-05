#include <iostream>

#include "sequence.hpp"
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
        cout << " What files would you want to work with : \n " <<endl;
        cout << " - Fasta files (write '1') : this will give you a Position Weight Matrix \n" <<endl;
        cout << " - Genomic fasta files and Bed files (write '2'): this could give you a Position Weight Matrix and a list of motif \n" <<endl;
        cout << " -  Genomic fasta files and BedGraph files (write '3') : this could give you a Position Weight Matrix and a list of motif completed with a bedgraph score \n" <<endl;
        cin >> answer;
        ++nbr;
    } while ( (answer != '1') and (answer != '2') and (answer != '3') and  (nbr <= 10));
    
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
        
        //COUGAR
        //Genom = "/Users/oriane/Desktop/chr11.fa";
        string chr("chr11.fa");
        chr.pop_back();
        chr.pop_back();
        chr.pop_back();
        
        cout << "Enter a Bed file \n";
        string Bed;
        cin >> Bed;
        //COUGAR
        //Bed = "/Users/oriane/Desktop/BMAL1_sites.bed";
        
        do {
            cout << "How long are the motifs you want ? (between 1 and 15)";
            cin >> nbr;
        } while ((nbr < 0) and (nbr > 10));
        
        
        vector<BedCoordinate> Bedcoordinates;
        Bedcoordinates = sequence.ReadBed(Bed, chr);
        vector<Coordinate> Coordinate_ (sequence.Convert_BedCToC(Bedcoordinates));
        vector<string> fasta_seq (sequence.scanFasta(Coordinate_, Genom, nbr));
        size_t sizein(fasta_seq.size());
        for (int k(0); k < sizein; ++k) {
            vector<char> vect (Protein.toVector(fasta_seq[k]));
            fasta_seq.push_back(sequence.toString(sequence.giveReverseComplementarySeq(vect)));
        }
        vector<int> debut;
        for (int z(0); z < Coordinate_.size(); ++z) {
            debut.push_back(Coordinate_[z].start);
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
        
        //COUGAR
        //Genom = "/Users/oriane/Desktop/chr11.fa";
        string chr("chr11.fa");
        chr.pop_back();
        chr.pop_back();
        chr.pop_back();
        
        cout << "Enter a GraphBed file \n";
        string Bed;
        cin >> Bed;
        //COUGAR
        //Bed = "/Users/oriane/Desktop/BMAL1_ZT06_selection.bedgraph";
        
        do {
            cout << "How long are the motifs you want ? (between 1 and 15)";
            cin >> nbr;
        } while ((nbr < 0) and (nbr > 10));
        
        
        vector<Coordinate> Coordinate_ (sequence.readBedGraph(Bed, chr));
        vector<Coordinate> relevantCoordinate;
        for (auto& n : Coordinate_) {
            if (n.score >= 2) {
                relevantCoordinate.push_back(n);
            }
        }
        vector<string> fasta_seq (sequence.scanFasta(relevantCoordinate, Genom, nbr));
        size_t sizein(fasta_seq.size());
        for (int k(0); k < sizein; ++k) {
            vector<char> vect (Protein.toVector(fasta_seq[k]));
            fasta_seq.push_back(sequence.toString(sequence.giveReverseComplementarySeq(vect)));
        }
        vector<int> debut;
        for (int z(0); z < relevantCoordinate.size(); ++z) {
            debut.push_back(relevantCoordinate[z].start);
        }
        
        Protein.EMalgorithm(nbr,fasta_seq, debut, sizein );
        sequence.loadMatrixOnFile("Matrix_Output", Protein.getmx());
        cout << "The Matrix has been saved on the Matrix_Output file on the test folder \n The List of site has been saved on the Motif_Output on the test folder";
        sequence.fillPosDir(Protein, chr);
        vector<double> Sommes;
        for (int x(0); x < sequence.getMotifs4Output().size(); ++x) {
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
    
    return 0;
}
