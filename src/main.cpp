#include <iostream>

#include "sequence.hpp"
#include "Logo.hpp"


using namespace std;

int main() 
{
    cout << "Hello ! Welcome to the Genom-1 DNA binding site analysis package ! \n" << endl;
    cout << "This program provides a few fonctionalities that will help you analyse and work on genomic sequences. \n" << endl;
    //info to documentation
    
    do {
        cout << " What files would you want to work with : \n " <<endl;
        cout << " - Fasta files (write '1') : this will give you a Position Weight Matrix \n" <<endl;
        cout << " - Genomic fasta files and Bed files (write '2'): this could give you a Position Weight Matrix or a list of motif \n" <<endl;
        cout << " -  Genomic fasta files and BedGraph files (write '3') : this could give you a Position Weight Matrix or a list of motif completed with a bedgraph score \n" <<endl;
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
        cout << sequences_[0];
        
        Protein.EMalgorithm(nbr,sequences_);
        
        sequence.loadMatrixOnFile("Output", Protein.getmx());
        cout << "The Matrix has been saved on the Output file on the test folder ";
        
        
        
    }
    
    
    else if (answer == '2') {
        
        // demander une matrice (et la convertir si nécessaire) et la séquence + fonction calculant et affichant les affinity scores
        cout << "Is the Matrix a PWM [type 0] or a PSSM[type 1] \n";
        bool type(0);
        cin >> type;
        cout << "Matrix : Enter the name of a file \n";
        string name;
        cin >> name;
        cout << name;
        MatrixProtein Prot;
        try {
            Prot.loadmatrix_fromfile(name);
        } catch (runtime_error message) {
            cout << message.what();
        }
        
        
    A: cout << "Sequence: enter a sequence \n";
        string sequence_;
        cin >> sequence_;
        while (sequence_.size() > Prot.getpwm_abs().size()) {
            cout << "The sequence is longer than the given Matrix allows. The sequence can be up to " <<Prot.getpwm_abs().size()<<" base long. Please reenter a sequence \n";
            cin >> sequence_;
        }
        vector<char> sequence;
        for (size_t w(0); w < sequence_.size() ; ++w) {
            if ((sequence_[w] != 'A') and (sequence_[w] != 'C') and (sequence_[w] != 'G') and (sequence_[w] != 'T')) {
                cout << "The sequence contains incorrect letters, Please check that all letters are bases and capital.\n";
                goto A;
            }
        }
        for (size_t w(0); w < sequence_.size() ; ++w) {
            sequence.push_back(sequence_[w]);
        }
        if (type)
        {
            cout << "The affinity score for the sequence "<<sequence_<<" is of "<<Prot.get_affinity_score_from_matrix(Prot.getpwm_abs(), sequence)<<"\n";
        } else {
            cout << "The affinity score for the sequence "<<sequence_<<" is of "<<Prot.get_affinity_score_from_matrix(Prot.getpssm_abs(), sequence)<<"\n";
        }
        
    }
    
    
    
    else if (answer == '3') {
        
        char choice;
        int nb(0);
        
        do {
            cout << " Would you like to get a PWM (write '1') or a PSSM (write '2') ? " <<endl;
            cin >> choice;
            ++nb;
        } while ( (choice != '1') and (choice != '2')  and  (nb <= 5));
        
        
        
        if ( choice == '1' or choice == '2')
            
            
        {
            vector<Pattern> patterns_;
            
            int nbr(0);
            
            cout << "How many Sequences do you want to enter ? \n";
            cin >> nbr;
            
            int length(0);
            cout << "How long are the Sequences ? \n";
            cin >> length;
            
            cout << "Enter the list of sequence and for each its binding score : \n" << endl;
            for (int j(0); j < nbr; ++j)
            {
                Pattern motif_;
                cout << "Sequence "<<j+1;
                string x;
                cin >> x;
                while (x.size() != length)
                {
                    
                    cout << "The size of the sequence doesn't correspond to the indicated length, please enter it again \n";
                    cin >> x;
                }
                for (size_t w(0); w < x.size() ; ++w)
                {
                    motif_.site.push_back(x[w]);
                }
                
                cout << "Binding Score"<<j+1;
                cin >> motif_.bScore;
                patterns_.push_back(motif_);
            }
            
            MatrixProtein protein_;
            //protein_.setPatterns(patterns_);getMotifs4Output()
            
            
            
            if (choice == '1')
            {
                protein_.display_PWM_rel();
            }
            
            
            if ( choice == '2')
                
            {
                
            }
        }
        
        
        
        if(nb > 5) { cout << " Abandon ! " <<endl; }
        
        
        // demander fichier .mat et liste de sites + fonction qui crée puis affiche la PWM ou PSSM (faire 2 conditions)
        
        
        
        
    }
    
    
    else if (answer == '4') {
        
        // demander la séquence et les binding scores + fonction pour afficher le logo
        vector<Motif> sequences;
        
        int nbr(0);
        cout << "How many Sequences do you want to enter ? \n";
        cin >> nbr;
        
        int length(0);
        cout << "How long are the Sequences ? \n";
        cin >> length;
        
        cout << "Enter the list of sequence and for each its binding score : \n" << endl;
        for (int j(0); j < nbr; ++j) {
            Motif motif_;
            cout << "Sequence "<<j+1;
            string x;
            cin >> x;
            while (x.size() != length)
            {
                cout << "The size of the sequence doesn't correspond to the indicated length, please enter it again \n";
                cin >> x;
            }
            for (size_t w(0); w < x.size() ; ++w) {
                motif_.lettres.push_back(x[w]);
            }
            
            cout << "Binding Score \n"<<j+1;
            cin >> motif_.binding_score;
            sequences.push_back(motif_);
        }
        
        
        
        Logo Logo_;
        
        Logo_.afficher_logo(sequences);
        
    }   
    return 0;
}
