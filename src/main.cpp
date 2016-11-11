//
//  main1.cpp
//  GENOM_1
//
#include <iostream>
#include "Sequence.hpp"
#include "Logo.hpp"
#include "matrix_reader.cpp"
#include "get_afinity_score_from_martix.cpp"

using namespace std;

int main() {
    
    char answer;
    int nbr(0);
    
    do {
        cout << " What would you like to do : " <<endl;
        cout << " - Outputting the list of sites (write '1')" <<endl;
        cout << " - Getting the affinity scores (write '2')" <<endl;
        cout << " - Outputting the motif description : PWM or PSSM (write '3')" <<endl;
        cout << " - Outputting the motif logo (write '4') \n" <<endl;
        cin >> answer;
        ++nbr;
    } while ( (answer != '1') and (answer != '2') and (answer != '3') and (answer != '4') and  (nbr <= 10));
    
    if(nbr > 10) { cout << " Abort ! " <<endl; }
    
    else if (answer == '1') {
        
            // demander le .fasta et les motifs + fonctions pour trouver puis afficher les sites 
			 Sequence seq;
			 Protein protein;
             string fastaFile;
             string matFile;
             int nbr(0);
            do { 
				cout << "Please type the file name (.fasta format) which you would like to use." << endl;
				cin >> fastaFile;
				
				cout << "Are you interested in finding a particular motif (write '1') or just displaying the list of possibles motifs from a PWM matrix (write '2') ?"
				char task;
				cin >> task;
				
				if(task == '1')
				{
					vector<char> motif;
					do {
					
						cout << "Please type the motif of interest (must be 7 bases long) composed of the A, C, T, G nucleotides : " << endl;
						cin >> motif;
					
						} while (motif.size() != 7);
					
					cout << "The " << motif << " motif will be searched in the file's sequences." << endl;
					
					vector<PosDir> info = seq.motifRecognition(motif, fastaFile);
					seq.outputSites(info);
 				
				} else if (task == '2') {
				
					cout >> "Please type the file name (.mat format) which you would like to use to obtain several motifs." << endl;
					cin >> matFile;
					// HERE WRITE METHOD THAT FINDS THE MOTIFS
					//
					//
					//
					
					
					do {
						cout >> "Please enter an affinity score as threshold, above which motifs with a higher score will be considered (must be lower than -1)." << endl;
						cout >> "By default, the programm offers a threshold." << endl; //maybe show the default threshold ?
						double threshold;
						cin << threshold;
					} while ((threshold < 1)); //check definition of score HERE !
					
					//motifRecognition is overloaded according to the subtask
					
					vector<PosDir> info = seq.motifRecognition(protein.getPatterns(), threshold);
					seq.outputSites(info);
				}
			
				++nbr
			} while ((task != '1') and (task != '2') and (nbr <= 10));
			if (nbr > 10) { cout << " Abort ! " << endl; }			
            
            
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
         vector<vector<double>> Matrix;
         try {
             Matrix = loadmatrix(name);
         } catch (runtime_error message) {
             cout << message.what();
         }
         
         
        A: cout << "Sequence: enter a sequence \n";
         string sequence_;
         cin >> sequence_;
         while (sequence_.size() > Matrix.size()) {
             cout << "The sequence is longer than the given Matrix allows. The sequence can be up to " <<Matrix.size()<<" base long. Please reenter a sequence \n";
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
             cout << "The affinity score for the sequence "<<sequence_<<" is of "<<get_afinity_score_from_matrix(Matrix, sequence, 1)<<"\n";
         } else {
             cout << "The affinity score for the sequence "<<sequence_<<" is of "<<get_afinity_score_from_matrix(Matrix, sequence)<<"\n";
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
           vector<Motif> patterns_;
         
				 int nbr(0);

				 cout << "How many Sequences do you want to enter ? \n";
				 cin >> nbr;
				 
				 int length(0);
				 cout << "How long are the Sequences ? \n";
				 cin >> length;
				 
				 cout << "Enter the list of sequence and for each its binding score : \n" << endl;
				 for (int j(0); j < nbr; ++j)
				 {
					 Motif motif_;
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
						 motif_.listOfsites.push_back(x[w]);
					 }
					 
					 cout << "Binding Score"<<j+1;
					 cin >> motif_.bScore;
					 patterns_.push_back(motif_);
				 }
				 
				 Protein protein_;
				 protein_.patterns(patterns_);
				 protein_.loadmatrix_fromscore(); 
				 
				 

				 if (choice == '1')
				  {
					protein_.display_PWM(); 
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
