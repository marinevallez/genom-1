//
//  main1.cpp
//  GENOM_1
//
#include <iostream>

using namespace std;

int main() {
    
    char answer;
    int nbr(0);
    
    do {
        cout << " What would you like to do : " <<endl;
        cout << " - Outputting the list of sites (write '1')" <<endl;
        cout << " - Getting the affinity scores (write '2')" <<endl;
        cout << " - Outputting the motif description : PWM or PSSM (write '3')" <<endl;
        cout << " - Outputting the motif logo (write '4')" <<endl;
        cin >> answer;
        ++nbr;
    } while ( (answer != '1') and (answer != '2') and (answer != '3') and (answer != '4') and  (nbr <= 10));
    
    if(nbr > 10) { cout << " Abandon ! " <<endl; }
    
    else if (answer == '1') {
            
        
            // demander le .fasta et les motifs + fonctions pour trouver puis afficher les sites
            
    }
    
     else if (answer == '2') {
         
            // demander une matrice (et la convertir si nécessaire) et la séquence + fonction calculant et affichant les affinity scores
            
     }
            
            
      else if (answer == '3') {
          
          char choice;
          int nb(0);
          
          do {
              cout << " Would you like to get a PWM (write '1') or a PSSM (write '2') ? " <<endl;
              cin >> choice;
              ++nb;
          } while ( (choice != '1') and (choice != '2')  and  (nb <= 5));
          
          if(nb > 5) { cout << " Abandon ! " <<endl; }
          
          // demander fichier .mat et liste de sites + fonction qui crée puis affiche la PWM ou PSSM (faire 2 conditions)
          
          
            
     }
            
            
     else if (answer == '4') {
         
          // demander la séquence et les binding scores + fonction pour afficher le logo
         
      }
            
            
    
    
    return 0;
    
}
