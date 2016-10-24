

#include <SFML/Graphics.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Letters.hpp"
#include "Axes.hpp"

sf::Texture texture1;
sf::Texture texture2;
sf::Texture texture3;
sf::Texture texture4;
struct Motif
{
    std::vector<char> lettres;
    double binding_score;
};

void test(double& L)
{
    if (L == 0) {
        L = pow(10, -100);
    }
}

vector<vector<double>> Processing(vector<Motif> motif_)
{
    
    vector<vector<double>> Matrice;
    double compteurA(0);
    double tailleA(0);

    double compteurT(0);
    double tailleT(0);
  
    double compteurG(0);
    double tailleG(0);

    double compteurC(0);
    double tailleC(0);

    double nombreTot(0);
    for (size_t j(0) ; j < 7 ; ++j)
    {
        for ( size_t i(0) ; i < motif_.size() ; ++i)
        {
            
            nombreTot = nombreTot + motif_[i].binding_score ;
            if ( motif_[i].lettres[j] == 'A')
            {
                compteurA = compteurA + motif_[i].binding_score*log2(motif_[i].binding_score) ;
                tailleA += (2 - compteurA)*motif_[i].binding_score;

            }
            if ( motif_[i].lettres[j] == 'T')
            {
                compteurT = compteurT + motif_[i].binding_score*log2(motif_[i].binding_score) ;
                tailleT += (2 - compteurT)*motif_[i].binding_score;

            }
            if ( motif_[i].lettres[j] == 'G')
            {
                compteurG = compteurG + motif_[i].binding_score*log2(motif_[i].binding_score) ;
                tailleG += (2 - compteurG)*motif_[i].binding_score;

            }
            if ( motif_[i].lettres[j] == 'C')
            {
                compteurC = compteurC + motif_[i].binding_score*log2(motif_[i].binding_score) ;
                tailleC += (2 - compteurC)*motif_[i].binding_score;

            }
            
        }
        test(tailleA);
        test(tailleT);
        test(tailleG);
        test(tailleC);
        
        Matrice.push_back({tailleA ,tailleC, tailleG ,tailleT});
        compteurA = 0;
        tailleA = 0;
        compteurT = 0;
      tailleT = 0;
        compteurG = 0;
      tailleG = 0;
        compteurC = 0;
        tailleC = 0;
        nombreTot = 0;

}
    return Matrice;
}

double To_double( const string& string ) // allows a convertion from string to double
{
    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}





sf::RectangleShape set_text( double taille, int lettre, double emplacement, int largeur)
{
    sf::Texture* ptexture;


    switch (lettre) {
        case 0 :
        {
            
            texture1.loadFromFile("../res/A.png");
            ptexture = &texture1;
            break;
        }
        case 1 :
        {
            
            texture2.loadFromFile("../res/C.png");
            ptexture = &texture2;
            break;
        }
        case 2 :
        {
            
            texture3.loadFromFile("../res/G.png");
            ptexture = &texture3;
            break;

        }
        case 3 :
        {
            
            texture4.loadFromFile("../res/T.png");
            ptexture = &texture4;
            break;
        }
            
            
        default:
            break;
    }
    sf::RectangleShape rectangle;
    rectangle.setSize(sf::Vector2f(100,pow(6,taille)*10));
    rectangle.setFillColor(sf::Color(250,250,250,250));
    rectangle.setPosition(largeur, emplacement);
    rectangle.setTexture(ptexture);
   
    
    return rectangle;
}
double size( int ltr, Letters tentative)
{
    return pow(6,tentative.get_lettre_size()[ltr])*10;
}

int afficher_logo(vector<vector<double>> Matrice)
{
    sf::RenderWindow window(sf::VideoMode(1400, 500), "PWM logo");
    
    //Chargement fonte
    sf::Font font;
    if (!font.loadFromFile("../res/arial.ttf")) return EXIT_FAILURE;
     window.clear(sf::Color(250,250,250));
    Axes Axes_;
    Axes_.draw(window);
   
    Letters tentative;
    int largeur(95);
    
    for (size_t j(0); j < Matrice.size(); ++j) {
    tentative.calculate_size(Matrice[j]);
    double nbr1 (0.0);
    int lt1 (0);
    double nbr2(0.0);
    int lt2 (0);
    double nbr3(0.0);
    int lt3 (0);
    double nbr4(0.0);
    int lt4 (0);
    for (size_t i(0); i < 4; ++i) {
        if (Matrice[j][i] > nbr1) {
            nbr4 = nbr3;
            lt4 = lt3;
            nbr3 = nbr2;
            lt3 = lt2;
            nbr2 = nbr1;
            lt2 = lt1;
            nbr1 = Matrice[j][i];
            lt1 = i;
        } else if (Matrice[j][i] > nbr2) {
            nbr4 = nbr3;
            lt4 = lt3;
            nbr3 = nbr2;
            lt3 = lt2;
            nbr2 = Matrice[j][i];
            lt2 = i;
        } else if (Matrice[j][i] > nbr3)
        {
            nbr4 = nbr3;
            lt4 = lt3;
            nbr3 = Matrice[j][i];
            lt3 = i;
        } else  {
            nbr4 = Matrice[j][i];
            lt4 = i;
        }
    }
    
    int position_axe_x(348);
    
    sf::RectangleShape rectangle1(set_text(tentative.get_lettre_size()[lt1], lt1, position_axe_x - size(lt4, tentative) - size(lt3, tentative) - size(lt2, tentative) - size(lt1, tentative), largeur));
    
    window.draw(rectangle1);
    
    sf::RectangleShape rectangle2(set_text(tentative.get_lettre_size()[lt2], lt2, position_axe_x - size(lt4, tentative) - size(lt3, tentative) - size(lt2, tentative), largeur));
    window.draw(rectangle2);
    sf::RectangleShape rectangle3(set_text( tentative.get_lettre_size()[lt3], lt3, position_axe_x - size(lt4, tentative) - size(lt3, tentative), largeur));
    window.draw(rectangle3);
    sf::RectangleShape rectangle4(set_text( tentative.get_lettre_size()[lt4], lt4,position_axe_x - size(lt4, tentative) , largeur));
    window.draw(rectangle4);
        
        largeur+=75;
    }

    window.display();
    while (window.isOpen())
    {
        //Evénements
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
        }
        
    }
    return EXIT_SUCCESS;

}

int main()
{
    vector<Motif> sequences;
    Motif m1 = {{'A', 'A', 'A', 'C', 'T', 'G', 'T'}, 0.1};
    sequences.push_back( m1) ;
    Motif m2 = {{'T', 'A', 'A', 'C', 'T','C','C'}, 0.13};
    sequences.push_back( m2) ;
    Motif m3 = {{'T', 'A', 'G', 'C' , 'C','A', 'C'}, 0.19};
    sequences.push_back( m3) ;
    Motif m4 = {{ 'A' , 'A', 'C', 'G', 'G','T', 'C'}, 0.187};
    sequences.push_back( m4) ;
    vector<vector<double>> test2 (Processing(sequences));
    afficher_logo(test2);
    return 0;
}
