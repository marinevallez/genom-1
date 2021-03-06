

#ifndef Letters_hpp
#define Letters_hpp

#include <SFML/Graphics.hpp>
#include <vector>
#include <stdio.h>
#include "matrixprotein.hpp"


struct Motif
{
    std::vector<char> lettres;
    double binding_score;
};

using namespace std;

class Letters
{
public:
    void draw(sf::RenderTarget& target);
    Letters(matrix m); // the constructor process the size that should be given to a lettre given the probability of apperance at a certain position of the sequence.
    vector<vector<double>> get_lettre_size(); // returns  Letters_size
    sf::RectangleShape set_text( double taille, int lettre, double emplacement, int largeur); // determine the textrure the position and the size of a lettre
    double size( int ltr, Letters tentative);
    void set_lettre_size(vector<vector<double>> tl);
    
    private :
    sf::Texture texture;
    sf::Sprite sprite;
    vector<vector<double>> Letters_size;
    
    
};

#endif /* Letters_hpp */
