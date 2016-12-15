

#ifndef Logo_hpp
#define Logo_hpp

#include <stdio.h>
#include <array>
#include <SFML/Graphics.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Letters.hpp"
#include "Axes.hpp"
#include "matrixprotein.hpp"

using namespace std;

class Logo {
    array<sf::Texture, 4> lettre_texture; // contain the texture of each letter filled with the emplacment of the image by the constructor
   
    
public:
    Logo(); // fills lettre_texture
    int afficher_logo(vector<Pattern> sequences, matrix); // Display the logo
    double size( double size_b); // calculates the size that must be displayed based on the calculation of the size depending on its probability (this function calculates a abritrary size and might result in  a lost of information !)
    sf::RectangleShape set_text( double taille, size_t lettre, double emplacement, int largeur); // for a given letter gives the appropriate size, texture and emplacment.
    vector<size_t> Letters_order(Letters Letters_, size_t j); // determine the order of probability for a given place on the sequence.
    
};

#endif /* Logo_hpp */
 