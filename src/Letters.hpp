//
//  Letters.hpp
//  partie 5
//
//  Created by Oriane Peter on 24.10.16.
//  Copyright © 2016 Oriane Peter. All rights reserved.
//

#ifndef Letters_hpp
#define Letters_hpp

#include <SFML/Graphics.hpp>
#include <vector>
#include <stdio.h>
#include <array>

using namespace std;

class Letters
{
public:
    void draw(sf::RenderTarget& target);
    Letters();
    void calculate_size(vector<double> tailles);
    array<double, 4> get_lettre_size();

private :
    sf::Texture texture;
    sf::Sprite sprite;
    array<double, 4> tailles_lettre; // 4 case qui représente les tailles de dans cette ordre ACGT
};

#endif /* Letters_hpp */
