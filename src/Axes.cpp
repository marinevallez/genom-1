//
//  Axes.cpp
//  partie 5
//
//  Created by Oriane Peter on 24.10.16.
//  Copyright © 2016 Oriane Peter. All rights reserved.
//

#include "Axes.hpp"

Axes::Axes()
{
    texture.loadFromFile("../res/graphPWM.png");
    sprite.setTexture(texture);
}
void Axes::draw(sf::RenderTarget& target)
{
    sprite.setOrigin(0, -50);
    
    target.draw(sprite);
}
