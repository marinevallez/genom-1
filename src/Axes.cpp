

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
