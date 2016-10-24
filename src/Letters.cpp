//
//  Letters.cpp
//  partie 5
//
//  Created by Oriane Peter on 24.10.16.
//  Copyright © 2016 Oriane Peter. All rights reserved.
//

#include "Letters.hpp"

Letters::Letters()
{
    
}
void Letters::calculate_size( vector<double> tailles )
{
    for (size_t i(0); i < tailles_lettre.size(); ++i) {
        tailles_lettre[i] = tailles[i];
    };
}
array<double, 4> Letters::get_lettre_size()
{
    return tailles_lettre;
}
