//
//  Rect.cpp
//  forceModel
//
//  Created by Ycl's Pro on 2021/6/12.
//

#include <stdio.h>
#include "Rect.h"

using namespace std;

// default constructor
Rect::Rect(){
    
}

Rect::Rect(int l, int t, int r, int b) {
    left = l;
    top = t;
    right = r;
    bottom = b;
}

void Rect::Assign(int l, int t, int r, int b) {
    left = l;
    top = t;
    right = r;
    bottom = b;
}
